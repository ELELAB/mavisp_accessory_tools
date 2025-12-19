#!/usr/bin/env python3
import os
import re
import pandas as pd
import json
import argparse
import yaml
import logging
import sys
from Bio import PDB
from arpeggio.core import InteractionComplex
from arpeggio.core.utils import max_mem_usage

# ──────────────────────────────── AMINO ACID DICTIONARY ────────────────────────────────
AA_1TO3 = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"
}

# ──────────────────────────────── LOGGER ────────────────────────────────
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


# ──────────────────────────────── UTILITIES ────────────────────────────────

def from_pdb_to_cif(pdb_path, output_file):
    """Convert a PDB file to CIF format.

    Args:
        pdb_path (str): Path to the input PDB file.
        output_file (str): Path to the output CIF file.
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("MyStructure", pdb_path)
    io = PDB.MMCIFIO()
    io.set_structure(structure)
    io.save(output_file)


def convert_residue(res):
    """Convert a residue name from 1-letter to 3-letter format if necessary.

    Args:
        res (str): Residue name (1–3 letters).

    Returns:
        str: Residue name in 3-letter format.
    """
    res = str(res).strip().upper()
    if len(res) == 1 and res in AA_1TO3:
        return AA_1TO3[res]
    elif 1 <= len(res) <= 3:
        return res.upper()
    else:
        sys.exit(f"Invalid residue name '{res}'. Must be 1–3 letters.")


def check_input_files(file_information, mode):
    """
    Validate the structure and format of residue input files.

    Ensures that required columns are present depending on 'mode', validates formats,
    and returns a set of tuples (chain_id, residue_3letter, position_str).

    Args:
        file_information (pd.DataFrame): Input DataFrame containing residue information.
        mode (str): Either 'single_molecule' or 'pairs_of_molecules'.

    Returns:
        set: Set of tuples (chain_id, residue_name, residue_position_as_str).
    """
    position_pattern = r'^\d+$'
    chain_id_pattern = r'^[A-Za-z]$'
    residue_pattern = r'^[A-Za-z]{1,3}$'

    if mode == "single_molecule":
        required_cols = {"first_molecule", "first_molecule_position", "first_molecule_chain_id"}
    elif mode == "pairs_of_molecules":
        required_cols = {
            "first_molecule", "first_molecule_position", "first_molecule_chain_id",
            "second_molecule", "second_molecule_position", "second_molecule_chain_id"
        }
    else:
        sys.exit(f"Unknown mode '{mode}'. Expected 'single_molecule' or 'pairs_of_molecules'.")

    missing_cols = [col for col in required_cols if col not in file_information.columns]
    if missing_cols:
        sys.exit(f"Missing required columns for mode '{mode}': {', '.join(missing_cols)}")

    chain_id_residues = set()
    prefixes = ["first"]
    if mode == "pairs_of_molecules":
        prefixes.append("second")

    for prefix in prefixes:
        mol_col = f"{prefix}_molecule"
        pos_col = f"{prefix}_molecule_position"
        chain_col = f"{prefix}_molecule_chain_id"

        df_sub = file_information[[mol_col, pos_col, chain_col]].copy()

        # Validate column formats
        valid_position = df_sub[pos_col].astype(str).str.match(position_pattern).all()
        valid_chain = df_sub[chain_col].astype(str).str.match(chain_id_pattern).all()
        valid_residue = df_sub[mol_col].astype(str).str.match(residue_pattern).all()

        if not (valid_position and valid_chain and valid_residue):
            sys.exit(
                f"Invalid format detected in '{prefix}' residue set.\n"
                f"Expected columns: {mol_col},{pos_col},{chain_col}\n"
                f"Example lines: ASP,699,B\nGLN,700,B\n"
                f"Rules:\n"
                f"- molecule: 1–3 letters (e.g., THR, GLY, ATP)\n"
                f"- position: numeric only (e.g., 187)\n"
                f"- chain_id: single letter (e.g., A, B, C)"
            )

        for _, row in df_sub.iterrows():
            chain_id_residues.add(
                (str(row[chain_col]), convert_residue(row[mol_col]), str(row[pos_col]))
            )

    logger.info(f"Validated {len(chain_id_residues)} residues from input file ({mode}).")
    return chain_id_residues


def check_residues_in_pdb(pdb_file_path, residue_list):
    """Verify that all residues from the input list exist in the PDB structure.

    Args:
        pdb_file_path (str): Path to the PDB file.
        residue_list (iterable): Iterable of (chain_id, res_name, res_pos) tuples.
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("pdb_structure", pdb_file_path)
    pdb_residues = set()
    for model in structure:
        for chain in model:
            for residue in chain:
                chain_id = chain.id
                res_name = residue.get_resname()
                res_id = str(residue.id[1])
                pdb_residues.add((str(chain_id), str(res_name), str(res_id)))
    missing = [res for res in residue_list if res not in pdb_residues]
    if missing:
        raise ValueError("Residues not found:\n" + "\n".join([f"{r[1]}{r[2]} (Chain {r[0]})" for r in missing]))


def run_arpeggio(cif_file, residue, write_hydrogenated, minimise_hydrogens,
                 minimisation_steps, minimisation_forcefield, minimisation_method,
                 vdw_comp, interacting_cutoff, ph, include_sequence_adjacent,
                 use_ambiguities, output):
    """Run Arpeggio on a CIF file to compute atomic contacts.

    This wrapper prepares InteractionComplex, runs the analysis and writes the JSON output.

    Returns:
        list: contacts (raw Arpeggio records)
    """
    i_complex = InteractionComplex(cif_file, vdw_comp, interacting_cutoff, ph)
    i_complex.structure_checks()

    if use_ambiguities:
        i_complex.address_ambiguities()
    if minimise_hydrogens:
        i_complex.minimize_hydrogens(minimisation_forcefield, minimisation_method, minimisation_steps)
    if write_hydrogenated:
        dire = os.path.dirname(output)
        file_name = os.path.basename(output).split(".")[0]
        i_complex.write_hydrogenated(dire, file_name)

    i_complex.initialize()
    i_complex.run_arpeggio(residue, interacting_cutoff, vdw_comp, include_sequence_adjacent)
    contacts = i_complex.get_contacts()

    with open(output, "w") as fp:
        json.dump(contacts, fp, indent=4, sort_keys=True)

    logger.info(f"Arpeggio run completed. Max memory usage: {max_mem_usage()}")
    return contacts


def flatten_contacts_data(record):
    """Flatten a nested Arpeggio contact record into a single-level dictionary."""
    flat = {}
    for k, v in record.items():
        if isinstance(v, dict):
            for subk, subv in v.items():
                flat[f"{k}_{subk}"] = subv
        else:
            flat[k] = v
    return flat


def parsed_arpeggio_output_to_df(residues_contact_information, positions, config):
    """Convert raw Arpeggio JSON contacts into a structured DataFrame.

    Args:
        residues_contact_information (list/df): flattened contacts (list of dicts or DataFrame)
        positions (iterable): positions of interest (strings or ints)
        config (dict): configuration containing output column names under
                       config["output_columns_customization"]

    Returns:
        pd.DataFrame: parsed and filtered DataFrame with standardized column names.
    """
    # Ensure we have a DataFrame
    if isinstance(residues_contact_information, list):
        df_filt = pd.DataFrame(residues_contact_information)
    else:
        df_filt = residues_contact_information.copy()

    positions = set(str(p) for p in positions)

    # initialize configured columns
    for col in config["output_columns_customization"].values():
        df_filt[col] = ""

    # iterate rows and populate configured output fields
    for idx, row in df_filt.iterrows():
        bgn_pos = row.get("bgn_auth_seq_id")
        end_pos = row.get("end_auth_seq_id")

        # case: target in begin
        if bgn_pos is not None and str(bgn_pos) in positions:
            df_filt.at[idx, config["output_columns_customization"]["col_input_residue"]] = f"{row.get('bgn_label_comp_id')}"
            df_filt.at[idx, config["output_columns_customization"]["col_input_residue_pos"]] = bgn_pos
            df_filt.at[idx, config["output_columns_customization"]["col_input_residue_chain_id"]] = row.get("bgn_auth_asym_id")
            df_filt.at[idx, config["output_columns_customization"]["col_input_residue_atom"]] = row.get("bgn_auth_atom_id")

            df_filt.at[idx, config["output_columns_customization"]["col_residue_in_contact"]] = f"{row.get('end_label_comp_id')}"
            df_filt.at[idx, config["output_columns_customization"]["col_residue_in_contact_pos"]] = row.get("end_auth_seq_id")
            df_filt.at[idx, config["output_columns_customization"]["col_residue_in_contact_atom"]] = row.get("end_auth_atom_id")
            df_filt.at[idx, config["output_columns_customization"]["col_residue_in_contact_chain_id"]] = row.get("end_auth_asym_id")

        # case: target in end
        elif end_pos is not None and str(end_pos) in positions:
            df_filt.at[idx, config["output_columns_customization"]["col_input_residue"]] = f"{row.get('end_label_comp_id')}"
            df_filt.at[idx, config["output_columns_customization"]["col_input_residue_pos"]] = end_pos
            df_filt.at[idx, config["output_columns_customization"]["col_input_residue_chain_id"]] = row.get("end_auth_asym_id")
            df_filt.at[idx, config["output_columns_customization"]["col_input_residue_atom"]] = row.get("end_auth_atom_id")

            df_filt.at[idx, config["output_columns_customization"]["col_residue_in_contact"]] = f"{row.get('bgn_label_comp_id')}"
            df_filt.at[idx, config["output_columns_customization"]["col_residue_in_contact_pos"]] = row.get("bgn_auth_seq_id")
            df_filt.at[idx, config["output_columns_customization"]["col_residue_in_contact_atom"]] = row.get("bgn_auth_atom_id")
            df_filt.at[idx, config["output_columns_customization"]["col_residue_in_contact_chain_id"]] = row.get("bgn_auth_asym_id")

    # final columns order
    cols_final = [
        config["output_columns_customization"]["col_input_residue"],
        config["output_columns_customization"]["col_input_residue_pos"],
        config["output_columns_customization"]["col_input_residue_chain_id"],
        config["output_columns_customization"]["col_input_residue_atom"],
        config["output_columns_customization"]["col_residue_in_contact"],
        config["output_columns_customization"]["col_residue_in_contact_pos"],
        config["output_columns_customization"]["col_residue_in_contact_atom"],
        config["output_columns_customization"]["col_residue_in_contact_chain_id"],
        "contact", "distance", "interacting_entities", "type"
    ]

    df_out = df_filt[cols_final]
    df_out = df_out[
        df_out[config["output_columns_customization"]["col_input_residue"]].astype(bool) &
        df_out[config["output_columns_customization"]["col_residue_in_contact"]].astype(bool)
    ]
    df_out["contact"] = df_out["contact"].apply(lambda x: ", ".join(x) if isinstance(x, list) else x)
    return df_out


def contacts_quality_filter(df, output_file, keep_proximal=False, logger=None):
    """Filter and export Arpeggio contacts by removing clashes or optionally proximal contacts.

    Writes:
        - clashes_contacts.csv (all clash entries)
        - output_file (filtered contacts)
    """
    log = logger.info if logger else print

    log("Filtering out clashing contacts...")
    clashes_df = df[df['contact'].astype(str).str.contains('clash', case=False, na=False)]
    clashes_df.to_csv("clashes_contacts.csv", sep=",", index=False)

    filtered_df = df[~df['contact'].astype(str).str.contains('clash', case=False, na=False)]

    if not keep_proximal:
        log("Removing proximal contacts from output file.")
        filtered_df = filtered_df[~(filtered_df['contact'] == "proximal")]
    else:
        log("Keeping proximal contacts in output file.")

    filtered_df.to_csv(output_file, sep=",", index=False)
    log(f"Filtered contacts written to {output_file}")
    return filtered_df


def generate_mutation_lists(df, config, three_one_letter_annotation, logger=None):
    """Generate single-point and saturation mutation lists from filtered contacts.

    Produces:
      - interacting_entities_list.txt  (AA1 + position)
      - interacting_entities_saturation_mutlist.txt (saturation mutagenesis list)
    """
    # invert provided mapping (1-letter -> 3-letter) to 3-letter -> 1-letter
    three_to_one = {v: k for k, v in three_one_letter_annotation.items()}

    log = logger.info if logger else print
    log("Generating mutation lists...")
    df['first_residue_pos'] = df[config["output_columns_customization"]["col_input_residue"]]+df[config["output_columns_customization"]["col_input_residue_pos"]].astype(str)
    df['second_residue_pos'] = df[config["output_columns_customization"]["col_residue_in_contact"]]+df[config["output_columns_customization"]["col_residue_in_contact_pos"]].astype(str)

    residues = list(set(
        df["first_residue_pos"].dropna().astype(str).tolist() +
        df["second_residue_pos"].dropna().astype(str).tolist()
    ))
    residue_pos_dict = {}

    with open("interacting_entities_list.txt", "w") as f:
        for residue in residues:
            match = re.search(r"([A-Z]{3})([0-9]+)", residue)
            if match is not None and len(match.groups()) == 2:
                aa3, pos = match.groups()
                if aa3 in three_to_one:
                    aa1 = three_to_one[aa3]
                    f.write(f"{aa1}{pos}\n")
                    residue_pos_dict[pos] = aa1
                else:
                    log(f"Skipping residue with unknown 3-letter code: {aa3}{pos}")
            else:
                log(f"Skipping non-amino acid interactor '{residue}'.")

    mutation_list = ["A","C","D","E","F","G","H","I","L","M","N","P","Q","R","S","T","V","Y","W","K"]

    with open("interacting_entities_saturation_mutlist.txt", "w") as f:
        for pos, wt_aa in residue_pos_dict.items():
            tmp_mut_list = [m for m in mutation_list if m != wt_aa]
            for mut_aa in tmp_mut_list:
                f.write(f"{wt_aa}{pos}{mut_aa}\n")

    log("Mutation list generation completed successfully.")


# ──────────────────────────────── MAIN ────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Generalized Arpeggio Runner (proteins, ligands, DNA). "
                    "Compute contacts and generate mutation lists for selected residues or pairs."
    )
    parser.add_argument("-p", "--pdb_file", required=True, help="Path to input PDB file.")
    parser.add_argument("-y", "--config_file", required=True, help="YAML config file for Arpeggio parameters.")
    parser.add_argument("-o", "--output_file", default="arpeggio_output.csv",
                        help="Output CSV file for filtered contacts.")
    parser.add_argument("-i", "--interface", action="store_true", default=False,
                        help="Keep only inter-chain interactions (interface between different chain IDs).")
    # Create mutually exclusive group
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-sm", "--single_molecule_file",
                       help="CSV for 'single_molecule' mode (first_molecule,first_molecule_position,first_molecule_chain_id).")
    group.add_argument("-pm", "--pairs_of_molecules_file",
                       help="CSV for 'pairs_of_molecules' mode (first_...,second_...).")
    parser.add_argument("-k", "--keep_proximal", action="store_true", default=False,
                        help="Keep 'proximal' contacts in output (default: remove them).")
    args = parser.parse_args()

    # Basic checks
    if not os.path.exists(args.pdb_file):
        raise FileNotFoundError(f"PDB file not found: {args.pdb_file}")
    if not os.path.exists(args.config_file):
        raise FileNotFoundError(f"Config file not found: {args.config_file}")

    with open(args.config_file, "r") as f:
        config = yaml.safe_load(f)

    chem_comp = config["arpeggio_run_parameters"].pop("chem_comp", None)
    if chem_comp is None:
        raise ValueError("Missing 'chem_comp' in YAML config")

    # Convert PDB -> CIF
    cif_file = args.pdb_file.replace(".pdb", ".cif")
    pdb_basename = os.path.splitext(os.path.basename(args.pdb_file))[0]
    from_pdb_to_cif(args.pdb_file, cif_file)

    # Inject chem_comp into CIF (if needed)
    with open(cif_file, "r") as f:
        lines = f.readlines()
    for idx, line in enumerate(lines):
        if line.startswith("loop_"):
            lines.insert(idx, chem_comp + "\n")
            break
    with open(cif_file, "w") as f:
        f.writelines(lines)

    # Ensure output directory exists
    os.makedirs(config["arpeggio_run_parameters"]["output"], exist_ok=True)
    config["arpeggio_run_parameters"]["cif_file"] = cif_file
    output_dir = config["arpeggio_run_parameters"]['output']

    # Process input files and validate residues against PDB
    chain_id_residues = set()
    if args.single_molecule_file:
        single_molecule_df = pd.read_csv(args.single_molecule_file)
        chain_id_residues = check_input_files(single_molecule_df, "single_molecule")
    if args.pairs_of_molecules_file:
        pairs_of_molecules_df = pd.read_csv(args.pairs_of_molecules_file)
        chain_id_residues = check_input_files(pairs_of_molecules_df, "pairs_of_molecules")

    if args.single_molecule_file or args.pairs_of_molecules_file:
        try:
            check_residues_in_pdb(args.pdb_file, chain_id_residues)
        except ValueError as e:
            logger.error(f"Residue validation failed: {e}")
            sys.exit(1)

    # ---------------- Mode 1: run on whole structure (all interactions) ----------------
    if not args.single_molecule_file and not args.pairs_of_molecules_file:
        logger.info("Running in mode: ALL interactions")
        # rename mapping to standardize Arpeggio columns to your configured output columns
        rename_map = {
            "bgn_label_comp_id": config["output_columns_customization"]["col_input_residue"],
            "bgn_auth_seq_id": config["output_columns_customization"]["col_input_residue_pos"],
            "bgn_auth_asym_id": config["output_columns_customization"]["col_input_residue_chain_id"],
            "bgn_auth_atom_id": config["output_columns_customization"]["col_input_residue_atom"],
            "end_label_comp_id": config["output_columns_customization"]["col_residue_in_contact"],
            "end_auth_seq_id": config["output_columns_customization"]["col_residue_in_contact_pos"],
            "end_auth_atom_id": config["output_columns_customization"]["col_residue_in_contact_atom"],
            "end_auth_asym_id": config["output_columns_customization"]["col_residue_in_contact_chain_id"],
        }

        config["arpeggio_run_parameters"]["residue"] = []  # analyze entire structure
        config["arpeggio_run_parameters"]["output"] = os.path.join(output_dir, "all_contacts.json")
        contacts = run_arpeggio(**config["arpeggio_run_parameters"])
        flatten_contacts = [flatten_contacts_data(contact) for contact in contacts]
        contacts_df = pd.DataFrame(flatten_contacts)
        contacts_df = contacts_df.rename(columns=rename_map)
        contacts_df["contact"] = contacts_df["contact"].apply(lambda x: ", ".join(x) if isinstance(x, list) else x)

        if args.interface:
            # keep only inter-chain contacts
            contacts_df = contacts_df.loc[
                contacts_df[config["output_columns_customization"]["col_residue_in_contact_chain_id"]] !=
                contacts_df[config["output_columns_customization"]["col_input_residue_chain_id"]]
            ]
            if contacts_df.empty:
                logger.info("No contacts were detected between distinct molecular entities (different chain IDs). No output will be generated.")
                sys.exit(0)

        quality_contacts_filtered_df = contacts_quality_filter(contacts_df, args.output_file, args.keep_proximal, logger)
        generate_mutation_lists(quality_contacts_filtered_df, config, AA_1TO3, logger)

    # ---------------- Mode 2: specific residues (single_molecule_file) ----------------
    elif args.single_molecule_file:
        logger.info(f"Running for {len(single_molecule_df)} target residues...")
        residues_contact_information = []
        positions = []
        for ch_id, res, pos in chain_id_residues:
            positions.append(pos)
            output_file = f"{pdb_basename}_{res}_{pos}.json"
            config["arpeggio_run_parameters"]['output'] = os.path.join(output_dir, output_file)
            config["arpeggio_run_parameters"]['residue'] = [f"/{ch_id}/{pos}/"]
            logger.info(f"Running Arpeggio for {res}{pos}...")
            contacts = run_arpeggio(**config["arpeggio_run_parameters"])
            residues_contact_information.extend(contacts)

        flatten_contacts = [flatten_contacts_data(contact) for contact in residues_contact_information]
        contacts_df = pd.DataFrame(flatten_contacts)
        parsed_contact_df = parsed_arpeggio_output_to_df(contacts_df, positions, config)

        if args.interface:
            parsed_contact_df = parsed_contact_df.loc[
                parsed_contact_df[config["output_columns_customization"]["col_residue_in_contact_chain_id"]] !=
                parsed_contact_df[config["output_columns_customization"]["col_input_residue_chain_id"]]
            ]
            if parsed_contact_df.empty:
                logger.info("No intermolecular contacts between specified residues were detected. No output will be generated.")
                sys.exit(0)

        quality_contacts_filtered_df = contacts_quality_filter(parsed_contact_df, args.output_file, args.keep_proximal, logger)
        generate_mutation_lists(quality_contacts_filtered_df, config, AA_1TO3, logger)

    # ---------------- Mode 3: pairs of residues ----------------
    elif args.pairs_of_molecules_file:
        logger.info(f"Running pairwise mode for {len(pairs_of_molecules_df)} pairs...")
        all_contacts = []
        positions = []
        for _, row in pairs_of_molecules_df.iterrows():
            r1 = row["first_molecule"]
            p1 = str(row["first_molecule_position"])
            c1 = row["first_molecule_chain_id"]
            r2 = row["second_molecule"]
            p2 = str(row["second_molecule_position"])
            c2 = row["second_molecule_chain_id"]

            positions.extend([p1, p2])
            residue_sel = [f"/{c1}/{p1}/", f"/{c2}/{p2}/"]
            config["arpeggio_run_parameters"]["residue"] = residue_sel
            out_file = os.path.join(output_dir, f"{r1}{p1}_{r2}{p2}_contacts.json")
            config["arpeggio_run_parameters"]["output"] = out_file

            logger.info(f"Running Arpeggio for pair {r1}{p1}-{r2}{p2}...")
            contacts = run_arpeggio(**config["arpeggio_run_parameters"])
            for c in contacts:
                if (str(c["bgn"]["auth_seq_id"]) in [p1, p2]) and (str(c["end"]["auth_seq_id"]) in [p1, p2]):
                    all_contacts.append(c)

        if not all_contacts:
            logger.info(
                f"None of the molecule pairs specified in the input file {args.pairs_of_molecules_file} "
                f"show direct contacts between the listed first and second residues. No output will be generated."
            )
            sys.exit(0)
        else:
            flatten_contacts = [flatten_contacts_data(contact) for contact in all_contacts]
            contacts_df = pd.DataFrame(flatten_contacts)
            parsed_contact_df = parsed_arpeggio_output_to_df(contacts_df, positions, config)
            quality_contacts_filtered_df = contacts_quality_filter(parsed_contact_df, args.output_file, args.keep_proximal, logger)
            generate_mutation_lists(quality_contacts_filtered_df, config, AA_1TO3, logger)
            logger.info(f"Pairwise analysis completed. Results saved to {args.output_file}")

    logger.info("Analysis completed successfully.")

if __name__ == "__main__":
    main()
