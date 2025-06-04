import os
import sys
import argparse
import warnings
import pandas as pd
import logging as log
import MDAnalysis as mda
from Bio.PDB import PDBParser
from cath_alphaflow.models.domains import GeneralDomainID, ChoppingPdbResLabel
from cath_alphaflow.commands.measure_globularity import (
    calculate_normed_radius_of_gyration,calculate_packing_density)

# Basic logging configuration
log.basicConfig(level=log.INFO, \
                format='%(levelname)s - %(message)s')

def validate_residue_range(pdb_file, trimming_range):
    """Validates that the provided residue trimming range is 
    within the actual residue range present in the PDB structure.

    Args:
        pdb_file (str): Path to the input PDB
        trimming_range (str): Residue range in format "start-end" ("1-122")
    ----------
    Returns:
        start_res (int): starting residue of the structure
        end_res (int): final residue of the structure 
    """
    if not os.path.isfile(pdb_file):
        log.error(f"❌ PDB file '{pdb_file}' does not exist.")
        sys.exit(1)

    u = mda.Universe(pdb_file)

    if len(u.segments) != 1:
        log.error(f"❌ Structure contains {len(u.segments)} chains. Input single chain only.")
        sys.exit(1)

    res_ids = sorted(set(res.resid for res in u.residues))
    if trimming_range:
        try:
            start_res, end_res = map(int, trimming_range.split('-'))
        except ValueError:
            log.error(f"❌ Invalid format for trimming range '{trimming_range}'. Use format e.g. '1-122'.")
            sys.exit(1)

        res_ids = set(res.resid for res in u.residues)
        requested_ids = set(range(start_res, end_res + 1))
        missing_ids = requested_ids - res_ids

        if missing_ids:
            log.error(f"❌ Residues not found in PDB: {sorted(missing_ids)}")
            sys.exit(1)
    else:
        start_res = res_ids[0]
        end_res = res_ids[-1]
        trimming_range = f"{start_res}-{end_res}"
        log.info(f"No residue range provided. Using full PDB.")

    log.info(f"✅ Residue range {trimming_range} and structure are valid.")
    return start_res, end_res

def calculate_asphericity(pdb_file, start_res, end_res):
    """Calculates the asphericity of the structure using MDAnalysis.

    Asphericity is a measure of how non-spherical a structure is,
    with 0 being a perfect sphere.

    Args:
        pdb_file (str):Path to the input PDB
        start_res (int): Starting residue of the structure
        end_res (int): Final residue of the structure

    Returns:
        float: Asphericity value for the structure
    """
    u = mda.Universe(pdb_file)

    selection = f"resid {start_res}:{end_res}"
    atoms = u.select_atoms(selection)
    asphericity = atoms.asphericity(pbc=False)

    log.info(f"✅ Asphericity = {asphericity:.3f}")

    return asphericity

def calculate_globularity_metrics(pdb_file, start_res, end_res):
    """Calculates globularity metrics for a trimmed domain using CATH AlphaFlow.

    The metrics include:
    - Normalized radius of gyration (Rg)
    - Packing density (PD)

    Args:
        pdb_file (str): Path to the input PDB
        start_res (int): Starting residue of the structure
        end_res (int): Final residue of the structure
    Returns:
        tuple: A tuple (normed_rg, packing_density), both floats
    """
    chopping_str = f"{start_res}-{end_res}"
    chopping = ChoppingPdbResLabel.from_str(f"{chopping_str}")
    model_id = os.path.splitext(os.path.basename(pdb_file))[0]
    domain_id = GeneralDomainID(raw_id=model_id, chopping=chopping)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(model_id, pdb_file)

    normed_rg = calculate_normed_radius_of_gyration(domain_id, structure, volume_resolution=5)
    packing_density = calculate_packing_density(domain_id, structure, include_all_atoms=False, distance_cutoff=5)
    log.info(f"✅ Normalized Rg = {normed_rg:.3f}, Packing Density = {packing_density:.3f}")

    return normed_rg, packing_density

def process_results(model_id, asphericity, normed_rg, packing_density, output_file,  asphericity_threshold, rg_threshold, pd_threshold):
    """Processes results to finalise an output
    file in CSV format. Adding column asserting
    Pass/Fail of thresholds.

    Args:
        model_id (str): Identifier for the model
        asphericity (flot): Asphericity value
        normed_rg (float): Normalised rg value
        packing_density (float): Packing Density
        output_file (str): Path to output file 
        asphericity_threshold (float): Asphericity threshold 
        rg_threshold (float): Normalised Rg threshold
        pd_threshold (float): Packing density threshold
    """
    result = {
        "model_id": model_id,
        "asphericity": asphericity,
        "normalized_radius_gyration": normed_rg,
        "packing_density": packing_density,}

    df = pd.DataFrame([result])

    df['satisfies_thresholds'] = (
    (df['asphericity'] < asphericity_threshold).astype(int) +
    (df['normalized_radius_gyration'] < rg_threshold).astype(int) +
    (df['packing_density'] > pd_threshold).astype(int)) >= 2

    df.to_csv(output_file, index=False, float_format="%.3f")
    log.info(f"✅ Final results written to {output_file}")

    row = df.iloc[0]
    status = "✅ PASS" if row['satisfies_thresholds'] else "❌ FAIL"
    log.info(f"{status}: {model_id} | Asp = {row['asphericity']:.3f}, Rg = {row['normalized_radius_gyration']:.3f}, PD = {row['packing_density']:.3f}")

def main():
    warnings.filterwarnings("ignore", category=UserWarning)
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    parser = argparse.ArgumentParser(description="Run globularity metrics using CATH AlphaFlow and MDAnalysis.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p", "--pdb", required=True, help="Input trimmed PDB file")
    parser.add_argument("-r", "--range", help="Residue trimming range (e.g., 160-214)")
    parser.add_argument("-o", "--output", help="Output CSV file")
    parser.add_argument("-a", "--asphericity", type=float, default=0.1, help=f"Threshold for asphericity")
    parser.add_argument("-pd", "--packing", type=float, default=10.333, help=f"Threshold for packing density")
    parser.add_argument("-rg", "--radius", type=float, default=0.356, help=f"Threshold for normalised radius of gyration")
    args = parser.parse_args()

    model_id = os.path.splitext(os.path.basename(args.pdb))[0]

    output_file = args.output or f"{model_id}_globularity_recap.csv"

    start_res, end_res = validate_residue_range(args.pdb, args.range)
    dom_asphericity = calculate_asphericity(args.pdb, start_res, end_res)
    normed_rg, packing_density = calculate_globularity_metrics(args.pdb, start_res, end_res)
    process_results(model_id, dom_asphericity, normed_rg, packing_density, output_file, args.asphericity, args.radius, args.packing)

if __name__ == "__main__":
    main()
