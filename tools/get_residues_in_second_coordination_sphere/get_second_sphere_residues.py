import os
import re
import pandas as pd
import json
import argparse
from Bio import PDB
import yaml
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # Set to INFO or DEBUG to capture more detailed logs

# Create a console handler and set the level
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)

# Create a formatter and attach it to the handler
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)

# Add the handler to the logger
logger.addHandler(console_handler)

from arpeggio.core import InteractionComplex
from arpeggio.core.utils import max_mem_usage

def check_residues_in_pdb(pdb_file_path, residue_list):
    """
    Checks whether a list of specific residues is present in a PDB file, accounting for chain information.
    
    Parameters:
    ----------
    pdb_file_path : str
        The file path to the PDB file to be analyzed.
        
    residue_list : list of tuples
        A list containing residues to check, where each residue is specified by a tuple 
        in the format (chain_id, residue_name, residue_id). 
        - chain_id : str : Identifier for the chain the residue belongs to (e.g., "A", "B").
        - residue_name : str : The three-letter code of the residue (e.g., "ALA" for alanine).
        - residue_id : int : The residue number, as it appears in the PDB file.
    
    Returns:
    -------
    None
        The function prints a message indicating whether all specified residues are found 
        in the PDB file or lists the residues that are missing.
    """

    # Create a parser
    parser = PDB.PDBParser(QUIET=True)
    
    # Parse the structure
    structure = parser.get_structure("pdb_structure", pdb_file_path)
    
    # Extract residues with chain information from the structure
    pdb_residues = set()
    for model in structure:
        for chain in model:
            for residue in chain:
                # Get chain ID, residue name, and residue ID for unique identification
                chain_id = chain.id
                res_name = residue.get_resname()
                res_id = str(residue.id[1])  # Get the residue number (ignoring insertion code)
                
                # Add to the set of residues with chain ID
                pdb_residues.add((chain_id, res_name, res_id))


    # Check each residue in the list with chain information
    missing_residues = [res for res in residue_list if res not in pdb_residues]

    # Display results
    if missing_residues:
        raise ValueError("The following residues are not present in the PDB file:\n" +
                     "\n".join([f"Residue: {res[1]} with ID: {res[2]} in Chain: {res[0]}" for res in missing_residues]))

def from_pdb_to_cif(pdb,output_file):
    # Crea un parser PDB
    parser = PDB.PDBParser(QUIET=True)

    # Carica la struttura dal file PDB
    structure = parser.get_structure("MyStructure", pdb)

    # Scrive la struttura nel formato CIF
    io = PDB.MMCIFIO()
    io.set_structure(structure)
    io.save(output_file)

def run_arpeggio(cif_file,
                 residue,
                 write_hydrogenated,
                 minimise_hydrogens,
                 minimisation_steps,
                 minimisation_forcefield,
                 minimisation_method,
                 vdw_comp,
                 interacting_cutoff,
                 ph,
                 include_sequence_adjacent,
                 use_ambiguities,
                 output):

    i_complex = InteractionComplex(
        cif_file, vdw_comp, interacting_cutoff, ph
    )
    i_complex.structure_checks()

    if use_ambiguities:
        i_complex.address_ambiguities()

    if minimise_hydrogens:
        i_complex.minimize_hydrogens(
            minimisation_forcefield,
            minimisation_method,
            minimisation_steps,
        )

    if write_hydrogenated:
        dire = output.split("/")[0]
        file_name = output.split("/")[1].split(".")[0]
        i_complex.write_hydrogenated(dire, file_name)

    i_complex.initialize()
    i_complex.run_arpeggio(
         residue, interacting_cutoff, vdw_comp, include_sequence_adjacent
    )

    contacts = i_complex.get_contacts()

    with open(output, "w") as fp:
        json.dump(contacts, fp, indent=4, sort_keys=True)

    logger.info(f"Program End. Maximum memory usage was {max_mem_usage()}.")

    return contacts


three_one_letter_annotation = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                               'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                               'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                               'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}



parser = argparse.ArgumentParser(description = 'The get_second_sphere_residues.py script takes a PDB file '\
                                               'and a list of target residues as input, '
                                               'conducting an Arpeggio analysis for each specified residue. '\
                                               'This analysis identifies all interacting residues within a '\
                                               'user-defined cutoff distance, as specified in the configuration '\
                                               'file. Interactions are categorized by interaction type, and results '\
                                               'are outputted to a CSV file that details each target residue, '\
                                               'its interacting residues, the atoms involved in each interaction, '\
                                               'and the specific nature of each contact. The output file is further'\
                                               ' processed to remove clashing contacts and "proximal" contacts in case" '\
                                               '-k flag is not specified.')

parser.add_argument("-p","--pdb_file",
                       dest = "pdb_file",
                       required = True,
                       type = str,
                       help = "Path to the PDB file containing the structure for "\
                              "which interacting residues will be identified.")

parser.add_argument("-cr","--catalytic_residues",
                       dest="catalytic_residues",
                       required = True,
                       type = str,
                       help = "CSV file specifying residues of interest, "\
                              "with chain ID and residue position.")

parser.add_argument("-o","--output",
                       dest = "output_file",
                       required = False,
                       type = str,
                       default = "second_sphere_catalytic_residues.csv",
                       help = "Output file name. Default: "\
                              "second_sphere_catalytic_residues.csv")

parser.add_argument("-k","--keep_proximal",
                       action = "store_true",
                       required = False,
                       default = False,
                       help = "keep the 'proximal' contacts in the"\
                              " output file. Default (they are removed "\
                              "from the output).")


parser.add_argument("-y","--config_file",
                       dest = "config_file",
                       required = True,
                       type = str,
                       help = "Configuration file in yaml format defining "\
                              "Arpeggio parameters and the chem_comp filed "\
                              "to be added to the cif file.")


args = parser.parse_args()


# ---------- Read the csv and process the catalytic residue list ---------- #

catalytic_residues = pd.read_csv(args.catalytic_residues)
pdb_file_path = args.pdb_file

# Check the presence of PDB 
if not os.path.exists(pdb_file_path):
    print(f"The PDB file specified in the input {pdb_file_path} was not found. Exiting...")
    exit(1)

# Check the presence of config file
if not os.path.exists(args.config_file):
    print(f"The config file specified in the input {args.config_file} was not found. Exiting...")
    exit(1)

with open(args.config_file, 'r') as file:
    config = yaml.safe_load(file)

chem_comp = config.get('chem_comp', None)
if chem_comp is None:
    print("Error: 'chem_comp' not found in the configuration file.")
    exit(1)

config.pop('chem_comp')

#############################################################################
#                                                                           #
#                            ERROR HANDLING                                 #
#                                                                           #
#############################################################################

# -------------------- check the config file format ----------------------- #

expected_params = {
    "write_hydrogenated": bool,
    "minimise_hydrogens": bool,
    "minimisation_steps": (int, lambda x: x > 0),
    "minimisation_forcefield": (str, lambda x: x in ["MMFF94", "UFF", "Ghemical"]),
    "minimisation_method": (str, lambda x: x in ["DistanceGeometry", "SteepestDescent", "ConjugateGradients"]),
    "vdw_comp": (float, lambda x: x >= 0),
    "interacting_cutoff": (float, lambda x: x > 0),
    "ph": (float, lambda x: 0 <= x <= 14),
    "include_sequence_adjacent": bool,
    "use_ambiguities": bool,
    "output": str
}

for param, rules in expected_params.items():
    # Check if the parameter is in the config
    if param not in config:
        print(f"Error: The parameter '{param}' is missing in the {args.config_file}.")
        exit(1)
    
    # Extract expected type and condition function, if any
    expected_type = rules if isinstance(rules, type) else rules[0]
    condition = None if isinstance(rules, type) else rules[1]
    
    # Type check
    if not isinstance(config[param], expected_type):
        print(f"Error: The parameter '{param}' must be of type {expected_type.__name__}.")
        exit(1)
    
    # Specific boolean check for exact matches
    if expected_type is bool and config[param] not in [True, False]:
        print(f"Error: The parameter '{param}' must be exactly True or False.")
        exit(1)
    
    # Condition check, if any
    if condition and not condition(config[param]):
        print(f"Error: The parameter '{param}' has an invalid value: {config[param]}")
        exit(1)

# ------------- check column name in catalytic residue list --------------- #

expected_columns = {"catalytic_residue", "chain_id"}

# Check if the DataFrame has exactly the expected columns
if not set(catalytic_residues.columns) == expected_columns:
    print(f"The name of the columns in {args.catalytic_residues} is in the "
          f"wrong format. The columns name must be: \n"\
          f"catalytic_residue \n"\
          f"chain_id")

# --- check the format of the information in the catalytic residue list --- #

# Define the patterns for each column
catalytic_residue_pattern = r'^[A-Z][a-z]{2}\d+$' 
chain_id_pattern = r'^[A-Z]$'

valid_catalytic_residue = catalytic_residues['catalytic_residue'].str.match(catalytic_residue_pattern).all()
valid_chain_id = catalytic_residues['chain_id'].str.match(chain_id_pattern).all()

if not valid_catalytic_residue or not valid_chain_id:
    print(f"The information in the {args.catalytic_residues} are in the"\
          f"wrong format. The correct format is the following:"\
          f"catalytic_residue: Thr187 ecc \n "\
          f"chain_id: A,B,C ecc")
    exit(1)

# ----------- Check that the catalytic residues are in the pdb ----------- #

catalytic_residues[['catalytic_residue_capital_letter', 'residue_position']] = catalytic_residues['catalytic_residue'].str.extract(r'([A-Za-z]+)(\d+)')
catalytic_residues['catalytic_residue_capital_letter'] = catalytic_residues['catalytic_residue_capital_letter'].str.upper()

chain_id_residues = set(zip(catalytic_residues['chain_id'],
                            catalytic_residues['catalytic_residue_capital_letter'], 
                            catalytic_residues['residue_position']))

try:
    check_residues_in_pdb(pdb_file_path,chain_id_residues)
except ValueError as e:
    print(f"Error: {e}. Please verify that the residues are correctly "\
          f"formatted and exist in the PDB file.")
    exit(1)

#############################################################################
#                                                                           #
#                               RUN ARPEGGIO                                #
#                                                                           #
#############################################################################

residues_second_sphere = []
columns = ['catalytic_residue',
           'catalytic_residue_atom',
           'second_sphere_residue',
           'second_sphere_residue_atom',
           'contact',
           'distance',
           'interacting_entities',
           'type']

# --------------------- Create cif file from pdb -------------------------- #

cif_file = pdb_file_path.split(".")[0] + ".cif"
pdb_file = pdb_file_path.split(".")[0]

from_pdb_to_cif(pdb_file_path,cif_file)


# ------------------------- Process cif file  ----------------------------- #

with open(cif_file, "r") as file:
        existing_content = file.read()

# Add chem_comp fields in the cif files
updated_content = chem_comp + "\n" + existing_content


# Read the content of the CIF file
with open(cif_file, "r") as file:
    lines = file.readlines()

# Find the first line starting with "loop_" and insert `chem_comp` before it
for index, line in enumerate(lines):
    if line.startswith("loop_"):
        lines.insert(index, chem_comp + "\n")
        break  # Stop after inserting once

# Write the updated content back into the file
with open(cif_file, "w") as file:
    file.writelines(lines)

# Create output folder for arpeggio outputs

config['cif_file'] = cif_file
output_dir = config['output']
os.makedirs(output_dir, exist_ok=True)

# ----------------- Run arpeggio for catalytic residues ------------------- #

for res,pos,ch_id in zip(catalytic_residues['catalytic_residue_capital_letter'].to_list(),
                         catalytic_residues['residue_position'].to_list(),
                         catalytic_residues['chain_id'].to_list()):

    output_file = f"{pdb_file}_{res}_{pos}.json"
    config['output'] = "/".join([output_dir,output_file])
    config['residue'] = [f"/{ch_id}/{pos}/"]
    print(f"Running Arpeggio for {res}{pos} residue...")

    contacts = run_arpeggio(**config)

#############################################################################
#                                                                           #
#                         PROCESS ARPEGGIO OUTPUT                           #
#                                                                           #
#############################################################################

    for element in contacts:
        if element["bgn"]["auth_seq_id"] != element["end"]["auth_seq_id"]:
            catalytic_residue = ""
            catalytic_residue_atom = ""
            second_sphere_residue = ""
            second_sphere_residue_atom = ""
            if element["bgn"]["auth_seq_id"] == int(pos):
                catalytic_residue = element["bgn"]["label_comp_id"] + str(element["bgn"]["auth_seq_id"])
                catalytic_residue_atom = element["bgn"]["auth_atom_id"]
                second_sphere_residue = element["end"]["label_comp_id"] + str(element["end"]["auth_seq_id"])
                second_sphere_residue_atom = element["end"]["auth_atom_id"]
            if element["end"]["auth_seq_id"] == int(pos):
                second_sphere_residue = element["bgn"]["label_comp_id"] + str(element["bgn"]["auth_seq_id"])
                second_sphere_residue_atom= element["bgn"]["auth_atom_id"]
                catalytic_residue = element["end"]["label_comp_id"] + str(element["end"]["auth_seq_id"])
                catalytic_residue_atom = element["end"]["auth_atom_id"]

            if catalytic_residue and catalytic_residue_atom and second_sphere_residue and second_sphere_residue_atom:
                contact = "_".join(element['contact'])
                distance = element["distance"]
                interacting_entities = element["interacting_entities"]
                inter_type = element['type']

                interaction_information = [catalytic_residue,
                                           catalytic_residue_atom,
                                           second_sphere_residue,
                                           second_sphere_residue_atom,
                                           contact,
                                           distance,
                                           interacting_entities,
                                           inter_type]

                residues_second_sphere.append(interaction_information)

#############################################################################
#                                                                           #
#                               WRITE OUTPUT                                #
#                                                                           #
#############################################################################

df = pd.DataFrame(residues_second_sphere,columns=columns)
df.sort_values(by=['catalytic_residue',
                   'second_sphere_residue',
                   'catalytic_residue_atom',
                   'second_sphere_residue_atom'],inplace=True)

# ------------------------- filter dataframe ----------------------------- #

# Clashing contacts dataframe

clashes_df = df[df['contact'].str.contains('clash', case=False, na=False)]
clashes_df.to_csv("clashes_contacts.csv", sep = ",", index = False)

# filtered dataframe without Clashing contacts

filtered_df = df[~df['contact'].str.contains('clash', case=False, na=False)]

if args.keep_proximal:
    # filtered dataframe with "proximal" contacts
    filtered_df.to_csv(args.output_file, sep=",", index = False)
else:
    # filtered dataframe without "proximal" contacts
    filtered_df = filtered_df[~(filtered_df['contact'] == "proximal")]
    filtered_df.to_csv(args.output_file, sep=",", index = False)

# ------------------------ write mutation lists -------------------------- #

residues = list(set(filtered_df['catalytic_residue'].to_list() + \
                    filtered_df['second_sphere_residue'].to_list()))
residue_pos_dict={}
with open("catalytic_and_second_sphere_residues.txt","w") as f:
    for residue in residues:
        match = re.search("([A-Z][A-Z][A-Z])([0-9]+)", residue)
        if match is not None and len(match.groups()) == 2:
            if match.group(1) in three_one_letter_annotation.keys():
                substring_list = list(match.groups())
                f.write(f"{three_one_letter_annotation[substring_list[0]]}{substring_list[1]}"+"\n")
                residue_pos_dict[substring_list[1]] = three_one_letter_annotation[substring_list[0]]

        else:
            print(f"the following interactor {residue} will not be listed in the"
                  f" output txt file since it's not an aminoacid")

mutation_list=["A","C","D","E","F","G","H","I","L","M","N","P",
               "Q","R","S","T","V","Y","W","K"]

with open("catalytic_and_second_sphere_residues_saturation_mutlist.txt","w") as f:
    for key,value in residue_pos_dict.items():
        tmp_mutation_list = mutation_list.copy()
        tmp_mutation_list.remove(value)
        for residue in tmp_mutation_list:
            mutation=value+str((int(key)+1))+residue
            f.write(mutation+"\n")













