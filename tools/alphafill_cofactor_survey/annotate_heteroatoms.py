import json
import argparse
import csv

# --- Arguments ---
parser = argparse.ArgumentParser(description="Annotate ligand IDs with cofactor names if these exist.")
parser.add_argument("-d", "--dictionary", required=True, help="Path to JSON dictionary file")
parser.add_argument("-c", "--compounds", required=True, help="Path to unique list of heteroatoms file (one per line)")
parser.add_argument("-o", "--output", default="all_heteroatoms_annotated.csv", help="Output CSV file")
args = parser.parse_args()

# --- Output files ---
annotated_output = "all_heteroatoms_annotated.csv"
cofactor_only_output = "cofactor_only.txt"


# --- Read dictionary ---
# from https://www.ebi.ac.uk/pdbe/api/doc/compounds.html
with open(args.dictionary, "r") as f:
    cofactor_dict = json.load(f)

# --- Create lookup: compound_id -> cofactor name ---
ligand_to_cofactor = {}
for name, entries in cofactor_dict.items():
    for entry in entries:
        for lid in entry.get("cofactors", []):
            if lid:
                ligand_to_cofactor[lid.upper()] = name


# --- read unique compound list from alphafill ---
with open(args.compounds, "r") as f:
    ligand_ids = [line.strip().upper() for line in f if line.strip()]

# --- Annotate compounds ---
annotated = []
cofactor_only = []

for ligand in ligand_ids:
    cofactor_name = ligand_to_cofactor.get(ligand)
    if cofactor_name:
        annotated.append((ligand, cofactor_name))
        cofactor_only.append(ligand)  # <-- tilføjet
    else:
        annotated.append((ligand, "not a cofactor"))


# --- Write  output ---
with open(annotated_output, "w", newline="") as out:
    writer = csv.writer(out)
    writer.writerow(["Heteroatom ID", "Annotation"])
    writer.writerows(annotated)

with open(cofactor_only_output, "w") as f:
    for ligand in cofactor_only:
        f.write(f"{ligand}\n")


print(f"Annotated {len(ligand_ids)} heteroatoms -> saved to '{annotated_output}'")
