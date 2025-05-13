import json
import argparse
import csv
import pandas as pd
import requests

# --- Arguments ---
parser = argparse.ArgumentParser(description="Annotate ligand IDs with cofactor names if these exist.")
parser.add_argument("-d", "--dictionary", required=True, help="Path to JSON dictionary file")
parser.add_argument("-c", "--compounds", required=True, help="Path to unique list of heteroatoms file (one per line)")
parser.add_argument("-o", "--output", default="all_heteroatoms_annotated.csv", help="Output CSV file")
args = parser.parse_args()

# --- RCSB GraphQL API config for metadata ---
GRAPHQL_ENDPOINT = "https://data.rcsb.org/graphql"

GRAPHQL_QUERY = """
query molecule($id: String!) {
    chem_comp(comp_id: $id) {
        chem_comp {
            id
            name
            formula
            type
        }
        pdbx_reference_molecule {
            chem_comp_id
            compound_details
            description
            class
            name
            represent_as
            type
        }
        drugbank {
            drugbank_info {
                drugbank_id
                drug_categories
                synonyms
                name
                description
            }
        }
    }
}
"""

# --- Functions for metadata ---
def query_pdb(compound_id):
    response = requests.post(
        GRAPHQL_ENDPOINT,
        json={"query": GRAPHQL_QUERY, "variables": {"id": compound_id}}
    )
    if response.status_code == 200:
        return response.json()
    else:
        print(f"Failed to fetch data for {compound_id}: {response.status_code}")
        return None

def extract_info(compound_id, data):
    info = {"Compound ID": compound_id}
    chem_comp = data.get("data", {}).get("chem_comp", {})

    if chem_comp.get("chem_comp"):
        cc = chem_comp["chem_comp"]
        info.update({
            "chem_name": cc.get("name"),
            "chem_formula": cc.get("formula"),
            "chem_type": cc.get("type")
        })

    if chem_comp.get("pdbx_reference_molecule"):
        ref = chem_comp["pdbx_reference_molecule"]
        info.update({
            "compound_details": ref.get("compound_details"),
            "compound_description": ref.get("description"),
            "compound_class": ref.get("class"),
            "compound_name": ref.get("name"),
            "represent_as": ref.get("represent_as"),
            "compound_type": ref.get("type")
        })

    if chem_comp.get("drugbank") and chem_comp["drugbank"].get("drugbank_info"):
        db = chem_comp["drugbank"]["drugbank_info"]
        info.update({
            "drugbank_id": db.get("drugbank_id"),
            "drug_categories": "; ".join(db.get("drug_categories", [])) if db.get("drug_categories") else None,
            "drug_synonyms": "; ".join(db.get("synonyms", [])) if db.get("synonyms") else None,
            "drug_name": db.get("name"),
            "drug_description": db.get("description")
        })

    return info



# --- Output files ---

cofactor_only_output = "cofactor_only.txt"


# --- Read dictionary ---
# from https://www.ebi.ac.uk/pdbe/api/doc/compounds.html
with open(args.dictionary, "r") as f:
    cofactor_dict = json.load(f)

# --- Create lookup: compound_id -> cofactor name ---
ligand_to_cofactor = {}
all_cofactors_from_json = set()
for name, entries in cofactor_dict.items():
    for entry in entries:
        for lid in entry.get("cofactors", []):
            if lid:
                ligand_to_cofactor[lid.upper()] = name
                all_cofactors_from_json.add(lid.upper()) 

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
        cofactor_only.append(ligand)
    else:
        annotated.append((ligand, "not a cofactor"))

annotated_df = pd.DataFrame(annotated, columns=["Compound ID", "Annotation"])

expanded_records = []

# --- Add metadata to annotated list ---
for _, row in annotated_df.iterrows():
    compound_id = row["Compound ID"]
    print(f"Querying PDB for {compound_id}...")
    data = query_pdb(compound_id)
    if data:
        info = extract_info(compound_id, data)
        info["Annotation"] = row["Annotation"]
        expanded_records.append(info)
    else:
        print(f"[WARNING] Failed to retrieve data for {compound_id}. No PDB info added.")
        with open("query_errors.log", "a") as log:
            log.write(f"{compound_id}\n")
        expanded_records.append({
            "Compound ID": compound_id,
            "Annotation": row["Annotation"],
            "chem_name": None,
            "chem_formula": None,
            "chem_type": None,
            "compound_details": None,
            "compound_description": None,
            "compound_class": None,
            "compound_name": None,
            "represent_as": None,
            "compound_type": None,
            "drugbank_id": None,
            "drug_categories": None,
            "drug_synonyms": None,
            "drug_name": None,
            "drug_description": None
            })

# --- Write output ---
expanded_df = pd.DataFrame(expanded_records)
expanded_df.insert(1, 'URL', expanded_df['Compound ID'].apply(lambda x: f"https://www.rcsb.org/ligand/{x}"))
anno = expanded_df.pop("Annotation")
expanded_df.insert(1, "Annotation", anno)
expanded_df.to_csv(args.output, index=False)
print(f"Saved enriched annotation to {args.output}")

with open(cofactor_only_output, "w") as f:
    for lid in sorted(all_cofactors_from_json):
        f.write(f"{lid}\n")


print(f"Annotated {len(ligand_ids)} heteroatoms -> saved to '{args.output}'")
