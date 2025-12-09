#!/usr/bin/env python3
import pandas as pd
import requests
import time

INPUT_FILE = "mega_train.csv"
OUTPUT_FILE = "mega_cleaned_full.csv"
GRAPHQL_URL = "https://data.rcsb.org/graphql"

# GraphQL query to get UniProt for a PDB ID
GRAPHQL_QUERY = """
query structure($id: String!) {
  entry(entry_id: $id) {
    rcsb_id
    polymer_entities {
      rcsb_polymer_entity_container_identifiers {
        reference_sequence_identifiers {
          database_name
          database_accession
        }
      }
    }
  }
}
"""

# Cache for already queried PDBs
pdb_cache = {}

def get_uniprot_rcsb(pdb_id):
    """Query RCSB GraphQL API to get UniProt ID(s) for a given PDB ID."""
    if pdb_id in pdb_cache:
        return pdb_cache[pdb_id]

    variables = {"id": pdb_id}
    try:
        response = requests.post(GRAPHQL_URL, json={"query": GRAPHQL_QUERY, "variables": variables})
        response.raise_for_status()
        data = response.json().get("data", {}).get("entry", {})
        if not data:
            print(f"[RCSB ERROR] {pdb_id}: No data returned")
            pdb_cache[pdb_id] = []
            return []

        uniprots = []
        for entity in data.get("polymer_entities", []):
            refs = entity.get("rcsb_polymer_entity_container_identifiers", {}).get("reference_sequence_identifiers", [])
            for ref in refs:
                if ref.get("database_name") == "UniProt":
                    uniprots.append(ref.get("database_accession"))
        pdb_cache[pdb_id] = uniprots
        return uniprots
    except Exception as e:
        print(f"[RCSB ERROR] {pdb_id}: {e}")
        pdb_cache[pdb_id] = []
        return []

# Load dataset
df = pd.read_csv(INPUT_FILE)

# Keep only the first row per WT_name
df_unique = df.groupby("WT_name", as_index=False).first()

# Prepare cleaned dataset in long format
cleaned_rows = []
for idx, row in df_unique.iterrows():
    wt_name = row["WT_name"]
    pdb_like = wt_name.endswith(".pdb") and len(wt_name) == 8
    pdb_id = wt_name.replace(".pdb", "") if pdb_like else None

    if pdb_like:
        uniprot_ids = get_uniprot_rcsb(pdb_id)
        if not uniprot_ids:
            cleaned_rows.append({
                "WT_name": wt_name,
                "pdb_id": pdb_id,
                "uniprot_id": None,
                "aa_seq": row["aa_seq"],
                "protein_name": row["name"]
            })
        else:
            for uid in uniprot_ids:
                cleaned_rows.append({
                    "WT_name": wt_name,
                    "pdb_id": pdb_id,
                    "uniprot_id": uid,
                    "aa_seq": row["aa_seq"],
                    "protein_name": row["name"]
                })
        # Sleep to avoid overloading API
        time.sleep(0.1)
    else:
        cleaned_rows.append({
            "WT_name": wt_name,
            "pdb_id": None,
            "uniprot_id": None,
            "aa_seq": row["aa_seq"],
            "protein_name": row["name"]
        })

# Create DataFrame and save
cleaned_df = pd.DataFrame(cleaned_rows)
cleaned_df.reset_index(drop=True, inplace=True)
cleaned_df.to_csv(OUTPUT_FILE, index=False)

print(f"Full cleaned dataset saved to {OUTPUT_FILE}")
print(f"Rows in output: {len(cleaned_df)}")

