#!/usr/bin/env python3
import pandas as pd
import requests

INPUT_FILE = "mega_train.csv" 
OUTPUT_FILE = "mega_cleaned.csv"
GRAPHQL_URL = "https://data.rcsb.org/graphql"

# GraphQL query to get UniProt for a PDB ID
GRAPHQL_QUERY = """
query structure($id: String!) {
  entry(entry_id: $id) {
    rcsb_id
    database_2 {
      database_id
      pdbx_database_accession
    }
  }
}
"""

def get_uniprot_rcsb(pdb_id):
    """Query RCSB GraphQL API to get UniProt ID(s) for a given PDB ID."""
    variables = {"id": pdb_id}
    try:
        response = requests.post(GRAPHQL_URL, json={"query": GRAPHQL_QUERY, "variables": variables})
        response.raise_for_status()
        data = response.json()
        databases = data.get("data", {}).get("entry", {}).get("database_2", [])
        for db in databases:
            if db["database_id"] == "UniProt":
                return db["pdbx_database_accession"]
    except Exception as e:
        print(f"Error fetching UniProt for {pdb_id}: {e}")
    return None

# Load the dataset
df = pd.read_csv(INPUT_FILE)

# Keep only the first row per WT_name
df_unique = df.groupby("WT_name", as_index=False).first()

# Prepare columns for cleaned dataset
cleaned_rows = []
for _, row in df_unique.iterrows():
    wt_name = row["WT_name"]
    pdb_like = wt_name.endswith(".pdb") and len(wt_name) == 8
    pdb_id = wt_name.replace(".pdb", "") if pdb_like else None
    uniprot_id = get_uniprot_rcsb(pdb_id) if pdb_like else None
    cleaned_rows.append({
        "WT_name": wt_name,
        "pdb_id": pdb_id,
        "uniprot_id": uniprot_id,
        "aa_seq": row["aa_seq"],
        "protein_name": row["name"]
    })
cleaned_df = pd.DataFrame(cleaned_rows)

# Reset index and save
cleaned_df.reset_index(drop=True, inplace=True)
cleaned_df.to_csv(OUTPUT_FILE, index=False)

print(f"Cleaned dataset saved to {OUTPUT_FILE}")

