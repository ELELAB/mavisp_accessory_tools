#!/usr/bin/env python3
import pandas as pd
import requests

# ---------- Configuration ----------
INPUT_FILE = "mega_train.csv"  # Replace with your actual dataset path
TOP_N = 5
GRAPHQL_URL = "https://data.rcsb.org/graphql"

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

# ---------- Function to get UniProt ID ----------
def get_uniprot_rcsb(pdb_id: str) -> str | None:
    variables = {"id": pdb_id}
    try:
        response = requests.post(GRAPHQL_URL, json={"query": GRAPHQL_QUERY, "variables": variables})
        response.raise_for_status()
        data = response.json()
        print(f"DEBUG {pdb_id}: {data}")  # <-- see what the API actually returns
        dbs = data.get("data", {}).get("entry", {}).get("database_2", [])
        for db in dbs:
            if db.get("database_id") == "UniProt":
                return db.get("pdbx_database_accession")
    except Exception as e:
        print(f"Failed to fetch UniProt ID for {pdb_id}: {e}")
    return None

# ---------- Load dataset ----------
df = pd.read_csv(INPUT_FILE)

# Filter to WT_names that look like PDB IDs
df_pdb_like = df[df["WT_name"].str.endswith(".pdb") & (df["WT_name"].str.len() == 8)]

# Drop duplicates so we only have one per WT_name
df_unique = df_pdb_like.drop_duplicates(subset="WT_name", keep="first")

# Take top N
df_test = df_unique.head(TOP_N)

# ---------- Process ----------
cleaned_rows = []

for _, row in df_test.iterrows():
    wt_name = row["WT_name"]
    pdb_id = wt_name.replace(".pdb", "")
    uniprot_id = get_uniprot_rcsb(pdb_id)
    
    cleaned_rows.append({
        "WT_name": wt_name,
        "uniprot_id": uniprot_id,
        "pdb_id": pdb_id,
        "aa_seq": row["aa_seq"],
        "protein_name": row["name"]
    })

# ---------- Output ----------
df_cleaned = pd.DataFrame(cleaned_rows)
print(df_cleaned)

# Optional: save to CSV
df_cleaned.to_csv("test_mega_cleaned_top5.csv", index=False)

