#!/usr/bin/env python3
import pandas as pd
import requests
import time

# -------------------------------
# Functions
# -------------------------------

def map_pdb_to_uniprot(pdb_ids):
    """Map a list of PDB IDs to UniProt IDs using UniProt REST API."""
    pdb_ids_str = ",".join(pdb_ids)
    url = "https://rest.uniprot.org/idmapping/run"
    data = {"from": "PDB", "to": "UniProtKB", "ids": pdb_ids_str}
    
    # submit mapping job
    response = requests.post(url, data=data)
    response.raise_for_status()
    job_id = response.json()["jobId"]
    
    # poll until finished
    status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
    while True:
        status = requests.get(status_url).json()
        if status.get("jobStatus") == "FINISHED":
            break
        time.sleep(1)
    
    # retrieve results
    result_url = f"https://rest.uniprot.org/idmapping/uniprotkb/results/{job_id}"
    results = requests.get(result_url).json()
    
    # build mapping dict
    mapping = {r['from']: r['to'] for r in results.get('results', [])}
    return mapping

# -------------------------------
# Load your dataset
# -------------------------------

input_file = "trainingdata.csv"  # change to your file
df = pd.read_csv(input_file)

# -------------------------------
# Keep only first row per WT_name
# -------------------------------

df_unique = df.drop_duplicates(subset="WT_name", keep="first").copy()

# -------------------------------
# Detect PDB-like names
# -------------------------------

df_unique["pdb_like"] = df_unique["WT_name"].str.endswith(".pdb") & (df_unique["WT_name"].str.len() == 8)
df_unique["pdb_id"] = df_unique["WT_name"].str.replace(".pdb", "", regex=False)
pdb_ids = df_unique.loc[df_unique["pdb_like"], "pdb_id"].tolist()

# -------------------------------
# Map PDB IDs to UniProt IDs
# -------------------------------

if pdb_ids:
    mapping = map_pdb_to_uniprot(pdb_ids)
    df_unique["uniprot_id"] = df_unique["pdb_id"].map(mapping)
else:
    df_unique["uniprot_id"] = None

# -------------------------------
# Keep only relevant columns
# -------------------------------

df_cleaned = df_unique[["WT_name", "uniprot_id", "pdb_id", "aa_seq", "name"]].rename(columns={"name": "protein_name"})

# -------------------------------
# Save cleaned dataset
# -------------------------------

output_file = "trainingdata_cleaned.csv"
df_cleaned.to_csv(output_file, index=False)
print(f"Cleaned data saved to {output_file}")

