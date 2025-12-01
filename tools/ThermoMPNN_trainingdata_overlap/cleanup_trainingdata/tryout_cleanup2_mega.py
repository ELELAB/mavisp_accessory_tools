#!/usr/bin/env python3
import pandas as pd
import requests
import time


def map_pdb_to_uniprot(pdb_ids):
    """Map a list of PDB IDs to UniProt IDs using UniProt REST API."""
    pdb_ids_str = ",".join(pdb_ids)
    url = "https://rest.uniprot.org/idmapping/run"
    data = {"from": "PDB", "to": "UniProtKB", "ids": pdb_ids_str}

    # submit mapping job
    response = requests.post(url, data=data)
    response.raise_for_status()
    job_id = response.json()["jobId"]
    print(f"Submitted job {job_id} for PDB IDs: {pdb_ids}")

    # poll until finished
    status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
    while True:
        status = requests.get(status_url).json()
        print(f"Job status: {status.get('jobStatus')}")
        if status.get("jobStatus") == "FINISHED":
            break
        time.sleep(1)  # wait 1 second before polling again

    # retrieve results
    result_url = f"https://rest.uniprot.org/idmapping/uniprotkb/results/{job_id}"
    results = requests.get(result_url).json()
    print(f"Received mapping results for {len(results.get('results', []))} entries")

    # build mapping dict
    mapping = {r['from']: r['to'] for r in results.get('results', [])}
    return mapping


input_file = "mega_train.csv"
df = pd.read_csv(input_file)

df_unique = df.drop_duplicates(subset="WT_name", keep="first").copy()

df_pdb = df_unique[df_unique["WT_name"].str.endswith(".pdb") & (df_unique["WT_name"].str.len() == 8)].copy()
df_top5 = df_pdb.head(5).copy()  # only first 5 PDB-like WT names

df_top5.loc[:, "pdb_id"] = df_top5["WT_name"].str.replace(".pdb", "", regex=False)
pdb_ids = df_top5["pdb_id"].tolist()
#Query uniprot mapping API
#Query UniProt mapping API one-by-one
mapping = {}
for pdb_id in pdb_ids:
    print(f"Looking up UniProt ID for PDB: {pdb_id} ...")
    try:
        single_mapping = map_pdb_to_uniprot([pdb_id])
        mapping.update(single_mapping)
        print(f"  {pdb_id} → {single_mapping.get(pdb_id)}")
    except Exception as e:
        print(f"  Error mapping {pdb_id}: {e}")
    time.sleep(1)  # be polite to the API

df_top5["uniprot_id"] = df_top5["pdb_id"].map(mapping)

# Keep only relevant columns
df_cleaned = df_top5[["WT_name", "uniprot_id", "pdb_id", "aa_seq", "name"]].rename(columns={"name": "protein_name"})

# Save cleaned dataset
output_file = "tryout2_cleaned_top5.csv"
df_cleaned.to_csv(output_file, index=False)
print(f"Top 5 PDB-like entries saved to {output_file}")

