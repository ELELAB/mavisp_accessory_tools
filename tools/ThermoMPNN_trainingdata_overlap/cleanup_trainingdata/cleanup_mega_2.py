#!/usr/bin/env python3
import pandas as pd
import requests
import time


def map_pdb_to_uniprot(pdb_ids, poll_interval=1, poll_timeout=60):
    """
    Map a list of PDB IDs to UniProt primary accession(s) using UniProt REST API.

    Returns a dict mapping PDB ID (e.g. "1LP1") -> list of UniProt accessions (e.g. ["P38507", "P12345"]).
    """
    if not pdb_ids:
        return {}

    pdb_ids_str = ",".join(pdb_ids)
    run_url = "https://rest.uniprot.org/idmapping/run"
    data = {"from": "PDB", "to": "UniProtKB", "ids": pdb_ids_str}

    # submit mapping job
    response = requests.post(run_url, data=data)
    response.raise_for_status()
    job_id = response.json().get("jobId")
    if not job_id:
        raise RuntimeError("No jobId received from UniProt run endpoint")
    print(f"Submitted job {job_id} for PDB IDs: {pdb_ids}")

    status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"

    elapsed = 0
    first_status = True
    results_json = None

    while True:
        status_resp = requests.get(status_url)
        status_resp.raise_for_status()
        status = status_resp.json()

        if first_status:
            print("FULL STATUS RESPONSE:")
            print(status)
            first_status = False

        job_status = status.get("jobStatus")
        print(f"Job status field: {job_status}")

        # Case A: status contains results directly
        if "results" in status and status.get("results") is not None:
            results_json = status["results"]
            break

        # Case B: redirectURL
        if "redirectURL" in status and status["redirectURL"]:
            redirect = status["redirectURL"]
            redirect_resp = requests.get(redirect)
            redirect_resp.raise_for_status()
            try:
                redirected = redirect_resp.json()
            except Exception:
                redirected = None

            if isinstance(redirected, dict) and "results" in redirected:
                results_json = redirected["results"]
                break
            break

        # Case C: job finished
        if job_status == "FINISHED":
            break
        if job_status == "FAILED":
            raise RuntimeError(f"UniProt mapping job {job_id} FAILED")

        time.sleep(poll_interval)
        elapsed += poll_interval
        if elapsed >= poll_timeout:
            raise TimeoutError(f"Timeout waiting for UniProt job {job_id} (waited {elapsed}s)")

    # Fetch from stream endpoint if needed
    if results_json is None:
        result_url = f"https://rest.uniprot.org/idmapping/uniprotkb/results/stream?jobId={job_id}"
        results_resp = requests.get(result_url)
        results_resp.raise_for_status()
        try:
            results_wrap = results_resp.json()
        except Exception:
            print("Could not decode JSON from results stream. Raw text:")
            print(results_resp.text[:2000])
            return {}
        results_json = results_wrap.get("results", results_wrap)

    # build mapping dict: PDB -> list of UniProt accessions
    mapping = {}
    for r in results_json:
        src = r.get("from")
        to_field = r.get("to")
        acc_list = []

        if isinstance(to_field, str):
            acc_list.append(to_field)
        elif isinstance(to_field, dict):
            acc = (
                to_field.get("primaryAccession")
                or to_field.get("uniProtkbId")
                or to_field.get("id")
                or to_field.get("accession")
            )
            if acc:
                acc_list.append(acc)
        elif isinstance(to_field, list):
            for item in to_field:
                if isinstance(item, str):
                    acc_list.append(item)
                elif isinstance(item, dict):
                    acc = (
                        item.get("primaryAccession")
                        or item.get("uniProtkbId")
                        or item.get("id")
                        or item.get("accession")
                    )
                    if acc:
                        acc_list.append(acc)

        if src:
            if src.upper() in mapping:
                mapping[src.upper()].extend(acc_list)
            else:
                mapping[src.upper()] = acc_list

    # deduplicate and sort
    for k in mapping:
        mapping[k] = list(sorted(set(mapping[k])))

    return mapping


# ==== main script ====

input_file = "mega_train.csv"
df = pd.read_csv(input_file)

# Only unique entries by WT_name, keep first
df_unique = df.drop_duplicates(subset="WT_name", keep="first").copy()

# Extract PDB-like entries (4 chars + ".pdb" = 8 length and ending in .pdb)
df_pdb = df_unique[
    df_unique["WT_name"].str.endswith(".pdb") & (df_unique["WT_name"].str.len() == 8)
].copy()

# Add pdb_id column (4-letter code)
df_pdb.loc[:, "pdb_id"] = df_pdb["WT_name"].str.replace(".pdb", "", regex=False)

# Unique PDB IDs from this set
pdb_ids = sorted(df_pdb["pdb_id"].str.upper().unique().tolist())
print(f"Found {len(pdb_ids)} unique PDB-like IDs to map.")

# Call UniProt mapping
mapping = {}
if pdb_ids:
    try:
        mapping = map_pdb_to_uniprot(pdb_ids)
    except Exception as e:
        print(f"Error during UniProt mapping: {e}")
        mapping = {}

# Attach uniprot_id for the PDB entries (comma-separated if multiple)
df_pdb["uniprot_id"] = df_pdb["pdb_id"].str.upper().map(
    lambda x: ",".join(mapping.get(x, [])) if x in mapping else ""
)

# Merge this info back to the full unique table (left-join on WT_name)
df_merged = df_unique.merge(
    df_pdb[["WT_name", "pdb_id", "uniprot_id"]],
    on="WT_name",
    how="left"
)

# Rename 'name' to 'protein_name' if exists
if "name" in df_merged.columns:
    df_merged = df_merged.rename(columns={"name": "protein_name"})

# Keep selected columns
cols = [c for c in ["WT_name", "pdb_id", "uniprot_id", "aa_seq", "protein_name"] if c in df_merged.columns]
df_cleaned = df_merged[cols].copy()

output_file = "mega_train_unique_with_uniprot.csv"
df_cleaned.to_csv(output_file, index=False)
print(f"All unique WT entries saved to {output_file}")
print(f"Rows in output: {len(df_cleaned)}")
