#!/usr/bin/env python3
import pandas as pd
import requests
import time


def map_pdb_to_uniprot(pdb_ids, poll_interval=1, poll_timeout=60):
    """
    Map a list of PDB IDs to UniProt primary accession(s) using UniProt REST API.

    Returns a dict mapping PDB ID (e.g. "1LP1") -> UniProt accession (e.g. "P38507").
    """
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

        # print the full status once for debugging
        if first_status:
            print("FULL STATUS RESPONSE:")
            print(status)
            first_status = False

        job_status = status.get("jobStatus")
        print(f"Job status field: {job_status}")

        # Case A: status contains results directly -> use them
        if "results" in status and status.get("results") is not None:
            print("Found 'results' directly in status response — using those.")
            results_json = status["results"]
            break

        # Case B: status provides a redirectURL where results can be fetched
        if "redirectURL" in status and status["redirectURL"]:
            redirect = status["redirectURL"]
            print(f"Found redirectURL in status: {redirect} — fetching.")
            # follow redirect (may return JSON or stream)
            redirect_resp = requests.get(redirect)
            redirect_resp.raise_for_status()
            try:
                redirected = redirect_resp.json()
            except Exception:
                # not JSON: print for debugging and try text fallback
                print("Redirect returned non-JSON content; raw text:")
                print(redirect_resp.text)
                redirected = None

            # If redirected JSON contains 'results' use that
            if isinstance(redirected, dict) and "results" in redirected:
                results_json = redirected["results"]
                break
            # Otherwise, attempt to fetch the recommended stream endpoint below
            break

        # Case C: explicit jobStatus FINISHED
        if job_status == "FINISHED":
            print("JobStatus == FINISHED — will fetch results from stream endpoint.")
            break
        if job_status == "FAILED":
            raise RuntimeError(f"UniProt mapping job {job_id} FAILED")

        # nothing useful yet -> wait/poll
        time.sleep(poll_interval)
        elapsed += poll_interval
        if elapsed >= poll_timeout:
            raise TimeoutError(f"Timeout waiting for UniProt job {job_id} (waited {elapsed}s)")

    # If we already extracted results_json from status/redirect, use it.
    if results_json is None:
        # Try recommended stream endpoint for final results
        result_url = f"https://rest.uniprot.org/idmapping/uniprotkb/results/stream?jobId={job_id}"
        print(f"Fetching results from stream endpoint: {result_url}")
        results_resp = requests.get(result_url)
        # If 404 or no content, raise so caller sees it
        results_resp.raise_for_status()
        try:
            results_wrap = results_resp.json()
        except Exception:
            # If stream is not JSON, print raw text for debugging and return empty mapping
            print("Could not decode JSON from results stream. Raw text:")
            print(results_resp.text[:2000])
            return {}
        # the results are usually under 'results' key
        results_json = results_wrap.get("results", results_wrap)

    # results_json should be an iterable of mappings
    mapping = {}
    for r in results_json:
        src = r.get("from")
        to_field = r.get("to")
        acc = None
        # to_field can be a string (accession) or a dict with 'primaryAccession' etc
        if isinstance(to_field, str):
            acc = to_field
        elif isinstance(to_field, dict):
            # prefer primaryAccession, then 'primaryAccession' typing, then uniProtkbId, then try 'id' or 'accession'
            acc = to_field.get("primaryAccession") or to_field.get("primaryAccession".lower()) \
                  or to_field.get("uniProtkbId") or to_field.get("id") or to_field.get("accession")
            # some responses nest accession under keys: try 'primaryAccession' or 'primaryAccession'
            if not acc:
                # If there is a nested 'primaryAccession' inside another dict field, try common names:
                for key in ("primaryAccession", "uniProtkbId", "id", "accession"):
                    if key in to_field:
                        acc = to_field[key]
                        break
        # final fallback: if nothing parseable, try stringifying 'to'
        if not acc:
            acc = str(to_field)
        if src:
            mapping[src.upper()] = acc  # PDB ids are case-insensitive; normalize to uppercase

    return mapping


input_file = "mega_train.csv"
df = pd.read_csv(input_file)

# Work on unique WT_name entries so we don't query the same PDB repeatedly
df_unique = df.drop_duplicates(subset="WT_name", keep="first").copy()

# Select PDB-like entries (4-letter ID + ".pdb")
df_pdb = df_unique[
    df_unique["WT_name"].str.endswith(".pdb") & (df_unique["WT_name"].str.len() == 8)
].copy()

# Take only the first 5 PDB-like WT names
df_top5 = df_pdb.head(5).copy()

# Extract PDB ID (strip ".pdb")
df_top5.loc[:, "pdb_id"] = df_top5["WT_name"].str.replace(".pdb", "", regex=False)

# Get list of PDB IDs to query
pdb_ids = df_top5["pdb_id"].tolist()
print(f"Top 5 PDB IDs to map: {pdb_ids}")

# Query UniProt mapping API one-by-one (nice for debugging)
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

# Attach UniProt IDs to the top 5 entries
df_top5["uniprot_id"] = df_top5["pdb_id"].map(mapping)

# Keep only relevant columns
df_cleaned = df_top5[["WT_name", "uniprot_id", "pdb_id", "aa_seq", "name"]].rename(
    columns={"name": "protein_name"}
)

# Save cleaned dataset – this *does* contain UniProt accessions now
output_file = "tryout2_cleaned_top5.csv"
df_cleaned.to_csv(output_file, index=False)
print(f"Top 5 PDB-like entries saved to {output_file}")

