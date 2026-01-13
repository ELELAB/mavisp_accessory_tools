#!/usr/bin/env python3

# Copyright (C) 2025 Eszter Toldi
# Technical University of Denmark, Danish Cancer Institute

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import pandas as pd
import requests
import time
import sys

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

    for attempt in range(3):
        try:
            response = requests.post(
                GRAPHQL_URL,
                json={"query": GRAPHQL_QUERY, "variables": variables},
                timeout=10
            )
            response.raise_for_status()
            break
        except Exception as e:
            if attempt == 2:
                raise RuntimeError(
                    f"[RCSB ERROR] {pdb_id}: API request failed after 3 attempts: {e}"
                )
            time.sleep(1)

    try:
        entry = response.json()["data"]["entry"]
        entities = entry["polymer_entities"]
    except KeyError as e:
        raise RuntimeError(f"{pdb_id}: unexpected GraphQL response, missing {e}")

    uniprots = []
    for entity in entities:
        container = entity.get("rcsb_polymer_entity_container_identifiers") or {}
        refs = container.get("reference_sequence_identifiers") or []

        for ref in refs:
            if ref.get("database_name") == "UniProt":
                uniprots.append(ref.get("database_accession"))

    if not uniprots:
        print(
            f"[RCSB INFO] {pdb_id}: No UniProt mappings found in RCSB. BLAST search will be attempted in step 2.")

    pdb_cache[pdb_id] = uniprots
    return uniprots

def main():
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
                    "wt_seq": row["wt_seq"],
                    "protein_name": row["name"]
                })
            else:
                for uid in uniprot_ids:
                    cleaned_rows.append({
                        "WT_name": wt_name,
                        "pdb_id": pdb_id,
                        "uniprot_id": uid,
                        "wt_seq": row["wt_seq"],
                        "protein_name": row["name"]
                    })
            # Sleep to avoid overloading API
            time.sleep(0.1)
        else:
            cleaned_rows.append({
                "WT_name": wt_name,
                "pdb_id": None,
                "uniprot_id": None,
                "wt_seq": row["wt_seq"],
                "protein_name": row["name"]
            })

    # Create DataFrame and save
    cleaned_df = pd.DataFrame(cleaned_rows)
    cleaned_df.reset_index(drop=True, inplace=True)
    cleaned_df.to_csv(OUTPUT_FILE, index=False)

    print(f"Full cleaned dataset saved to {OUTPUT_FILE}")
    print(f"Rows in output: {len(cleaned_df)}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(e)
        sys.exit(1)
