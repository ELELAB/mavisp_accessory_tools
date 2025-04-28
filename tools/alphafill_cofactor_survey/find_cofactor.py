#!/usr/bin/env python3
#Script to query alphafill for a Uniprot entry to list potential cofactors with identity >=30 %

import argparse
import requests
import sys

# -----------------------------------------
# Parse command-line argument (UniProt AC)
# -----------------------------------------
parser = argparse.ArgumentParser(description="Find cofactors from AlphaFill JSON (identity ≥ 30%)")
parser.add_argument('-u', type=str, required=True, help='UniProt ID')
args = parser.parse_args()
upid = args.u

# -----------------------------------------
# Set up request to AlphaFill JSON endpoint
# -----------------------------------------
json_url = f"https://alphafill.eu/v1/aff/{upid}/json"
unique_cofactors = set()

try:
     # -------------------------------------
    # Request metadata from AlphaFill
    # -------------------------------------
    r = requests.get(json_url)
    r.raise_for_status()
except requests.RequestException as e:
    print(f"Error fetching data for {upid}: {e}")
    sys.exit(1)
try:
    data = r.json()
except ValueError as e:
    print(f"Error parsing JSON: {e}")
    sys.exit(1)
    # -------------------------------------
    # Loop over each structural hit (keep only >30% identity)
    # -------------------------------------
hits = data.get("hits", [])
if not hits:
    print(f"{upid},None")
    sys.exit(0)

for hit in hits:
    try:
        identity = float(hit["alignment"]["identity"])
        if identity < 0.3:
            continue

        transplants = hit.get("transplants", [])
        if transplants is None:
            continue

        for ligand in transplants:
            cofactor = ligand.get("analogue_id")
            if cofactor:
                unique_cofactors.add(cofactor)

    except (KeyError, ValueError, TypeError) as e:
        print(f"Error processing hit: {e}")
        continue
    # -------------------------------------
    # Output results
    # -------------------------------------
if unique_cofactors:
    cofactors_str = ' '.join(sorted(unique_cofactors))
    print(f"{upid},{cofactors_str}")
else:
    print(f"{upid},None")

