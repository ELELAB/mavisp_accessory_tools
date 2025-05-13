#!/usr/bin/env python3

import argparse
import requests
import sys
import csv
import os

def fetch_cofactors(upid):
    """ Fetch cofactors from AlphaFill for one UniProt ID. Returns (status, list of cofactors or None) """
    json_url = f"https://alphafill.eu/v1/aff/{upid}/json"
    unique_cofactors = set()

    r = requests.get(json_url)
    r.raise_for_status()

    try:
        data = r.json()
    except ValueError as e:
        print(f"Error parsing JSON for {upid}: {e}")
        raise e

    hits = data.get("hits", [])
    if not hits:
        return "none", []

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
                    unique_cofactors.add(cofactor.upper())
        except (KeyError, ValueError, TypeError) as e:
            print(f"Error processing hit for {upid}: {e}")
            continue

    if unique_cofactors:
        return "ok", sorted(unique_cofactors)
    else:
        return "none", []

def load_filter_file(filepath):
    with open(filepath) as f:
        return set(line.strip().upper() for line in f if line.strip())

def main():
    parser = argparse.ArgumentParser(description="Query AlphaFill for cofactors, with optional filtering.")
    parser.add_argument("-u", help="Single UniProt ID")
    parser.add_argument("-i", help="File with UniProt IDs")
    parser.add_argument("-f", "--filter", help="Optional file with valid cofactor IDs")
    parser.add_argument("-o", "--output", default="summary_output", help="Output file prefix")

    args = parser.parse_args()

    if not args.u and not args.i:
        print("Error: Provide either -u or -i")
        sys.exit(1)

    if args.u and args.i:
        print("Error: Provide either -u or -i, not both.")
        sys.exit(1)

    upids = [args.u] if args.u else [line.strip() for line in open(args.i) if line.strip()]
    filter_set = load_filter_file(args.filter) if args.filter else None

    out_all = f"{args.output}.csv"
    out_filtered = f"{args.output}_filtered.csv" if filter_set else None
    out_unique = f"{args.output}_unique_heteroatoms.txt"

    all_rows = []
    filtered_rows = []
    unique_heteroatoms = set()

    for upid in upids:
        try:
            status, cofactors = fetch_cofactors(upid)
        except (requests.RequestException, ValueError) as e:
            print(f"Error fetching or parsing data for {upid}: {e}")
            all_rows.append([upid, "Error"])
            if filter_set:
                filtered_rows.append([upid, "Error"])
            continue

        if status == "none":
            all_rows.append([upid, "None"])
            if filter_set:
                filtered_rows.append([upid, "None"])
            continue

        cof_str = ' '.join(cofactors)
        all_rows.append([upid, cof_str])
        unique_heteroatoms.update(cofactors)

        if filter_set:
            matched = [c for c in cofactors if c in filter_set]
            filtered_rows.append([upid, ' '.join(matched) if matched else "None"])


    # Write output files
    with open(out_all, "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["protein", "heteroatoms"])
        writer.writerows(all_rows)

    if filter_set:
        with open(out_filtered, "w", newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["protein", "cofactors"])
            writer.writerows(filtered_rows)

    with open(out_unique, "w") as f:
        for c in sorted(unique_heteroatoms):
            f.write(c + "\n")

if __name__ == "__main__":
    main()


