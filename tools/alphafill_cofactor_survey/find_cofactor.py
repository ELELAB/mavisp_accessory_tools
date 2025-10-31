# find_cofactors.py
# Copyright (C) 2025 Matteo Tiberti & Laura Kappel,
# Danish Cancer Institute

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

#!/usr/bin/env python3

import argparse
import requests
import sys
import csv
import os
import pandas as pd

# ====== NEW IMPORTS ======
from Bio.PDB import PDBParser
import numpy as np

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
    parser.add_argument("-i", help="Index file with dataset, should have column 'UniProt AC'")
    parser.add_argument("-f", "--filter", help="Optional file with valid cofactor IDs")
    parser.add_argument("-o", "--output", default="summary_output", help="Output file prefix")

# ======== ADDING PDB ARGUMENTS ========
    parser.add_argument("--pdb", help="PDB file to analyze for ligands")
    parser.add_argument("--start", type=int, help="Start residue of trimmed region")
    parser.add_argument("--end", type=int, help="End residue of trimmed region")

    args = parser.parse_args()

    if not args.u and not args.i:
        print("Error: Provide either -u or -i")
        sys.exit(1)

    if args.u and args.i:
        print("Error: Provide either -u or -i, not both.")
        sys.exit(1)

# ======== PDB ANALYSIS ========
    if args.pdb:
        if not os.path.exists(args.pdb):
            print(f"Error: PDB file {args.pdb} not found")
            sys.exit(1)
        if args.start is None or args.end is None:
            print("Error: --start and --end required when using --pdb")
            sys.exit(1)

        print(f"Analyzing PDB file: {args.pdb}")
        print(f"Residue range: {args.start}-{args.end}")

# PDB reading - get structures
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', args.pdb)
        results = []

        model = structure[0]

        for chain in model:
            for residue in chain:
                if residue.id[0] != ' ' and residue.id[0] != 'W':
                    contacts = []
                    min_distance = 100.0

                    # Compare ligand atoms to all protein atoms
                    for ligand_atom in residue:
                        for chain2 in model:
                            for protein_res in chain2:
                                # Skip if it is the same
                                if protein_res == residue or protein_res.id[0] != ' ':
                                    continue
                                    
                                for protein_atom in protein_res:
                                    # Calculate distance
                                    distance = np.linalg.norm(ligand_atom.coord - protein_atom.coord)
                                    
                                    # Check if within 4Å cutoff
                                    if distance < 4.0:
                                        contact_str = f"{protein_res.resname}{protein_res.id[1]}"
                                        if contact_str not in contacts:  
                                            contacts.append(contact_str)
                                    
                                    # Track shortest distance
                                    if distance < min_distance:
                                        min_distance = distance
                    
                    # Categorize based on contact residues
                    inside_count = 0
                    for contact in contacts:
                        # Extract residue number
                        res_num = int(''.join(filter(str.isdigit, contact)))
                        if args.start <= res_num <= args.end:
                            inside_count += 1
                    
                    # Categorization
                    if inside_count == len(contacts) and inside_count > 0:
                        category = "inside trim"
                    elif inside_count > 0 and inside_count < len(contacts):
                        category = "cross boundary"
                    elif inside_count == 0 and len(contacts) > 0:
                        category = "outside trim"
                    else:
                        category = "no contacts"
                    
                    # Ligand types
                    ligand_name = residue.resname.upper()
                    if ligand_name in ['NA', 'K', 'CL', 'CA', 'MG', 'ZN', 'FE']:
                        ligand_type = "ion"
                    elif ligand_name in ['ATP', 'ADP', 'GMP', 'GDP', 'NAP', 'NAD', 'FAD']:
                        ligand_type = "cofactor"
                    else:
                        ligand_type = "other ligand"
                    
                    results.append({
                        'ligand_id': f"{residue.resname}_{chain.id}_{residue.id[1]}",
                        'resname': residue.resname,
                        'chain': chain.id, 
                        'residue_number': residue.id[1],
                        'contact_residues': ';'.join(contacts),
                        'contact_count': len(contacts),
                        'shortest_distance': round(min_distance, 2),
                        'location_category': category,
                        'biological_label': ligand_type
                    })

    upids =[args.u] if args.u else pd.read_csv(args.i)['Uniprot AC'].dropna().astype(str).tolist()
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

    # Write ligand results to CSV
    if args.pdb and 'results' in locals() and results:
        ligand_output = f"{args.output}_ligand_contacts.csv"
        with open(ligand_output, "w", newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["ligand_id", "resname", "chain", "residue_number", "contact_residues", 
                           "contact_count", "shortest_distance", "location_category", "biological_label"])
            for r in results:
                writer.writerow([r['ligand_id'], r['resname'], r['chain'], r['residue_number'], 
                               r['contact_residues'], r['contact_count'], r['shortest_distance'],
                               r['location_category'], r['biological_label']])

        print(f"Ligand analysis complete: {len(results)} ligands found")
        print(f"Results saved to: {ligand_output}")

if __name__ == "__main__":
    main()
