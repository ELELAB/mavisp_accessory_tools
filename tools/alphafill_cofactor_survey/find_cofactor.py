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
import os
import json
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser


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

def calculate_ligand_contacts(ligand, protein_residues, distance_cutoff):
    """ Calculate minimum distance and contacts for a ligand """
    min_distance = np.inf
    contacts = []

    for ligand_atom in ligand:
        if ligand_atom.element == 'H':
            continue
        for protein_res in protein_residues:
            contact_found = False
            for protein_atom in protein_res:
                if protein_atom.element == 'H':
                    continue
                distance = np.linalg.norm(ligand_atom.coord - protein_atom.coord)
                if distance < distance_cutoff:
                    contact_str = f"{protein_res.resname}{protein_res.id[1]}"
                    if contact_str not in contacts:
                        contacts.append(contact_str)
                        contact_found = True
                if distance < min_distance:
                    min_distance = distance
                if contact_found:
                    break
    return min_distance, contacts

def main():

    parser = argparse.ArgumentParser(description="Query AlphaFill for cofactors, with optional filtering.")
    parser.add_argument("-u", help="Single UniProt ID")
    parser.add_argument("-i", help="Index file with dataset, should have column 'UniProt AC'")
    parser.add_argument("-f", "--filter", help="Optional file with valid cofactor IDs")
    parser.add_argument("-o", "--output", default="summary_output", help="Output file prefix")
    parser.add_argument("-p", "--pdb", help="PDB file to analyze for ligands")
    parser.add_argument("-s", "--start", type=int, help="Start residue of trimmed region")
    parser.add_argument("-e", "--end", type=int, help="End residue of trimmed region")
    parser.add_argument("-c", "--cofactors", help="JSON file with cofactor definitions")
    parser.add_argument("-n", "--chain", help="Chain ID to analyze in PDB file")
    parser.add_argument("-d", "--distance", type=float, default=4.0, help="Distance cutoff for contacts in Angstroms (default: 4.0)")
    parser.add_argument("-l", "--ions", help="File with list of ion names (one per line)")


    args = parser.parse_args()

    cofactor_list = set()
    if args.cofactors:
        try:
            with open(args.cofactors, 'r') as f:
                cofactors_data = json.load(f)
            cofactor_list = set(cofactor.upper() for cofactor in cofactors_data.keys())
        except FileNotFoundError:
            print(f"Warning: {args.cofactors} not found, using basic classification")
            cofactor_list = set()
    if args.ions:
        ion_list = load_filter_file(args.ions)
    else:
        ion_list = set(['NA', 'K', 'CL', 'CA', 'MG', 'ZN', 'FE'])

    if not args.u and not args.i:
        print("Error: Provide either -u or -i")
        sys.exit(1)

    if args.u and args.i:
        print("Error: Provide either -u or -i, not both.")
        sys.exit(1)

    pdb_results = None

    if args.pdb:
        if not os.path.exists(args.pdb):
            print(f"Error: PDB file {args.pdb} not found")
            sys.exit(1)
        if args.start is None or args.end is None:
            print("Error: -s and -e required when using -p")
            sys.exit(1)
        if args.chain is None:
            print("Error: -n (chain ID) required when using -p")
            sys.exit(1)

        print(f"Analyzing PDB file: {args.pdb}")
        print(f"Chain: {args.chain}")
        print(f"Residue range: {args.start}-{args.end}")
        print(f"Distance cutoff: {args.distance} Å")

        pdb_parser = PDBParser(QUIET=True)
        structure =pdb_parser.get_structure('protein', args.pdb)
        model = structure[0]


        protein_residues = []
        ligands = []
        for chain in model:
            if chain.id != args.chain:
                continue
            for residue in chain:
                if residue.id[0] == ' ':
                    res_num = residue.id[1]
                    if args.start <= res_num <= args.end:
                        protein_residues.append(residue)
                elif residue.id[0] != 'W':
                    ligands.append(residue)
        pdb_results = []


        for ligand in ligands:
            min_distance, contacts = calculate_ligand_contacts(ligand, protein_residues, args.distance)

            if len(contacts) == 0:
                category = "no contacts"
            else:
                category = "inside trim"

            ligand_name = ligand.resname.upper()

            if ligand_name in cofactor_list:
                ligand_type = "cofactor"
            elif ligand_name in ion_list:
                ligand_type = "ion"
            else:
                ligand_type = "other ligand"

            pdb_results.append({
                'ligand_id': f"{ligand.resname}_{ligand.parent.id}_{ligand.id[1]}",
                'resname': ligand.resname,
                'chain': ligand.parent.id, 
                'residue_number': ligand.id[1],
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
    df_all = pd.DataFrame(all_rows, columns=["protein", "heteroatoms"])
    df_all.to_csv(out_all, index=False)

    if filter_set:
        df_filtered = pd.DataFrame(filtered_rows, columns=["protein", "cofactors"])
        df_filtered.to_csv(out_filtered, index=False)

    with open(out_unique, "w") as f:
        for c in sorted(unique_heteroatoms):
            f.write(c + "\n")

    if args.pdb and pdb_results is not None:
        ligand_output = f"{args.output}_ligand_contacts.csv"
        df_ligands = pd.DataFrame(pdb_results)
        df_ligands.to_csv(ligand_output, index=False)

        print(f"Ligand analysis complete: {len(pdb_results)} ligands found")
        print(f"Results saved to: {ligand_output}")

if __name__ == "__main__":
    main()
