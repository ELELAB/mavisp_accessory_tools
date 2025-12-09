#!/usr/bin/env python3
import argparse

parser = argparse.ArgumentParser(
    description="Rename DNA residues (A,T,C,G → DA,DT,DC,DG) in specific chains"
)
parser.add_argument("-i", "--input", required=True, help="Input PDB file")
parser.add_argument("-o", "--output", required=True, help="Output PDB file")
parser.add_argument(
    "-c", "--chains", nargs="+", default=["P", "T"],
    help="Chains to rename (default: P T)"
)
args = parser.parse_args()

rename_map = {"A": "DA", "T": "DT", "C": "DC", "G": "DG"}

with open(args.input) as fin, open(args.output, "w") as fout:
    for line in fin:
        if line.startswith(("ATOM", "HETATM")):
            chain_id = line[21]
            resname = line[17:20].strip()
            if chain_id in args.chains and resname in rename_map:
                new_resname = rename_map[resname].ljust(3)
                line = line[:17] + new_resname + line[20:]
        fout.write(line)

print(f"Done. Saved renamed PDB to: {args.output}")

