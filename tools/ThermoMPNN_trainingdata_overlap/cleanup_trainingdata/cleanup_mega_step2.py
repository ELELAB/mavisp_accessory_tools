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
from Bio.Blast import NCBIWWW, NCBIXML
import time
import argparse

INPUT_FILE = "mega_cleaned_full.csv"
OUTPUT_FILE = "mega_cleaned_full_with_blast.csv"

def run_blast(seq, identity_threshold=30, coverage_threshold=70):
    """Run BLASTp against SwissProt and return UniProt IDs + identities."""
    for attempt in range(3):
        try:
            result_handle = NCBIWWW.qblast("blastp", "swissprot", seq, entrez_query="txid9606[ORGN]")
            blast_record = NCBIXML.read(result_handle)
            break
        except Exception as e:
            if attempt == 2:
                raise RuntimeError(f"[BLAST ERROR] BLAST failed after 3 attempts: {e}")
            time.sleep(1)

    hits = []
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            identity_percent = (hsp.identities / hsp.align_length) * 100
            coverage_percent = hsp.align_length / len(seq) * 100

            # Print results
            print(f"  → Hit: {alignment.hit_id}")
            print(f"    Description: {alignment.hit_def}")
            print(f"    Identity: {identity_percent:.2f}%")
            print(f"    Coverage: {coverage_percent:.2f}%")

            if identity_percent >= identity_threshold and coverage_percent >= coverage_threshold:
                # Extract UniProt accession
                hit_id = alignment.hit_id
                if "|" in hit_id:
                    parts = hit_id.split("|")
                    if len(parts) >= 2:
                        uniprot = parts[1]   # e.g. sp|P12345|NAME
                    else:
                        uniprot = hit_id
                else:
                    uniprot = hit_id

                hits.append({
                    "uniprot_id": uniprot,
                    "identity": identity_percent,
                    "coverage": coverage_percent
                })

    return hits

def main():
    parser = argparse.ArgumentParser(description="Run BLASTp on sequences and add UniProt hits.")
    parser.add_argument("--identity", type=float, default=30, help="Identity threshold (%)")
    parser.add_argument("--coverage", type=float, default=70, help="Coverage threshold (%)")
    args = parser.parse_args()

    print(f"\nUsing thresholds: Identity >= {args.identity}%, Coverage >= {args.coverage}%")

    # Load input
    df = pd.read_csv(INPUT_FILE)
    df_missing = df[df['uniprot_id'].isna()]
    df_known = df[~df['uniprot_id'].isna()]

    blast_rows = []
    for idx, row in df_missing.iterrows():
        wt_name = row['WT_name']
        seq = row['wt_seq']
        print(f"\n=== Running BLAST for {wt_name} ===")

        hits = run_blast(seq, identity_threshold=args.identity, coverage_threshold=args.coverage)

        if hits:
            for hit in hits:
                blast_rows.append({
                    "WT_name": wt_name,
                    "pdb_id": row['pdb_id'],
                    "uniprot_id": hit["uniprot_id"],
                    "blast_identity": hit["identity"],
                    "blast_coverage": hit["coverage"],
                    "wt_seq": seq,
                    "protein_name": row['protein_name']
                })
        else:
            blast_rows.append({
                "WT_name": wt_name,
                "pdb_id": row['pdb_id'],
                "uniprot_id": None,
                "blast_identity": None,
                "blast_coverage": None,
                "wt_seq": seq,
                "protein_name": row['protein_name']
            })

        time.sleep(3)

    # Combine known and BLASTed
    df_final = pd.concat([df_known, pd.DataFrame(blast_rows)], ignore_index=True)
    df_final.to_csv(OUTPUT_FILE, index=False)
    print(f"\nSaved updated dataset to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
