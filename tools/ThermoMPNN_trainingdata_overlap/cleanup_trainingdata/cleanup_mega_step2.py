#!/usr/bin/env python3
import pandas as pd
from Bio.Blast import NCBIWWW, NCBIXML
import time

INPUT_FILE = "mega_cleaned_full.csv"
OUTPUT_FILE = "mega_cleaned_full_with_blast.csv"
IDENTITY_THRESHOLD = 30   # 30% identity threshold

def run_blast(seq, identity_threshold=IDENTITY_THRESHOLD):
    """Run BLASTp against SwissProt and return UniProt IDs + identities."""
    try:
        result_handle = NCBIWWW.qblast("blastp", "swissprot", seq)
        blast_record = NCBIXML.read(result_handle)
        hits = []

        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                identity_percent = (hsp.identities / hsp.align_length) * 100

                # Print details live
                print(f"  → Hit: {alignment.hit_id}")
                print(f"    Description: {alignment.hit_def}")
                print(f"    Identity: {identity_percent:.2f}%")

                if identity_percent >= identity_threshold:
                    # Extract UniProt accession safely
                    hit_id = alignment.hit_id
                    if "|" in hit_id:
                        parts = hit_id.split("|")
                        if len(parts) >= 2:
                            uniprot = parts[1]   # sp|P12345|NAME
                        else:
                            uniprot = hit_id
                    else:
                        uniprot = hit_id

                    hits.append({
                        "uniprot_id": uniprot,
                        "identity": identity_percent
                    })

        return hits

    except Exception as e:
        print(f"[BLAST ERROR] {e}")
        return []

# Load input
df = pd.read_csv(INPUT_FILE)

df_missing = df[df['uniprot_id'].isna()]
df_known = df[~df['uniprot_id'].isna()]

blast_rows = []
for idx, row in df_missing.iterrows():
    wt_name = row['WT_name']
    seq = row['aa_seq']
    print(f"\n=== Running BLAST for {wt_name} ===")

    hits = run_blast(seq)

    if hits:
        for hit in hits:
            blast_rows.append({
                "WT_name": wt_name,
                "pdb_id": row['pdb_id'],
                "uniprot_id": hit["uniprot_id"],
                "blast_identity": hit["identity"],
                "aa_seq": seq,
                "protein_name": row['protein_name']
            })
    else:
        blast_rows.append({
            "WT_name": wt_name,
            "pdb_id": row['pdb_id'],
            "uniprot_id": None,
            "blast_identity": None,
            "aa_seq": seq,
            "protein_name": row['protein_name']
        })

    time.sleep(3)  # be nice to NCBI

# Combine
df_final = pd.concat([df_known, pd.DataFrame(blast_rows)], ignore_index=True)
df_final.to_csv(OUTPUT_FILE, index=False)

print(f"\nSaved updated dataset to {OUTPUT_FILE}")

