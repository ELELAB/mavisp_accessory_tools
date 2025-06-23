#!/usr/bin/env python3

"""
This script reads a CoCoNat prediction output TSV file and summarises the coiled-coil segments.
For each segment, it reports the residue range, predicted oligomeric state, and its probability.

Usage:
    python summarise_coconat_segments.py --input_file=path/to/input.tsv --output_file=path/to/output.tsv
"""

import pandas as pd
import argparse

def summarise_segments(input_file, output_file):
    df = pd.read_csv(input_file, sep='\t')

    summaries = []
    start_pos = None
    oligo_state = None
    poligo = None
    protein_id = None

    for idx, row in df.iterrows():
        if row['CC_CLASS'] != 'i':
            if start_pos is None:
                start_pos = idx + 1  # 1-based indexing
                oligo_state = row['OligoState']
                poligo = row['POligo']
                protein_id = row['ID']
        else:
            if start_pos is not None:
                end_pos = idx
                summaries.append({
                    'ID': protein_id,
                    'RES_RANGE': f"{start_pos}-{end_pos}",
                    'OligoState': oligo_state,
                    'POligo': poligo
                })
                start_pos = None

    if start_pos is not None:
        end_pos = len(df)
        summaries.append({
            'ID': protein_id,
            'RES_RANGE': f"{start_pos}-{end_pos}",
            'OligoState': oligo_state,
            'POligo': poligo
        })

    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv(output_file, sep='\t', index=False)
    print(f"Summary written to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarise CoCoNat coiled-coil segments.")
    parser.add_argument('--input_file', required=True, help='Path to CoCoNat TSV output')
    parser.add_argument('--output_file', required=True, help='Path to save the summary TSV')

    args = parser.parse_args()

    summarise_segments(args.input_file, args.output_file)
