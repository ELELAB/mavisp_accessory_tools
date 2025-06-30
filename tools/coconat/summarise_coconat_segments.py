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
    try:
        df = pd.read_csv(input_file, sep='\t')
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        return
    except pd.errors.ParserError:
        print(f"Error: Failed to parse the input file '{input_file}'. Check if it's a valid TSV.")
        return
    except Exception as e:
        print(f"Unexpected error while reading '{input_file}': {e}")
        return

    # Check if file is empty
    if df.empty:
        print(f"Error: Input file '{input_file}' is empty.")
        return

    # Check for expected columns
    expected_columns = {'ID', 'CC_CLASS', 'OligoState', 'POligo'}
    missing_cols = expected_columns - set(df.columns)
    if missing_cols:
        print(f"Error: Input file '{input_file}' is missing expected columns: {missing_cols}")
        return

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
                    'ResRange': f"{start_pos}-{end_pos}",
                    'OligoState': oligo_state,
                    'POligo': poligo
                })
                start_pos = None

    if start_pos is not None:
        end_pos = len(df)
        summaries.append({
            'ID': protein_id,
            'ResRange': f"{start_pos}-{end_pos}",
            'OligoState': oligo_state,
            'POligo': poligo
     })


    # Check if any segments were found
    if not summaries:
        print(f"Warning: No coiled-coil segments found in '{input_file}'. Output file will not be created.")
        return

    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv(output_file, sep='\t', index=False)
    print(f"Summary written to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarise CoCoNat coiled-coil segments.")
    parser.add_argument('-i','--input_file', required=True, help='Path to CoCoNat TSV output')
    parser.add_argument('-o','--output_file', required=True, help='Path to save the summary TSV')

    args = parser.parse_args()

    summarise_segments(args.input_file, args.output_file)
