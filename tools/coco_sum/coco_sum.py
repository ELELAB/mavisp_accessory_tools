#!/usr/bin/env python

# Copyright (C) 2025 Eirini  Giannakopoulou
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

"""
This script reads a CoCoNat prediction output TSV file and summarises the coiled-coil segments.
For each segment, it reports the residue range, predicted oligomeric state, and its probability.

"""

import sys
import argparse
import pandas as pd

def summarise_segments(input_file, output_file):
    # Read input file
    try:
        df = pd.read_csv(input_file, sep='\t')
    except FileNotFoundError:
        sys.exit(f"Error: Input file '{input_file}' not found.")

    # Check if file is empty
    if df.empty:
        sys.exit(f"Error: Input file '{input_file}' is empty.")

    # Check for expected columns
    expected_columns = {'ID', 'CC_CLASS', 'OligoState', 'POligo'}
    if not expected_columns.issubset(df.columns):
        missing_cols = expected_columns - set(df.columns)
        sys.exit(f"Error: Input file '{input_file}' is missing expected columns: {missing_cols}")

    summaries = []
    start_pos = None
    oligo_state = None
    poligo = None
    protein_id = None

    # Identify CC segments 
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

    # Handle segments at end of file
    if start_pos is not None:
        end_pos = len(df)
        summaries.append({
            'ID': protein_id,
            'ResRange': f"{start_pos}-{end_pos}",
            'OligoState': oligo_state,
            'POligo': poligo
     })

    # Check if output is empty
    if not summaries:
        sys.exit(f"Warning: No coiled-coil segments found in '{input_file}'. \
                 Output file will not be created.")

    # Write results to file
    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv(output_file, sep='\t', index=False)
    print(f"Summary written to: {output_file}")

def parse_arguments():
    parser = argparse.ArgumentParser(description="Summarise CoCoNat coiled-coil segments.")
    parser.add_argument('-i','--input_file', required=True, help='Path to CoCoNat output')
    parser.add_argument('-o','--output_file', required=True, help='Path to save the summary file')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    summarise_segments(args.input_file, args.output_file)
