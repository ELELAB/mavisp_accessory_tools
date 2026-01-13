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

import os
import pandas as pd
import argparse

def load_csv(csv_path):
    """ Load csv file from user-specified path. """
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"File not found: {csv_path}")

    print(f"Using file:", csv_path)
    return pd.read_csv(csv_path)

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Compare index.csv with a training csv.")
    parser.add_argument(
        "-i", "--index",
        required=True,
        type=str,
        help="Path to the index.csv file from MAVISp."
    )

    parser.add_argument(
        "-t", "--training",
        required=True,
        type=str,
        help="Path to the training csv from ThermoMPNN."
    )
    args = parser.parse_args()

    # Load CSVs
    index_df = load_csv(args.index)
    training_df = load_csv(args.training)

    # Standardize column names
    combined_uniprot_col = "uniprot_id"
    index_uniprot_col = "Uniprot AC"

    # Compare UniProt IDs
    index_df["_in_training"] = index_df[index_uniprot_col].isin(training_df[combined_uniprot_col])
    not_in_training = index_df[~index_df["_in_training"]].drop(columns="_in_training")
    in_training = index_df[index_df["_in_training"]].drop(columns="_in_training")

    # Save results
    not_in_training.to_csv("benchmarking_dataset.csv", index=False)
    in_training.to_csv("overlap.csv", index=False)

    print(f"{len(not_in_training)} entries were not in training CSV. Saving to benchmarking_dataset.csv")
    print(f"{len(in_training)} entries overlap with training CSV. Saving to overlap.csv")

if __name__ == "__main__":
    main()
