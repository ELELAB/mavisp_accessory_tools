#!/usr/bin/env python3
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
        "index_csv_path",
        type=str,
        help="Path to the index.csv file from MAVISp."
    )
    parser.add_argument(
        "training_csv_path",
        type=str,
        help="Path to the training csv from ThermoMPNN to compare against."
    )
    args = parser.parse_args()

    # Load CSVs
    index_df = load_csv(args.index_csv_path)
    training_df = load_csv(args.training_csv_path)

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
