#!/usr/bin/env python3
import pandas as pd

# Define column scheme
TRAINING_COLUMNS = [
    "pdb_id",
    "uniprot_id",
    "wt_seq",
    "protein_name",
    "source"
]

def clean_uniprot(x):
    """Remove UniProt version suffixes (e.g. P12345.2 → P12345)."""
    if pd.isna(x):
        return x
    return str(x).split(".")[0]

def prepare_mega(path):
    """Prepare MEGA dataset for training."""
    mega = pd.read_csv(path)
    mega["uniprot_id"] = mega["uniprot_id"].apply(clean_uniprot)
    mega["source"] = "mega"

    # Use only defined columns
    mega = mega[TRAINING_COLUMNS]

    return mega

def prepare_fireprot(path):
    """Prepare FireProt dataset for training."""
    fireprot = pd.read_csv(path)
    fireprot["uniprot_id"] = fireprot["uniprot_id"].apply(clean_uniprot)

    # Rename columns
    fireprot = fireprot.rename(columns={
        "pdb_id_corrected": "pdb_id",
        "sequence": "wt_seq"
    })

    fireprot["source"] = "fireprot"

    # Use only defined columns
    fireprot = fireprot[TRAINING_COLUMNS]

    return fireprot


def main():
    # Input files
    MEGA_INPUT = "mega_cleaned_full_with_blast.csv"
    FIREPROT_INPUT = "fireprot_clean.csv"

    # Output files
    MEGA_OUTPUT = "mega_cleaned_training.csv"
    FIREPROT_OUTPUT = "fireprot_cleaned_training.csv"
    COMBINED_OUTPUT = "combined_mega_fireprot.csv"

    # Prepare datasets
    mega = prepare_mega(MEGA_INPUT)
    fireprot = prepare_fireprot(FIREPROT_INPUT)

    # Save cleaned individual datasets
    mega.to_csv(MEGA_OUTPUT, index=False)
    fireprot.to_csv(FIREPROT_OUTPUT, index=False)

    # Combine datasets
    combined = pd.concat([mega, fireprot], ignore_index=True)

    # Remove duplicates by UniProt ID, take first occurence
    combined = combined.drop_duplicates(subset="uniprot_id", keep="first")

    combined.to_csv(COMBINED_OUTPUT, index=False)

    # Summary
    print(f"MEGA entries saved: {len(mega)} → {MEGA_OUTPUT}")
    print(f"FireProt entries saved: {len(fireprot)} → {FIREPROT_OUTPUT}")
    print(f"Combined entries saved: {len(combined)} → {COMBINED_OUTPUT}")


if __name__ == "__main__":
    main()
