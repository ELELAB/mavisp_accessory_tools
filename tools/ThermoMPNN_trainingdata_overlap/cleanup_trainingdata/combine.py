#!/usr/bin/env python3
import pandas as pd

# Clean UniProt IDs
def clean_uniprot(x):
    if pd.isna(x):
        return x
    return str(x).split(".")[0]   # remove version suffix

# Load cleaned and blasted/uniprot ID matched mega dataset
mega = pd.read_csv("mega_cleaned_full_with_blast_B.csv")
mega["uniprot_id"] = mega["uniprot_id"].apply(clean_uniprot)

# Load cleaned FireProt dataset
fireprot = pd.read_csv("fireprot_clean.csv")
fireprot["uniprot_id"] = fireprot["uniprot_id"].apply(clean_uniprot)

# Standardize column names to match
fireprot = fireprot.rename(columns={
    "pdb_id_corrected": "pdb_id",
    "sequence": "aa_seq"
})

# Keep only relevant columns
fireprot = fireprot[["pdb_id", "uniprot_id", "aa_seq", "protein_name"]]

# Add dataset origin
mega["source"] = "mega"
fireprot["source"] = "fireprot"

# Combine
combined = pd.concat([mega, fireprot], ignore_index=True)

# Remove duplicates by UniProt
combined = combined.drop_duplicates(subset="uniprot_id", keep="first")

# Save final file
combined.to_csv("combined_mega_fireprot.csv", index=False)

print(f"Combined dataset saved, total entries: {len(combined)}")

