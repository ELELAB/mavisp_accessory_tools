#!/usr/bin/env python3
import pandas as pd

# Load fireprot dataset
df = pd.read_csv("fireprot_train.csv")

columns_keep = ["pdb_id_corrected", "uniprot_id", "sequence", "protein_name"]
df_small = df[columns_keep]
df_unique = df_small.drop_duplicates().reset_index(drop=True)

# Save to a new file
df_unique.to_csv("fireprot_clean.csv", index=False)

print(f"Saved {len(df_unique)} unique entries to fireprot_clean.csv")
