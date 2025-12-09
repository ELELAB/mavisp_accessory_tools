#!/usr/bin/env python3
import os
import re
import pandas as pd
from datetime import datetime

# CONFIG
all_root = "/data/raw_data/computational_data/mavisp_database_saturation/"  # directory with all MAVISp dataset folders
combined_csv_path = "./combined_mega_fireprot.csv"

# 1. Find newest ALL folder
date_pattern = re.compile(r"(\d{2})(\d{2})(\d{4})_")  # e.g. 01102025_ALL

candidates = []
for dirname in os.listdir(all_root):
    full_path = os.path.join(all_root, dirname)
    if os.path.isdir(full_path):
        m = date_pattern.match(dirname)
        if m:
            day, month, year = m.groups()
            date = datetime(int(year), int(month), int(day))
            candidates.append((date, full_path))

if not candidates:
    raise ValueError("No valid directories found in " + all_root)

newest_dir = max(candidates)[1]
print("Newest directory found:", newest_dir)

# 2. Find index.csv in newest dir
simple_dir = os.path.join(newest_dir, "simple_mode")
if not os.path.isdir(simple_dir):
    raise FileNotFoundError(f"No simple_mode directory found inside {newest_dir}")

index_csv_path = os.path.join(simple_dir, "index.csv")
if not os.path.exists(index_csv_path):
    raise FileNotFoundError(f"No index.csv found in {simple_dir}")

print("Using index file:", index_csv_path)

# 3. Load CSVs
index_df = pd.read_csv(index_csv_path)
combined_df = pd.read_csv(combined_csv_path)

# Standardize column names
combined_uniprot_col = "uniprot_id"
index_uniprot_col = "Uniprot AC" if "Uniprot AC" in index_df.columns else "uniprot_id"

# 4. Compare UniProt IDs
index_df["_in_combined"] = index_df[index_uniprot_col].isin(combined_df[combined_uniprot_col])

not_in_combined = index_df[~index_df["_in_combined"]].drop(columns="_in_combined")
in_combined = index_df[index_df["_in_combined"]].drop(columns="_in_combined")

# 5. Save results
not_in_combined.to_csv("benchmarking_dataset.csv", index=False)
in_combined.to_csv("overlap.csv", index=False)

print(f"{len(not_in_combined)} entries were not in combined ThermoMPNN training data. Saving to benchmarking_dataset.csv")
print(f"{len(in_combined)} entries overlap with ThermoMPNN training data. Saving to overlap.csv")

