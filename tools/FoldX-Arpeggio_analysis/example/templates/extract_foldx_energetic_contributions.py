#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob
import re
import pandas as pd
import argparse

# ============================================================
# ARGUMENT PARSER
# ============================================================
parser = argparse.ArgumentParser(description="Compute energetic contributions from FoldX outputs")
parser.add_argument("-i","--input", required=True, help="Main folder containing subfolders with .fxout files")
parser.add_argument("-m","--mutfile", required=True, help="Path to mutations.txt file")
parser.add_argument("-o","--outdir", required=True, help="Folder to store output CSV files")
parser.add_argument("--debug", action="store_true", help="Enable verbose debug output")

# A — optional block (positive_mean)
parser.add_argument("--enable-positive-mean", action="store_true",
                    help="Enable generation of positive_mean.csv (optional block A)")

args = parser.parse_args()

main_dir = args.input
mut_file_path = args.mutfile
out_dir = args.outdir
DEBUG = args.debug

# ============================================================
# DEBUG
# ============================================================
def debug(msg):
    print(f"[DEBUG] {msg}")

# ============================================================
# CHECK DIRECTORIES
# ============================================================
if not os.path.isdir(main_dir):
    raise ValueError(f"Input folder does not exist: {main_dir}")
if not os.path.isfile(mut_file_path):
    raise ValueError(f"Mutations file does not exist: {mut_file_path}")
os.makedirs(out_dir, exist_ok=True)

# ============================================================
# HELPERS
# ============================================================
def extract_key(filename):
    m = re.search(r"Repair_(\d+_\d+)", filename)
    return m.group(1) if m else None

def read_fxout(path):
    debug(f"Reading file: {path}")
    with open(path) as f:
        lines = f.readlines()

    header_idx = next((i for i, l in enumerate(lines) if l.startswith("Pdb")), None)
    if header_idx is None:
        raise ValueError(f"No header found in {path}")

    col_names = [c.replace(" ", "_") for c in lines[header_idx].strip().split("\t")]
    data_lines = lines[header_idx+1:]

    data = []
    for line in data_lines:
        if not line.strip():
            continue
        values = line.strip().split("\t")
        if len(values) < len(col_names):
            values += [None] * (len(col_names) - len(values))
        data.append(values)

    df = pd.DataFrame(data, columns=col_names)
    numeric_cols = df.columns[3:]
    df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors='coerce')

    return df, numeric_cols

def parse_delta(wt_file, mut_file):
    debug(f"Computing model delta: WT={wt_file}, MUT={mut_file}")
    wt, numeric_cols = read_fxout(wt_file)
    mut, _ = read_fxout(mut_file)
    delta_values = mut.iloc[0][numeric_cols].values - wt.iloc[0][numeric_cols].values
    return pd.DataFrame([delta_values], columns=numeric_cols)


# ============================================================
# REORDER COLUMNS: (wt_mean → wt_std → mut_mean → mut_std → delta)
# ============================================================
def reorder_columns(df):
    meta = ["Position", "Mutation"]
    
    base_cols = {}
    for col in df.columns:
        if col in meta:
            continue
        if "mean" in col:
            base = col.replace("_wt_mean", "")
            base_cols.setdefault(base, [])
        elif "std" in col:
            base = col.replace("_wt_std", "")
            base_cols.setdefault(base, [])
        else:
            base_cols.setdefault(col, [])

    ordered = meta[:]  # prima le colonne Position e Mutation
    for base in sorted(base_cols.keys()):
        for suffix in ["mean", "std", ""]:
            col = f"{base}{suffix}"
            if col in df.columns:
                ordered.append(col)

    return df[ordered]


# ============================================================
# READ MUTATIONS LIST
# ============================================================
with open(mut_file_path) as f:
    mutations_list = [line.strip() for line in f if line.strip()]

debug(f"Loaded {len(mutations_list)} mutations")

# ============================================================
# PER-MODEL DELTAS
# ============================================================
# ============================================================
# PER-MODEL DELTAS (WT ↔ MUT correctly paired)
# ============================================================
results = []

for pos_folder in sorted(os.listdir(main_dir)):
    pos_path = os.path.join(main_dir, pos_folder)
    if not os.path.isdir(pos_path):
        continue

    debug(f"Processing folder: {pos_folder}")

    # all files.fxout for this folder
    all_fx = sorted(glob.glob(os.path.join(pos_path, "Interaction*.fxout")))

    # Separate WT and MUT
    wt_files = [f for f in all_fx if "WT" in f]
    mut_files = [f for f in all_fx if "WT" not in f]

    if not wt_files or not mut_files:
        continue

    # create WT dictionary {model_key:file}
    wt_dict = {extract_key(f): f for f in wt_files}

    for idx, mutation in enumerate(mutations_list):
        mut_prefix = f"{idx+1}_"
        # all the mutants of the idx group
        mut_group_files = [f for f in mut_files if f"_Repair_{mut_prefix}" in f]
        if not mut_group_files:
            continue

        for mut_file in mut_group_files:
            # Extract key model from mutan, i.e '1_0', '1_1' ...
            mut_key = extract_key(mut_file)
            if mut_key in wt_dict:
                wt_file = wt_dict[mut_key]
            else:
                debug(f"No WT found for {mut_key}, using first WT as fallback")
                wt_file = wt_files[0]

            delta_df = parse_delta(wt_file, mut_file)
            delta_df["Position"] = pos_folder
            delta_df["Mutation"] = mutation
            delta_df["Model"] = mut_key
            results.append(delta_df)

# SAVE delta.csv
statistic_pattern_col = [
        "Interface_Residues",
        "Interface_Residues_Clashing",
        "Interface_Residues_VdW_Clashing",
        "Interface_Residues_BB_Clashing",
        "Number_of_Residues"
    ]
df = pd.concat(results)
statistic_columns = [c for c in df.columns if any(bad in c for bad in statistic_pattern_col)]
if statistic_columns:
    statistic_df = df[["Position", "Mutation"] + statistic_columns]
    dropping_col = [c for c in statistic_df.columns if "_std" in c]
    statistic_df_mean = statistic_df.drop(columns=dropping_col, errors="ignore")
    statistic_df_mean.to_csv(os.path.join(out_dir, "statistics.csv"),index=False)
    print("statistics.csv written.")
else:
    print("No statistic columns found; statistics.csv not written.")
df = df.drop(columns=statistic_columns, errors="ignore")
df = df.set_index(["Position", "Mutation", "Model"])
df.to_csv(os.path.join(out_dir, "delta.csv"))
print(f"delta.csv written.")

# ============================================================
# C + D — positives.csv (values> 0)
# ============================================================
pos_df = df.copy()
pos_df[pos_df <= 0] = pd.NA
pos_df.to_csv(os.path.join(out_dir, "positives.csv"))
print("positives.csv written.")

# ============================================================
# MEAN + STD across the 5 models (NO WT / MUT)
# ============================================================
df = df.reset_index()
energy_cols = [
    c for c in df.columns
    if c not in ["Position", "Mutation", "Model"]
]

agg_df = (
    df
    .groupby(["Position", "Mutation"], as_index=False)
    .agg(
        **{f"{c}_mean": (c, "mean") for c in energy_cols},
        **{f"{c}_std":  (c, lambda x: x.std(ddof=0)) for c in energy_cols},
    )
)

agg_df_ordered = reorder_columns(agg_df)
agg_df_ordered.to_csv(os.path.join(out_dir, "delta_mean.csv"), index=False)
print("delta_mean.csv written.")

# ============================================================
# POSITIVE MEAN ONLY (mean > 0)
# ============================================================

positive_df = agg_df.copy()

mean_cols = [c for c in positive_df.columns if c.endswith("_mean")]

for c in mean_cols:
    positive_df[c] = positive_df[c].where(positive_df[c] > 0, pd.NA)

std_cols = [c.replace("_mean", "_std") for c in mean_cols]
for c in std_cols:
    mean_col = c.replace("_std", "_mean")
    positive_df[c] = positive_df[c].where(positive_df[mean_col] > 0, pd.NA)
positive_df_ordered = reorder_columns(positive_df)
positive_df_ordered.to_csv(os.path.join(out_dir, "positive_mean.csv"))
print("positive_mean.csv written.")