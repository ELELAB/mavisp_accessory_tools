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
# MEAN / STD for each energetic contribution
# ============================================================
def compute_mean_delta(wt_files, mut_files):
    debug(f"Computing mean delta:\n WT files: {wt_files}\n MUT files: {mut_files}")

    wt_dfs, mut_dfs = [], []

    for f in wt_files:
        df, numeric_cols = read_fxout(f)
        wt_dfs.append(df[numeric_cols].iloc[0])

    for f in mut_files:
        df, _ = read_fxout(f)
        mut_dfs.append(df[numeric_cols].iloc[0])

    wt_df = pd.DataFrame(wt_dfs)
    mut_df = pd.DataFrame(mut_dfs)

    wt_mean = wt_df.mean()
    mut_mean = mut_df.mean()

    wt_std  = wt_df.std()
    mut_std = mut_df.std()

    out = {}

    for col in numeric_cols:
        out[f"{col}_wt_mean"] = wt_mean[col]
        out[f"{col}_wt_std"]  = wt_std[col]
        out[f"{col}_mut_mean"] = mut_mean[col]
        out[f"{col}_mut_std"]  = mut_std[col]
        out[f"{col}"] = mut_mean[col] - wt_mean[col]  # delta_mean

    return pd.DataFrame([out])

# ============================================================
# REORDER COLUMNS: (wt_mean → wt_std → mut_mean → mut_std → delta)
# ============================================================
def reorder_columns(df):
    base_cols = {}
    for col in df.columns:
        if "_wt_mean" in col:
            base = col.replace("_wt_mean", "")
            base_cols.setdefault(base, [])
        elif "_wt_std" in col:
            base = col.replace("_wt_std", "")
            base_cols.setdefault(base, [])
        elif "_mut_mean" in col:
            base = col.replace("_mut_mean", "")
            base_cols.setdefault(base, [])
        elif "_mut_std" in col:
            base = col.replace("_mut_std", "")
            base_cols.setdefault(base, [])
        else:
            # delta
            base_cols.setdefault(col, [])

    ordered = []
    for base in sorted(base_cols.keys()):
        for suffix in ["_wt_mean", "_wt_std", "_mut_mean", "_mut_std", ""]:
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
df = pd.concat(results).set_index(["Position", "Mutation", "Model"])
print(df)
statistic_columns = [c for c in df.columns if any(bad in c for bad in statistic_pattern_col)]
df = df.drop(columns=statistic_columns, errors="ignore")
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
# MEAN DELTA (with mean/std WT e MUT)
# ============================================================
mean_results = []

for pos_folder in sorted(os.listdir(main_dir)):
    pos_path = os.path.join(main_dir, pos_folder)
    if not os.path.isdir(pos_path):
        continue

    all_fx = sorted(glob.glob(os.path.join(pos_path, "Interaction*.fxout")))
    wt_files = [f for f in all_fx if "WT" in f]
    mut_files = [f for f in all_fx if "WT" not in f]

    if not wt_files or not mut_files:
        continue

    for idx, mutation in enumerate(mutations_list):
        mut_prefix = f"{idx+1}_"
        mut_group_files = [f for f in mut_files if f"_Repair_{mut_prefix}" in f]
        if not mut_group_files:
            continue

        wt_group_files = [f for f in wt_files if f"_Repair_{mut_prefix}" in f]
        if not wt_group_files:
            wt_group_files = wt_files

        mean_delta_df = compute_mean_delta(wt_group_files, mut_group_files)
        mean_delta_df["Position"] = pos_folder
        mean_delta_df["Mutation"] = mutation
        mean_results.append(mean_delta_df)

mean_df = pd.concat(mean_results).set_index(["Position", "Mutation"])
statistic_columns = [c for c in mean_df.columns if any(bad in c for bad in statistic_pattern_col)]
statistic_df_mean = mean_df[statistic_columns]
dropping_col = [c for c in mean_df.columns if "_std" in c]
statistic_df_mean = statistic_df_mean.drop(columns=dropping_col, errors="ignore")
statistic_df_mean.to_csv(os.path.join(out_dir, "statistics.csv"),index=False)
mean_df = mean_df.drop(columns=statistic_columns, errors="ignore")
mean_df = reorder_columns(mean_df)
mean_df.to_csv(os.path.join(out_dir, "delta_mean.csv"))
print("delta_mean.csv written.")

# ============================================================
# A — OPTIONAL positive_mean.csv
# ============================================================
if not args.enable_positive_mean:
    pos_mean = mean_df.copy()

contribs = sorted(set(c.rsplit("_", 2)[0] for c in pos_mean.columns if "_mut_mean" in c))
delta_dict = {}

contribs = sorted(c.replace("_mut_mean","") for c in pos_mean.columns if c.endswith("_mut_mean"))
delta_dict = {}

for c in contribs:
    wt_mean_col = f"{c}_wt_mean"
    mut_mean_col = f"{c}_mut_mean"
    wt_std_col  = f"{c}_wt_std"
    mut_std_col = f"{c}_mut_std"

    if wt_mean_col in pos_mean.columns and mut_mean_col in pos_mean.columns:
        delta = pos_mean[mut_mean_col] - pos_mean[wt_mean_col]
        # Put NaN where delta <= 0
        delta_filtered = delta.where(delta > 0, pd.NA)
        delta_dict[f"{c}_delta"] = delta_filtered
        # Associated std only if delta > 0
        wt_std_vals = pos_mean[wt_std_col].where(delta > 0, pd.NA)
        mut_std_vals = pos_mean[mut_std_col].where(delta > 0, pd.NA)
        delta_dict[f"{c}_wt_std"] = wt_std_vals
        delta_dict[f"{c}_mut_std"] = mut_std_vals

positive_delta_df = pd.DataFrame(delta_dict)
positive_delta_df.insert(0, "Mutation", pos_mean.index.get_level_values("Mutation"))
positive_delta_df.insert(0, "Position", pos_mean.index.get_level_values("Position"))
# write CSV
positive_delta_df.to_csv(os.path.join(out_dir, "positive_mean.csv"), index=False)
print(f"positive_mean.csv written: {os.path.join(out_dir, 'positive_mean.csv')}")
