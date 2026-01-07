#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot FoldX mean energies per mutation with propagated uncertainty")
    parser.add_argument("-i", "--input", required=True, help="CSV file with FoldX energies")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    parser.add_argument("--xlabelsize", type=int, default=10, help="Font size for x-axis labels (default: 10)")
    parser.add_argument("--ylabelsize", type=int, default=12, help="Font size for y-axis labels (default: 12)")
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.input)

    # Detect the whole features (columna with_mean)
    features = [c.replace("_mean","") for c in df.columns if c.endswith("_mean")]

    features_kept = []

    for feat in features:
        mean_col = f"{feat}_mean"
        std_col  = f"{feat}_std"

        # if the column is not 0 is kept
        if mean_col in df.columns:
            if not np.isclose(df[mean_col].astype(float), 0.0).all():
                features_kept.append(feat)

    # rebuild the df with valid columns
    cols_to_keep = ["Position", "Mutation"]

    for feat in features_kept:
        cols_to_keep.append(f"{feat}_mean")
        if f"{feat}_std" in df.columns:
            cols_to_keep.append(f"{feat}_std")

    df = df[cols_to_keep]
    features = features_kept

    # Bring Interaction_Energy as first
    if "Interaction_Energy" in features:
        features = (
            ["Interaction_Energy"]
            + [f for f in features if f != "Interaction_Energy"]
        )

    for idx, row in df.iterrows():
        mutation_name = f"{row['Position']}_{row['Mutation']}"
        values = []
        errors = []

        for feat in features:
            mean_col = f"{feat}_mean"
            std_col  = f"{feat}_std"
            if mean_col in df.columns and std_col in df.columns:
                mean_val = row[mean_col]
                std_val  = row[std_col]
                values.append(mean_val)
                errors.append(std_val)

        # Plot
        plt.figure(figsize=(14,6))
        x = np.arange(len(features))
        colors = ["#FDBE85"] + ["#A6CEE3"] * (len(values) - 1)
        plt.bar(x, values, yerr=errors, color=colors, edgecolor='black', capsize=3)
        plt.xticks(x, features, rotation=90, fontsize=args.xlabelsize)
        plt.yticks(fontsize=args.ylabelsize)
        plt.ylabel("Value", fontsize=args.ylabelsize)
        plt.title(mutation_name)
        plt.tight_layout()

        # Save
        out_file = os.path.join(args.outdir, f"{mutation_name}.png")
        plt.savefig(out_file, dpi=300)
        plt.close()
        print(f"Plot saved: {out_file}")
