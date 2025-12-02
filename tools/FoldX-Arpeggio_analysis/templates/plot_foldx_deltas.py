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
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.input)
    cols = df.columns.tolist()
    cols.insert(0, cols.pop(cols.index('Interaction_Energy')))
    df = df[cols]

    delta_cols =  [c for c in df.columns if not c.endswith("_mean") and not c.endswith("_std")]

    # Loop on every line
    for idx, row in df.iterrows():
        mutation_name = f"{row['Position']}_{row['Mutation']}"
        values = []
        errors = []
        features = []

        for col in delta_cols:
            wt_std_col = f"{col}_wt_std"
            mut_std_col = f"{col}_mut_std"
            if wt_std_col not in df.columns or mut_std_col not in df.columns:
                continue
            mean_val = float(row[col])
            if mean_val != 0:
                err_val = np.sqrt(row[wt_std_col]**2 + row[mut_std_col]**2)
                values.append(mean_val)
                errors.append(err_val)
                features.append(col)

        # Plot
        plt.figure(figsize=(14,6))
        x = np.arange(len(features))  # 0, 1, 2, ..., N-1
        colors = ["red"] + ["blue"] * (len(values)-1)
        plt.bar(x, values, yerr=errors, color=colors, capsize=3)
        plt.xticks(x, features, rotation=90)  # assign correct labels
        plt.ylabel("Value")
        plt.title(mutation_name)
        plt.tight_layout()

        # Save
        out_file = os.path.join(args.outdir, f"{mutation_name}.png")
        plt.savefig(out_file, dpi=300)
        plt.close()
        print(f"Plot saved: {out_file}")
