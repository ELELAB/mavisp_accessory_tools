# pLDDT analysis: are disagreements between FoldX5 and FoldX5.1 enriched in structurally uncertain regions?
# each mutation is treated as a separate data point (mutations at the same position share the same pLDDT)
import os
import glob
import argparse
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt
import seaborn as sns
warnings.filterwarnings("ignore")
mpl.rcParams.update({
    "font.family":       "sans-serif",
    "font.size":         12,
    "axes.titlesize":    13,
    "axes.labelsize":    12,
    "xtick.labelsize":   11,
    "ytick.labelsize":   11,
    "legend.fontsize":   11,
    "axes.spines.top":   False,
    "axes.spines.right": False,
})
color_agree    = "#4A90A4"
color_disagree = "#B22222"
color_dist     = "#6A9CC4"
col_ddg   = "Stability (FoldX5, alphafold, kcal/mol)"
col_plddt = "AlphaFold2 model pLDDT score"
col_mut   = "Mutation"


def classify_ddg(ddg):
    # MAVISp classification thresholds applied to FoldX ddG only
    # stabilizing: ddG <= -3, destabilizing: ddG >= 3, neutral: -2 < ddG < 2, uncertain: everything else
    if ddg <= -3:
        return "stabilizing"
    elif ddg >= 3:
        return "destabilizing"
    elif -2 < ddg < 2:
        return "neutral"
    else:
        return "uncertain"


def load_all_proteins(foldx5_dir, foldx51_dir):
    # load all csvs from both folders, match by filename
    # classify each mutation independently
    # label each mutation as agree or disagree between the two versions
    foldx5_files  = {Path(f).stem: f for f in glob.glob(os.path.join(foldx5_dir,  "*.csv"))}
    foldx51_files = {Path(f).stem: f for f in glob.glob(os.path.join(foldx51_dir, "*.csv"))}
    common_proteins = sorted(set(foldx5_files.keys()) & set(foldx51_files.keys()))
    print(f"proteins in foldx5:   {len(foldx5_files)}")
    print(f"proteins in foldx5.1: {len(foldx51_files)}")
    print(f"common proteins:      {len(common_proteins)}")
    all_protein_data = []
    for stem in common_proteins:
        protein_name = stem.replace("-simple_mode", "")
        try:
            df_foldx5  = pd.read_csv(foldx5_files[stem])
            df_foldx51 = pd.read_csv(foldx51_files[stem])
            df_foldx5  = df_foldx5.rename(columns={col_ddg: "ddg_foldx5"})
            df_foldx51 = df_foldx51.rename(columns={col_ddg: "ddg_foldx51"})
            df_merged = pd.merge(
                df_foldx51[[col_mut, "ddg_foldx51", col_plddt]],
                df_foldx5[[col_mut, "ddg_foldx5"]],
                on=col_mut
            )
            if df_merged.empty:
                continue
            df_merged["class_foldx5"]  = df_merged["ddg_foldx5"].apply(classify_ddg)
            df_merged["class_foldx51"] = df_merged["ddg_foldx51"].apply(classify_ddg)
            df_merged["agreement"]     = df_merged.apply(
                lambda row: "agree" if row["class_foldx5"] == row["class_foldx51"] else "disagree",
                axis=1
            )
            df_merged["ddg_diff"]     = abs(df_merged["ddg_foldx51"] - df_merged["ddg_foldx5"])
            df_merged["protein_name"] = protein_name
            all_protein_data.append(df_merged)
        except KeyError as e:
            print(f"  skipping {protein_name}: missing column {e}")
            continue
        except Exception as e:
            raise RuntimeError(f"unexpected error processing {protein_name}: {e}")
    if not all_protein_data:
        return pd.DataFrame()
    return pd.concat(all_protein_data, ignore_index=True)


def main():
    parser = argparse.ArgumentParser(description="plddt analysis: disagreements between foldx5 and foldx5.1")
    parser.add_argument("-f", "--foldx5_dir",  required=True, help="path to foldx5 dataset_tables folder")
    parser.add_argument("-i", "--foldx51_dir", required=True, help="path to foldx5.1 dataset_tables folder")
    parser.add_argument("-o", "--output_dir",  default="./plddt_results", help="where to save results")
    args = parser.parse_args()

    print("plddt analysis: foldx5 vs foldx5.1 classification agreement\n")
    print("[1] loading data...")
    df_mutations = load_all_proteins(args.foldx5_dir, args.foldx51_dir)
    if df_mutations.empty:
        print("error: no data loaded. check folder paths and csv structure.")
        return

    print(f"\ntotal mutations: {len(df_mutations)}")
    print(f"agree:    {(df_mutations['agreement'] == 'agree').sum()}")
    print(f"disagree: {(df_mutations['agreement'] == 'disagree').sum()}")

    os.makedirs(args.output_dir, exist_ok=True)
    csv_out = os.path.join(args.output_dir, "mutation_agreement_results.csv")
    df_mutations.to_csv(csv_out, index=False)
    print(f"\nsaved results table: {csv_out}")

    print("\n[2] global plots")

    # kde + barplot of plddt scores for agree vs disagree mutations
    os.makedirs(args.output_dir, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    ax = axes[0]
    ax_bar = axes[1]
    for group, color in [("agree", color_agree), ("disagree", color_disagree)]:
        subset = df_mutations[df_mutations["agreement"] == group][col_plddt]
        if subset.empty:
            continue
        sns.kdeplot(subset, ax=ax, color=color,
                    label=f"{group} (n={len(subset)})",
                    fill=True, alpha=0.3, linewidth=2)
        ax.axvline(subset.median(), color=color, linestyle="--", linewidth=1.5, alpha=0.8)
    ax.set_xlabel("pLDDT score")
    ax.set_ylabel("density")
    ax.set_title("pLDDT distribution by classification agreement\nFoldX5 vs FoldX5.1 (all proteins)")
    ax.legend(frameon=False)
    bin_edges = np.arange(0, 101, 3)
    for group, color in [("agree", color_agree), ("disagree", color_disagree)]:
        subset = df_mutations[df_mutations["agreement"] == group][col_plddt]
        if subset.empty:
            continue
        counts, edges = np.histogram(subset, bins=bin_edges)
        counts = counts / counts.sum()
        ax_bar.bar(edges[:-1], counts, width=3, color=color, alpha=0.6,
                   label=f"{group} (n={len(subset)})", align="edge", edgecolor="white")
    ax_bar.set_xlabel("pLDDT score")
    ax_bar.set_ylabel("fraction of mutations")
    ax_bar.set_title("pLDDT distribution by classification agreement\n(3 Å bins)")
    ax_bar.legend(frameon=False)
    plt.tight_layout()
    out = os.path.join(args.output_dir, "plddt_distribution_agree_disagree.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"saved: {out}")

    # scatter: plddt vs absolute ddg difference per mutation
    fig, ax = plt.subplots(figsize=(8, 6))
    for group, color in [("agree", color_agree), ("disagree", color_disagree)]:
        subset = df_mutations[df_mutations["agreement"] == group]
        if subset.empty:
            continue
        ax.scatter(subset[col_plddt], subset["ddg_diff"],
                   color=color, alpha=0.5, s=40,
                   label=f"{group} (n={len(subset)})",
                   edgecolors="none")
    ax.set_xlabel("pLDDT score")
    ax.set_ylabel("|ΔΔG FoldX5.1 − ΔΔG FoldX5| (kcal/mol)")
    ax.set_title("pLDDT vs prediction difference per mutation\nFoldX5 vs FoldX5.1 (all proteins)")
    ax.legend(frameon=False)
    plt.tight_layout()
    out = os.path.join(args.output_dir, "plddt_vs_ddg_diff_scatter.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"saved: {out}")

    print("\n[3] per-protein summary")

    # per-protein: median plddt of disagree - median plddt of agree
    per_protein_results = []
    for protein_name, protein_group in df_mutations.groupby("protein_name"):
        agree_plddt    = protein_group[protein_group["agreement"] == "agree"][col_plddt]
        disagree_plddt = protein_group[protein_group["agreement"] == "disagree"][col_plddt]
        if agree_plddt.empty or disagree_plddt.empty:
            continue
        delta_median = disagree_plddt.median() - agree_plddt.median()
        per_protein_results.append({
            "protein_name":          protein_name,
            "median_plddt_agree":    agree_plddt.median(),
            "median_plddt_disagree": disagree_plddt.median(),
            "delta_median_plddt":    delta_median,
            "n_agree_mutations":     len(agree_plddt),
            "n_disagree_mutations":  len(disagree_plddt),
        })
    if not per_protein_results:
        print("  no proteins with both agree and disagree mutations — skipping level 2 plot")
        return
    df_per_protein = pd.DataFrame(per_protein_results)
    out_csv = os.path.join(args.output_dir, "per_protein_plddt_delta.csv")
    df_per_protein.to_csv(out_csv, index=False)
    print(f"saved: {out_csv}")
    n_negative = (df_per_protein["delta_median_plddt"] < 0).sum()
    n_positive = (df_per_protein["delta_median_plddt"] > 0).sum()
    print(f"  proteins where disagree plddt < agree plddt: {n_negative}")
    print(f"  proteins where disagree plddt > agree plddt: {n_positive}")
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.histplot(df_per_protein["delta_median_plddt"], kde=True, ax=ax,
                 color=color_dist, edgecolor="white", bins=20, alpha=0.8)
    ax.axvline(0, color="#555555", linestyle="--", linewidth=1.5, label="no difference (0)")
    ax.axvline(df_per_protein["delta_median_plddt"].median(), color=color_disagree,
               linestyle="--", linewidth=1.5,
               label=f"median: {df_per_protein['delta_median_plddt'].median():.1f}")
    ax.set_xlabel("median pLDDT (disagree) − median pLDDT (agree)")
    ax.set_ylabel("number of proteins")
    ax.set_title(f"per-protein pLDDT difference: disagree vs agree mutations\n"
                 f"(n={len(df_per_protein)} proteins — "
                 f"{n_negative} negative, {n_positive} positive)")
    ax.legend(frameon=False)
    plt.tight_layout()
    out = os.path.join(args.output_dir, "per_protein_plddt_delta.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"saved: {out}")

    print("\ndone.")


if __name__ == "__main__":
    main()
