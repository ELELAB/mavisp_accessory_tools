# RMSD of initial structures: FoldX5 vs FoldX5.1
# For each protein/domain in both datasets, computes all-atom and CA RMSD
# Matching is done by UniProt ID + protein name
# RMSD is computed only on overlapping residues between the two structures
# Run: python rmsd_matrix.py -f /data/user/shared_projects/mavisp_ensemble_sim_length/foldx5.1_evaluation/foldx5_initial_structures -i /data/user/shared_projects/mavisp_ensemble_sim_length/foldx5.1_evaluation/data_collection_foldx5.1 -o ./rmsd_results
import os
import glob
import argparse
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
import MDAnalysis as mda
from MDAnalysis.analysis import align
warnings.filterwarnings("ignore")


COLOR_VIOLIN_ALL  = "#C4A882"   # pale brown
COLOR_VIOLIN_CA   = "#B8A0C8"   # pale purple
COLOR_BAR_ALL     = "#B22222"   # firebrick
COLOR_BAR_CA      = "#4A90A4"   # muted teal
COLOR_DIST        = "#6A9CC4"   # soft blue

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

# 1. PDB COLLECTION

def get_residue_range(pdb_path):
    residues = set()
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                resnum = int(line[22:26].strip())
                residues.add(resnum)
    if not residues:
        return None, None
    return min(residues), max(residues)


def collect_foldx5_pdbs(foldx5_dir):               # FoldX5 folder structure: foldx5_dir/PROTEINNAME/UNIPROTID_START-END.pdb
                                                   # If a protein has multiple regions, each is stored separately
    pdbs = {}
    pattern = os.path.join(foldx5_dir, "*", "*.pdb")
    for pdb_path in glob.glob(pattern):
        protein_name = Path(pdb_path).parent.name
        uniprot_id   = Path(pdb_path).stem.split("_")[0]
        key = f"{uniprot_id}_{protein_name}"
        res_start, res_end = get_residue_range(pdb_path)
        region_key = f"{key}_{res_start}-{res_end}"
        pdbs[region_key] = {"path": pdb_path, "res_start": res_start, "res_end": res_end, "uniprot": uniprot_id, "protein": protein_name}
    print(f"[FoldX5]   Found {len(pdbs)} PDB files in {foldx5_dir}")
    return pdbs

def collect_foldx51_pdbs(foldx51_dir):             # Path: foldx51_dir/USER/PROTEIN/free/stability/mutatex_input_preparation/AFDB_XXX-YYY/model_v6/saturation/UNIPROTID_trimmed.pdb
    pdbs = {}
    pattern = os.path.join(foldx51_dir, "**", "*_trimmed.pdb")
    for pdb_path in glob.glob(pattern, recursive=True):
        if "mutatex_input_preparation" not in pdb_path:
            continue
        parts = Path(pdb_path).parts
        try:
            foldx51_parts = Path(foldx51_dir).parts
            relative_parts = parts[len(foldx51_parts):]
            protein_name = relative_parts[1]
            uniprot_id   = Path(pdb_path).stem.replace("_trimmed", "")
            res_start, res_end = get_residue_range(pdb_path)
            region_key = f"{uniprot_id}_{protein_name}_{res_start}-{res_end}"
            if region_key not in pdbs:
                pdbs[region_key] = {"path": pdb_path, "res_start": res_start, "res_end": res_end, "uniprot": uniprot_id, "protein": protein_name}
        except Exception as e:
            print(f"  [WARNING] Could not parse path {pdb_path}: {e}")
            continue
    print(f"[FoldX5.1] Found {len(pdbs)} PDB files in {foldx51_dir}")
    return pdbs


def group_by_protein(pdbs):
    grouped = {}
    for region_key, info in pdbs.items():
        base_key = f"{info['uniprot']}_{info['protein']}"
        grouped.setdefault(base_key, []).append(info)
    return grouped


def find_common_structures(foldx5_pdbs, foldx51_pdbs):
    # For each match, compute the overlapping residue range -  if no overlap exists, skip and report
    pairs = []
    skipped_no_match = []
    skipped_no_overlap = []

    foldx5_grouped  = group_by_protein(foldx5_pdbs)
    foldx51_grouped = group_by_protein(foldx51_pdbs)

    all_proteins = set(foldx5_grouped.keys()) | set(foldx51_grouped.keys())
    for base_key in sorted(all_proteins):
        if base_key not in foldx5_grouped or base_key not in foldx51_grouped:
            skipped_no_match.append(base_key)
            continue

        for info5 in foldx5_grouped[base_key]:
            for info51 in foldx51_grouped[base_key]:

                overlap_start = max(info5["res_start"], info51["res_start"])   # compute overlap
                overlap_end   = min(info5["res_end"],   info51["res_end"])
                if overlap_start > overlap_end:
                    skipped_no_overlap.append(
                        f"{base_key} ({info5['res_start']}-{info5['res_end']} vs {info51['res_start']}-{info51['res_end']})"
                    )
                    continue
                pairs.append({
                    "label":         f"{base_key}_{overlap_start}-{overlap_end}",
                    "pdb_foldx5":    info5["path"],
                    "pdb_foldx51":   info51["path"],
                    "overlap_start": overlap_start,
                    "overlap_end":   overlap_end,
                    "range_foldx5":  f"{info5['res_start']}-{info5['res_end']}",
                    "range_foldx51": f"{info51['res_start']}-{info51['res_end']}",
                })
    print(f"\nPairs with overlapping residues:  {len(pairs)}")
    print(f"Skipped - no match across datasets: {len(skipped_no_match)}")
    print(f"Skipped - no residue overlap:       {len(skipped_no_overlap)}")
    if skipped_no_match:
        print(f"  [No match examples]: {skipped_no_match[:5]}")
    if skipped_no_overlap:
        print(f"  [No overlap examples]: {skipped_no_overlap[:5]}")
    return pairs


# 2. RMSD CALCULATION

def compute_rmsd_pair(pair):
    # Load both PDBs into MDAnalysis
    # align.alignto runs QCP algorithm and returns RMSD
    # FoldX5 = mobile, FoldX5.1 = reference
    label         = pair["label"]
    pdb_mobile    = pair["pdb_foldx5"]
    pdb_ref       = pair["pdb_foldx51"]
    overlap_start = pair["overlap_start"]
    overlap_end   = pair["overlap_end"]

    sel = f"protein and resid {overlap_start}:{overlap_end}"
    sel_ca = f"protein and name CA and resid {overlap_start}:{overlap_end}"
    try:
        mobile = mda.Universe(pdb_mobile)
        ref    = mda.Universe(pdb_ref)
        # All-atom RMSD on overlapping residues
        mobile_sel = mobile.select_atoms(sel)
        ref_sel    = ref.select_atoms(sel)
        if len(mobile_sel) != len(ref_sel):
            print(f"  [SKIP all-atom] {label}: atom count mismatch ({len(mobile_sel)} vs {len(ref_sel)})")
            rmsd_all = np.nan
        else:

            _, rmsd_all = align.alignto(mobile, ref, select=sel, match_atoms=False)
        # CA-only RMSD on overlapping residues
        mobile_ca = mobile.select_atoms(sel_ca)
        ref_ca    = ref.select_atoms(sel_ca)
        if len(mobile_ca) != len(ref_ca):
            print(f"  [SKIP CA] {label}: CA count mismatch ({len(mobile_ca)} vs {len(ref_ca)})")
            rmsd_ca = np.nan
        else:
            _, rmsd_ca = align.alignto(mobile, ref, select=sel_ca, match_atoms=False)
        return rmsd_all, rmsd_ca
    except Exception as e:
        print(f"  [ERROR] {label}: {e}")
        return np.nan, np.nan

def run_rmsd_analysis(pairs):
    # Loop over all matched pairs and compute RMSD for each
    results = []
    n = len(pairs)
    for i, pair in enumerate(pairs):
        print(f"  [{i+1}/{n}] {pair['label']}", end=" ... ")
        rmsd_all, rmsd_ca = compute_rmsd_pair(pair)
        print(f"all-atom: {rmsd_all:.3f} Å | CA: {rmsd_ca:.3f} Å" if not np.isnan(rmsd_all) else "FAILED")
        results.append({
            "label":          pair["label"],
            "rmsd_all_atoms": rmsd_all,
            "rmsd_ca":        rmsd_ca,
            "overlap_range":  f"{pair['overlap_start']}-{pair['overlap_end']}",
            "range_foldx5":   pair["range_foldx5"],
            "range_foldx51":  pair["range_foldx51"],
        })
    return pd.DataFrame(results)


# 3. VISUALIZATION

def plot_rmsd_distributions(df, output_dir):
    # Plot 1: histogram + KDE for all-atom and CA RMSD with mean and median lines
    # Plot 2: scatter of all-atom vs CA RMSD per structure
    # Plot 3: bar plot of all structures sorted by CA RMSD
    os.makedirs(output_dir, exist_ok=True)
    df_clean = df.dropna(subset=["rmsd_all_atoms", "rmsd_ca"])
    df_clean = df_clean.copy()
    df_clean["label_clean"] = df_clean["label"].apply(lambda x: "_".join(x.split("_")[1:]))

    # Plot 1: distributions
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("RMSD distributions: FoldX5 vs FoldX5.1 initial structures", fontsize=13, fontweight="bold", y=1.01)
    for ax, col, title in zip(axes, ["rmsd_all_atoms", "rmsd_ca"], ["All-atom RMSD (Å)", "Cα RMSD (Å)"]):
        bin_edges = np.arange(0, df_clean[col].max() + 1, 1)
        sns.histplot(df_clean[col], ax=ax, color=COLOR_DIST, edgecolor="white", bins=bin_edges, alpha=0.8)
        ax.axvline(df_clean[col].median(), color="#CC0000", linestyle="--", linewidth=1.5, label=f"Median: {df_clean[col].median():.2f} Å")
        ax.axvline(df_clean[col].mean(), color="#555555", linestyle=":", linewidth=1.5, label=f"Mean: {df_clean[col].mean():.2f} Å")
        ax.set_xlabel(title)
        ax.set_ylabel("Number of entries")
        ax.legend(frameon=False)
    plt.tight_layout()
    out = os.path.join(output_dir, "rmsd_distributions.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {out}")

    # Plot 2: all-atom vs CA scatter with diagonal reference line
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(df_clean["rmsd_all_atoms"], df_clean["rmsd_ca"], color=COLOR_DIST, alpha=0.8, edgecolors="white", s=80)
    lims = [0, max(df_clean[["rmsd_all_atoms", "rmsd_ca"]].max()) * 1.05]
    ax.plot(lims, lims, color="#999999", linestyle="--", linewidth=1, label="y = x")
    ax.set_xlabel("All-atom RMSD (Å)")
    ax.set_ylabel("Cα RMSD (Å)")
    ax.set_title("All-atom vs Cα RMSD per structure")
    ax.legend(frameon=False)
    out = os.path.join(output_dir, "rmsd_scatter_allatom_vs_ca.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")

    # Plot 3: bar plot sorted by CA RMSD
    df_sorted = df_clean.sort_values("rmsd_ca", ascending=False).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(max(10, len(df_sorted) * 0.5), 5))
    x = np.arange(len(df_sorted))
    ax.bar(x - 0.2, df_sorted["rmsd_all_atoms"], width=0.4, label="All-atom", color=COLOR_BAR_ALL, alpha=0.85)
    ax.bar(x + 0.2, df_sorted["rmsd_ca"], width=0.4, label="Cα", color=COLOR_BAR_CA, alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(df_sorted["label_clean"], rotation=45, ha="right", fontsize=11)
    ax.set_ylabel("RMSD (Å)")
    ax.set_title("RMSD per structure (sorted by Cα RMSD)")
    ax.legend(frameon=False)
    plt.tight_layout()
    out = os.path.join(output_dir, "rmsd_barplot_sorted.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")


def plot_rmsd_violin(df, output_dir):
    # Violin plot: distribution of all-atom and CA RMSD across all proteins
    df_clean = df.dropna(subset=["rmsd_all_atoms", "rmsd_ca"])
    df_clean = df_clean.copy()
    df_clean["label_clean"] = df_clean["label"].apply(lambda x: "_".join(x.split("_")[1:]))
    df_melt = df_clean.melt(
        id_vars="label_clean",
        value_vars=["rmsd_all_atoms", "rmsd_ca"],
        var_name="type",
        value_name="rmsd"
    )
    df_melt["type"] = df_melt["type"].replace({"rmsd_all_atoms": "All-atom", "rmsd_ca": "Cα"})
    fig, ax = plt.subplots(figsize=(6, 6))
    sns.violinplot(data=df_melt, x="type", y="rmsd", ax=ax,
                   palette={"All-atom": COLOR_VIOLIN_ALL, "Cα": COLOR_VIOLIN_CA},
                   inner="box", cut=0)
    ax.set_xlabel("")
    ax.set_ylabel("RMSD (Å)")
    ax.set_title("RMSD distribution across all proteins\nFoldX5 vs FoldX5.1")
    plt.tight_layout()
    out = os.path.join(output_dir, "rmsd_violin.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")


# 4. MAIN

def main():
    parser = argparse.ArgumentParser(description="Compute RMSD between FoldX5 and FoldX5.1 initial structures")
    parser.add_argument("-f", "--foldx5_dir", required=True, help="Path to FoldX5 folder (organized by protein name)")
    parser.add_argument("-i", "--foldx51_dir", required=True, help="Root of FoldX5.1 data collection folder")
    parser.add_argument("-o", "--output_dir", default="./rmsd_results", help="Where to save results (default: ./rmsd_results)")
    args = parser.parse_args()

    print("RMSD analysis: FoldX5 vs FoldX5.1 initial structures")
    print("\n[1] Collecting PDB files...")
    foldx5_pdbs  = collect_foldx5_pdbs(args.foldx5_dir)
    foldx51_pdbs = collect_foldx51_pdbs(args.foldx51_dir)
    print("\n[2] Matching structures and computing overlaps...")
    pairs = find_common_structures(foldx5_pdbs, foldx51_pdbs)
    if not pairs:
        print("\nERROR: No matching pairs found. Check folder structure and naming.")
        return

    # Compute RMSD for each pair
    print(f"\n[3] Computing RMSD for {len(pairs)} structure pairs...")
    df = run_rmsd_analysis(pairs)
    os.makedirs(args.output_dir, exist_ok=True)
    csv_out = os.path.join(args.output_dir, "rmsd_results.csv")
    df.to_csv(csv_out, index=False)
    print(f"\nSaved results table: {csv_out}")
    print("\n[4] Summary statistics:")
    print(df[["rmsd_all_atoms", "rmsd_ca"]].describe().round(3).to_string())
    print("\n[5] Generating plots...")
    plot_rmsd_distributions(df, args.output_dir)
    plot_rmsd_violin(df, args.output_dir)
    print("\nDone.")
if __name__ == "__main__":
    main()
