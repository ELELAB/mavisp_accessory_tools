# FoldX Version Comparison Tool

## What it does

Compares protein stability predictions between FoldX5 and FoldX5.1.

For each protein, it calculates:
- Pearson correlation between the two versions
- Classification accuracy (stabilizing/neutral/destabilizing)
- Generates scatter plots and confusion matrices

## How to run

python foldx_version_comparison.py \
  --version1-dir /data/user/shared_projects/mavisp_ensemble_sim_length/foldx5.1_evaluation/foldx5_csv_repo/28112025_foldx5_candidates \
  --version2-dir /data/user/shared_projects/mavisp_ensemble_sim_length/foldx5.1_evaluation/foldx5.1_csv_repo/1712025_folx5.1_candidates \
  --output-dir ~/foldx_comparison_results

## Run tests

python test_foldx_comparison.py

## View results

cd ~/foldx_comparison_results

# rmsd_matrix.py

## What it does

Computes all-atom and Cα RMSD between the initial PDB structures used in the
FoldX5 and FoldX5.1 datasets. Structures are matched by UniProt ID and protein
name. RMSD is computed only on overlapping residues when the two versions cover
different regions of the same protein. Alignment is performed using the QCP
algorithm via MDAnalysis.

## Requirements

Python 3.10
MDAnalysis >= 2.7.0
numpy, pandas, matplotlib, seaborn

## How to run

python rmsd_matrix.py \
    -f /path/to/foldx5_initial_structures \
    -i /path/to/data_collection_foldx5.1 \
    -o ./rmsd_results

# Output

rmsd_results.csv                table with all-atom and Cα RMSD per structure pair
rmsd_distributions.png          histogram of RMSD values with mean and median lines
rmsd_scatter_allatom_vs_ca.png  scatter of all-atom vs Cα RMSD per structure
rmsd_barplot_sorted.png         bar plot of all structures sorted by Cα RMSD
rmsd_violin.png                 violin plot of RMSD distribution across all protein




