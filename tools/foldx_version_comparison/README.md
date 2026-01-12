# FoldX Version Comparison Tool

# What it does

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



