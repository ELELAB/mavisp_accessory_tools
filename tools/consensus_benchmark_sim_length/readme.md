# Performance and confusion matrix for mavisp stability classification on ensembles with (Rosetta/FoldX) or (RaSP/FoldX) consensus depending on different simulation times

mavisp_sim_lengths_comparison.py is a Python script designed to generate confusion matrix plots and calculate various performance metrics for the stability classification of MAVISP across different ensemble datasets. These datasets are derived from simulations conducted at various time points and are compared against a gold-standard simulation time.

The script handles consensus classifications (Rosetta/FoldX and/or RaSP/FoldX), selecting the appropriate input files for each combination. The gold standard for each gene is the dataset corresponding to the simulation time specified via command-line flags, while other simulation times are used for comparison. For each gene, the script generates a confusion matrix and calculates various performance metrics.

The output includes confusion matrix plots, a CSV file containing the performance metrics, and the input CSV files used to generate these results. These are saved in directories named according to the simulation time used for comparison. Additionally, the script concatenates the datasets for all genes to create overall confusion matrix plots and performance metrics, which are saved in appropriately named folders.

## Requirements

The script requires the following Python version:

Python 3.10
The script also requires the following Python libraries:

pandas
csv
os
numpy (imported as np)
seaborn (imported as sn)
matplotlib.pyplot (imported as plt)
argparse
re
glob
sklearn

## Description

mavisp_sim_lengths_comparison.py takes as input a simulation time, which specifies the dataset to use as the gold standard, and processes all CSV files located in the same folder as the script. The script performs several checks and data preparation steps to ensure proper input handling and analysis.

### filtering step

First, the script verifies that the simulation time is in the correct format and that it appears in the name of at least one CSV file. It also checks that the CSV files adhere to the expected format. The files are then sorted based on the consensus classifications they contain (e.g., Rosetta/FoldX and/or RaSP/FoldX). CSV files that lack both consensus columns, or that belong to a gene set where the gold-standard dataset (corresponding to the specified simulation time) is missing both consensus columns, are discarded. Additionally, if any gene datasets are missing their gold standard, the entire set of datasets for that gene is removed.

After filtering, the script organizes the selected datasets into appropriate data structures for analysis. Two levels of comparison are performed:

Gene-specific comparison: Each gene's dataset is compared individually.
Comprehensive comparison: If more than one gene is analyzed, a combined dataset is created by concatenating the CSV files that passed the filtering step.
Datasets are parsed to remove any rows that contain NaN, Uncertain, or Stability classifications in the consensus columns.

### Confusion matrix and matrices 

For each dataset and consensus classification, the script generates a confusion matrix with the following categories:

True Positive (TP): The count of destabilizing effects in the gold-standard dataset classified as Destabilizing in datasets from other simulation times according to MAVISP consensus.
True Negative (TN): The count of neutral effects in the gold-standard dataset classified as neutral in datasets from other simulation times.
False Positive (FP): The count of neutral effects in the gold-standard dataset classified as destabilizing in datasets from other simulation times.
False Negative (FN): The count of destabilizing effects in the gold-standard dataset classified as neutral in datasets from other simulation times.

Along with the confusion matrix, the script calculates the follwoign performance metrics:

Sensitivity: TP / (TP + FN)
Specificity: TN / (TN + FP)
Accuracy: (TP + TN) / (TP + TN + FP + FN)
Precision: TP / (TP + FP)
F1 Score: 2 * (Precision * Sensitivity) / (Precision + Sensitivity)
These metrics provide insight into how well the other simulation times align with the gold standard.

Upon completion, the script generates the following output files:

Confusion matrix plots in both PDF and PNG formats.
A CSV file containing performance metrics.
The performance metrics CSV includes the following columns: sensitivity, specificity, accuracy, precision and F1 score.

An example of the CSV structure:

| sensitivity | specificity | accuracy | precision | F1 score | 
|:-----------:|:-----------:|:--------:|:---------:|:--------:|
|    1        |        1    |     1    |     1     |     1    |


## Input

The script processes CSV files located in the same folder as the script, expecting them to follow the naming format:
gene_name_simulation_time(ns)-ensemble_mode.csv

For example:
ACVR1B_1000ns-ensemble-mode.csv

### Required Flag
- **-t, --sim_time**
This flag specifies the simulation time (in nanoseconds) for the dataset to be used as the gold standard for generating the confusion matrix and performance metrics. The format should be {simtime}ns (e.g., 1000ns).


## Output

The script generates a folder for each consensus found in the csvs files. In each of these folders a directory indicating the type of comparison (gene-specific or chomprensive of all the datsets) is created. The output files are collected in appropiate subfolders indicating the simulation time of the dataset used as comparison with the golden standard.

Here is an example of the output structure, with each consensus (Rosetta/FoldX and RaSP/FoldX) containing datasets for three genes, each with simulation times of 100ns, 500ns, and 1000ns:
```
├── RaSP-FoldX_consensus
│   ├── gene_specific_comparison_1000ns
│   │   ├── ACVR1B
│   │   │   ├── 100ns
│   │   │   │   ├── ACVR1B_1000ns_VS_100ns.png
│   │   │   │   ├── ACVR1B_1000ns_VS_100ns_dataset.csv
│   │   │   │   └── performance_ACVR1B_1000ns_VS_100ns.csv
│   │   │   └── 500ns
│   │   │       ├── ACVR1B_1000ns_VS_500ns.png
│   │   │       ├── ACVR1B_1000ns_VS_500ns_dataset.csv
│   │   │       └── performance_ACVR1B_1000ns_VS_500ns.csv
│   │   ├── APOC1
│   │   │   ├── 100ns
│   │   │   │   ├── APOC1_1000ns_VS_100ns.png
│   │   │   │   ├── APOC1_1000ns_VS_100ns_dataset.csv
│   │   │   │   └── performance_APOC1_1000ns_VS_100ns.csv
│   │   │   └── 500ns
│   │   │       ├── APOC1_1000ns_VS_500ns.png
│   │   │       ├── APOC1_1000ns_VS_500ns_dataset.csv
│   │   │       └── performance_APOC1_1000ns_VS_500ns.csv
│   │   ├── B2M
│   │   │   ├── 100ns
│   │   │   │   ├── B2M_1000ns_VS_100ns.png
│   │   │   │   ├── B2M_1000ns_VS_100ns_dataset.csv
│   │   │   │   └── performance_B2M_1000ns_VS_100ns.csv
│   │   │   └── 500ns
│   │   │       ├── B2M_1000ns_VS_500ns.png
│   │   │       ├── B2M_1000ns_VS_500ns_dataset.csv
│   │   │       └── performance_B2M_1000ns_VS_500ns.csv
│   └── overall_comparison_1000ns
│       ├── 100ns
│       │   ├── 1000ns_VS_100ns_dataset.csv
│       │   ├── performance_1000ns_VS_100ns.csv
│       │   ├── whole datasets_1000ns_VS_100ns.png
│       │   └── whole_datasets_1000ns_VS_100ns.png
│       └── 500ns
│           ├── 1000ns_VS_500ns_dataset.csv
│           ├── performance_1000ns_VS_500ns.csv
│           ├── whole datasets_1000ns_VS_500ns.png
│           └── whole_datasets_1000ns_VS_500ns.png
├── Rosetta-FoldX_consensus
│   ├── gene_specific_comparison_1000ns
│   │   ├── ATG12
│   │   │   ├── 100ns
│   │   │   │   ├── ATG12_1000ns_VS_100ns.png
│   │   │   │   ├── ATG12_1000ns_VS_100ns_dataset.csv
│   │   │   │   └── performance_ATG12_1000ns_VS_100ns.csv
│   │   │   └── 500ns
│   │   │       ├── ATG12_1000ns_VS_500ns.png
│   │   │       ├── ATG12_1000ns_VS_500ns_dataset.csv
│   │   │       └── performance_ATG12_1000ns_VS_500ns.csv
│   │   ├── MAP1LC3B
│   │   │   ├── 100ns
│   │   │   │   ├── MAP1LC3B_1000ns_VS_100ns.png
│   │   │   │   ├── MAP1LC3B_1000ns_VS_100ns_dataset.csv
│   │   │   │   └── performance_MAP1LC3B_1000ns_VS_100ns.csv
│   │   │   └── 500ns
│   │   │       ├── MAP1LC3B_1000ns_VS_500ns.png
│   │   │       ├── MAP1LC3B_1000ns_VS_500ns_dataset.csv
│   │   │       └── performance_MAP1LC3B_1000ns_VS_500ns.csv
│   │   ├── MEFV
│   └── overall_comparison_1000ns
│       ├── 100ns
│       │   ├── 1000ns_VS_100ns_dataset.csv
│       │   ├── performance_1000ns_VS_100ns.csv
│       │   ├── whole datasets_1000ns_VS_100ns.png
│       │   └── whole_datasets_1000ns_VS_100ns.png
│       └── 500ns
│           ├── 1000ns_VS_500ns_dataset.csv
│           ├── performance_1000ns_VS_500ns.csv
│           ├── whole datasets_1000ns_VS_500ns.png
│           └── whole_datasets_1000ns_VS_500ns.png

```

## Usage 
```
module load python/3.10/modulefile
mavisp_sim_lenghts_comparison.py -t 1000ns -i csv/
```

N.B for an example run, try the run.sh script in the emample folder:
```
bash run.sh
```
(Change the simulation time in case it's needed)
