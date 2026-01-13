## Cleanup of ThermoMPNN training data for selection of a benchmarking protein set

This folder contains scripts to clean and preprocess the FireProt and Mega datasets that were used for training **ThermoMPNN**.

The scripts isolate **PDB IDs** where available and fetch the corresponding **UniProt accession numbers (ACs)** from the RCSB PDB.  If no PDB ID is found, or if no UniProt AC is associated with the PDB entry, **BLAST** is used to identify matching human proteins with:

- ≥ 30 % sequence identity  
- ≥ 70 % sequence coverage  

If multiple entries in the query datasets correspond to different variants of the same protein, **only the first entry is taken into account**.

After processing, the cleaned datasets are written out as **separate CSV files**, as well as **combined into a single CSV file**.

**Note:**  These scripts do **not** need to be rerun every time a benchmarking dataset for ThermoMPNN is selected. They are provided for **documentation and reproducibility** primarily.

### To run cleanup for fireprot dataset:
./cleanup_fireprot.py

### To run cleanup for mega dataset:
./cleanup_mega_step1.py
./cleanup_mega_step2.py

### Optional arguments for cleanup_mega_step2.py
 `--identity` BLAST identity threshold in percent (default: 30%)
 `--coverage` BLAST alignment coverage in percent (default: 70%)


### To write csvs and combine the datasets:
./combine.py

### Sources
The training datasets for ThermoMPNN were sourced from https://github.com/Kuhlman-Lab/ThermoMPNN/tree/main/data_all/training.
For information on the training of ThermoMPNN, see:
Dieckhaus, M. Brocidiacono, N.Z. Randolph, & B. Kuhlman, Transfer learning to leverage larger datasets for improved prediction of protein stability changes, Proc. Natl. Acad. Sci. U.S.A. 121 (6) e2314853121, https://doi.org/10.1073/pnas.2314853121 (2024).

