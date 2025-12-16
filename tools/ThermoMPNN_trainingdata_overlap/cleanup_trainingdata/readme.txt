## Cleanup of ThermoMPNN training data for selection of benchmarking protein set

# These scripts clean the datasets from fireprot and mega databases that were used 
# for the training of ThermoMPNN, they isolate the PDB IDs if they are given, and fetch
# the according Uniprot AC from RCSB PDB. If no PDB ID is found or nor Uniprot AC is connected,
# BLAST is run to identify all entries with over 30 % identity and 70 % coverage with the query sequence. 
# If in the query datasets multiple entries correspond to different variants of a protein, 
# only the first entry will be taken into account. After this step, the two datasets are written out 
# as separate csvs, and also combined into one csv.

# These scripts DO NOT have to be rerun every time to select a benchmarking dataset for ThermoMPNN,
# but are for documentation.  

# To run cleanup for fireprot dataset:
./cleanup_fireprot.py

# To run cleanup for mega dataset:
./cleanup_mega_step1.py
./cleanup_mega_step2.py
# Optional arguments for cleanup_mega_step2.py
# `--identity` BLAST identity threshold in percent (default: 30 %)
# `--coverage` BLAST alignment coverage in percent (default: 70 %)


# To write csvs and combine the datasets:
./combine.py
