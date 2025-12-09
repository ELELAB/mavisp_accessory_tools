## Cleanup of ThermoMPNN training data for selection of benchmarking protein set

# These scripts clean the datasets from fireprot and mega databases that were used 
# for the training of ThermoMPNN, they isolate the PDB IDs if they are given, and fetch
# the according Uniprot AC from RCSB PDB. If no PDB ID is found or nor Uniprot AC is connected,
# BLAST is run to identify all entries with over 30 % identity with the query sequence. 
# If in the query datasets multiple entries correspond to different variants of a protein, 
# only the first entry will be taken into account. After this step, the two datasets are combined.

# These scripts DO NOT have to be rerun every time to select a benchmarking dataset for ThermoMPNN,
# but are for documentation.  

# To run cleanup for fireprot dataset:
./cleanup_fireprot.py

# To run cleanup for mega dataset:
./cleanup_mega_step1.py
./cleanup_mega_step2.py

# To combine both:
./combine.py
