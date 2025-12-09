## Find overlap in ThermoMPNN training data and current MAVISp index.csv to select benchmarking proteins

# As a preliminary step, the datasets from ThermoMPNN were cleaned and Uniprot ACs added, as detailed in
# the directory cleanup_trainingdata. The resulting combined csv from the mega and fireprot datasets is
# here compared to the latest MAVISp index.csv, to find overlaps vs. not overlapping proteins, 
# which can be used for benchmarking.

# To run: 
./compare_ThermoMPNN_trainingdata.py

# Check if overlap.csv and benchmarking_dataset.csv were created
