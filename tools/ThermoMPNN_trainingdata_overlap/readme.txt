## Find overlap in ThermoMPNN training data and current MAVISp index.csv to select benchmarking proteins

# As a preliminary step, the datasets from ThermoMPNN were cleaned and Uniprot ACs added, as detailed in
# the directory cleanup_trainingdata. The resulting combined csv from the mega and fireprot datasets is
# here compared to the latest MAVISp index.csv, to find overlaps vs. not overlapping proteins, 
# which can be used for benchmarking.

# To run: 
./compare_ThermoMPNN_trainingdata.py <index_csv_path> <training_csv_path>
# example: 
./compare_ThermoMPNN_trainingdata.py /data/raw_data/computational_data/mavisp_database_saturation/09122025_ALL/simple_mode/index.csv combined_mega_fireprot.csv 


# Check if overlap.csv and benchmarking_dataset.csv were created
