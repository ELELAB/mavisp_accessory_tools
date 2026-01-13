## Find overlap in ThermoMPNN training data and current MAVISp index.csv to select benchmarking proteins

This script identifies proteins that **overlap** between the ThermoMPNN training dataset and a **MAVISp `index.csv`**, as well as proteins that are **not present in the ThermoMPNN training data**. The resulting non-overlapping proteins can be used as a benchmarking set.

As a preliminary step, the ThermoMPNN training datasets were cleaned and annotated with UniProt accession numbers, as described in the `cleanup_trainingdata` directory. The combined FireProt and Mega dataset generated there is used as the training reference in this comparison.


### To run:
./compare_ThermoMPNN_trainingdata.py -i <index_csv_path> -t <training_csv_path>

### Example:
./compare_ThermoMPNN_trainingdata.py -i /data/raw_data/computational_data/mavisp_database_saturation/09122025_ALL/simple_mode/index.csv -t cleanup_trainingdata/combined_mega_fireprot.csv


Check if the output files overlap.csv and benchmarking_dataset.csv were created.
