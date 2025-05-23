Structure of folders to run the script should always be:
directory_name/Gene_name/AF_XX-YY/simulation_time/energies.csv

Requirements:
python/3.10

Run script:
bash run.sh input_directory simulation_time output_directory


Run script with python:
python ensemble_confusion_matrix.py input_directory simulation_time output_directory


Example:

bash run.sh runs 1000 output


or 



python ensemble_confusion_matrix.py runs 1000 output


The order of the arguments should always be like mentioned above.

NB: simulation_time argument is the gold standard simulation time with respect to which the output will be generated. 

