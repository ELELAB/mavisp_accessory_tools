input_directory=$1
simulation_time=$2
output_directory=$3
module load python/3.10/modulefile
python ../ensemble_confusion_matrix.py $input_directory $simulation_time $output_directory 
