#activate the python environment (module load python)
module load python
uniprot=$1
#starting=$2
#ending=$3
ln -s ../pdbminer/results/$1/$1_all.csv
#to run for mavisp - replace uniprot ID with the one of your protein
python find_PDBminer_complexes.py -i $1_all.csv --binding_interface -d 10 -o $1_filtered.csv
#python find_PDBminer_complexes.py -i $1_all.csv --binding_interface -d 10 -s $2 -e $3 -o $1_filtered.csv
