module load python/3.10/modulefile
uniprot=$1
./string2pdb $uniprot -a
