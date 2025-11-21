module load python
pdb=$1
pdb_element $pdb.pdb > tmp.pdb
pdb_delelem -H tmp.pdb > noH.pdb
rm tmp.pdb
rm $pdb.pdb

