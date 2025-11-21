#requirements: 
# 1) directories cl1, cl2, cl3
# 2) The cluster pdb file in each.

module load python/3.10/modulefile  

#WHAT DOES IT DO:
# for each cluster - 
# 1) Add chain if missing
# 2) Remove solvent and hydrogens
# 3) Fix aminoacid specific names: HIE/HIS, HID/HIS, HIP/HIS, LYN/LYS, ASH/ASP, GLH/GLU, CYX/CYS).
# 4) Ensure that each line contains the correct number of characters. 
# 5) Remove END/TER/ENDMOL lines. 

#if the script doesn't work one way would be to prepare the pdb manually in the cl folder:
pdb=$1
pdb_element $pdb.pdb > tmp.pdb
pdb_delelem -H tmp.pdb > noH.pdb
rm tmp.pdb
rm $pdb.pdb

tsp -N 4 bash run_analysis.sh chain cores
#example
#tsp -N 4 bash run_analysis.sh  A 4
