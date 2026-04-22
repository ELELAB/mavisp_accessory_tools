#We need to add a flag for isoform support
#the -s for isoform number
source /usr/local/envs/PDBminer/bin/activate
PDBminer -g SUMO3 -u P55854 -s 2 -n 1 -f csv
