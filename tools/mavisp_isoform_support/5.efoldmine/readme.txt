#Requirement:
# - fasta file of protein sequence


#Activate environment:
. /usr/local/envs/efoldmine/bin/activate


#Run efoldmine:
b2bTools -i protein.fasta -t protein.tabular -o protein.json --efoldmine


