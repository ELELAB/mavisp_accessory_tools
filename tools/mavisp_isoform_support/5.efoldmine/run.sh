. /usr/local/envs/efoldmine/bin/activate
uniprot=$1
cp ../4.demask/$uniprot.fasta .
#Run efoldmine:
b2bTools -i $uniprot.fasta -t $uniprot.tabular -o $uniprot.json --efoldmine
