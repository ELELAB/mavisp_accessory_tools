#have a fasta file with your full protein sequence 
#download using the isoform extension after the uniprot ac (eg. P55854-2)
wget https://www.uniprot.org/uniprot/$uniprot.fasta

#run with tsp

tsp -N 4 bash do.sh file.fasta
