# in order to run AlphaFolder-Multimer you will need a FASTA file containing
# all the sequences of your proteins that will take part in the complex.
# this file is usually called input.fasta

# once this is ready, edit and read carefully the comment in the afmulti.sh
# script and edit what you need. Then run:

tsp -N 1 bash afmulti.sh

# one the run is finished run the bash script pdockq2.sh.  The script runs pDockQ2_i.py and 
# calculates for each model the pDockQ scores, in the output file pDockQ_scores.txt, and other metrics 
# (ptm, iptm, pLDDT, PAE, ranking confidence), in the output file detailed_scores.csv .
# N.B. run the bash script with tsp if your complex includes more than two chains

tsp -N 1 bash pdockq2.sh 

# after the run is finished you can run the clean.sh to extract PAE scores, compress files and remove pkl files

bash clean.sh
