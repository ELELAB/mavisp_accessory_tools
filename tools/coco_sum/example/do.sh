#!/bin/bash

# Load python 
module load python

# Requirements 
# Copy protein fasta file 

#Ab initio coiled-coil prediction 
python ../run_coconat_abinitio_docker.py --fasta_file $1 --output_file $2.tsv --plm_dir /usr/local/coconat/
# Load python 
module load python
# Summarise results
python ../coco_sum.py -i $2.tsv -o $2_summary.tsv

