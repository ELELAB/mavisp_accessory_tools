#!/bin/bash

# Activate the environment
source /data/user/karok/coconat/coconat_env/bin/activate

# Requirements 
# Copy protein fasta file 

#Ab initio coiled-coil prediction 
python run_coconat_abinitio_docker.py --fasta_file $1 --output_file $2.tsv --plm_dir /data/user/EiriniGian/devel/coconat/env/coconat/coconat-plms

# Load python 
module load python
# Summarise results
python coco_sum.py -i $2.tsv -o $2_summary.tsv

