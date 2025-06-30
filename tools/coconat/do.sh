#!/bin/bash

# Activate the environment
conda activate coconat

# Assign command-line arguments to descriptive variables
INPUT_FASTA=$1
OUTPUT=$2
ANNOTATED_OUTPUT=$3

# Run CoCoNat on example FASTA
python run_coconat_abinitio_docker.py \
  --fasta_file=example_annotated/${INPUT_FASTA}.fasta \
  --output_file=example_annotated/${OUTPUT}.tsv \
  --plm_dir=/data/user/EiriniGian/devel/coconat/env/coconat/coconat-plms

# Annotate CoCoNat output with segment info
python summarise_coconat_segments.py \
  --input_file=example_annotated/${OUTPUT}.tsv \
  --output_file=example_annotated/${ANNOTATED_OUTPUT}.tsv

