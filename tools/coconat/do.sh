#!/bin/bash

# Run CoCoNat on example FASTA
python run_coconat_abinitio_docker.py --fasta_file=example_annotated/Q6ZNE5.fasta \
--output_file=example_annotated/ATG14_Q6ZNE5.tsv --plm_dir=/data/user/EiriniGian/devel/coconat/env/coconat/coconat-plms

# Annotate CoCoNat output with segment info
python summarise_coconat_segments.py --input_file=example_annotated/ATG14_Q6ZNE5.tsv \
--output_file=example_annotated/ATG14_Q6ZNE5_annotated.tsv
