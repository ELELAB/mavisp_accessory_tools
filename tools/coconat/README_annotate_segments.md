# CoCoNat Coiled-Coil Prediction & Annotation Pipeline

It provides scripts to predict coiled-coil regions in protein sequences and annotate the output with segment positions, registers, and oligomeric states.
---

## Setup

1. Ensure you have Python 3.8+ and required dependencies (e.g., pandas) installed. 
2. Make sure you have the CoCoNat PLM models downloaded and available at `/path/to/coconat-plms` or your preferred directory. 
3. Place your protein FASTA files in a convenient folder, e.g., `example/` or your working directory.

---

## Running the pipeline manually (protein by protein)

### Step 1: Run CoCoNat ab initio prediction

python run_coconat_abinitio_docker.py --fasta_file=example/YOUR_PROTEIN.fasta \
--output_file=results/YOUR_PROTEIN.tsv --plm_dir=/path/to/coconat-plms

# Replace YOUR_PROTEIN with the protein ID of interest.
# This will generate residue-level coiled-coil predictions in a TSV file.


### Step 2: Annotate coiled-coil segments

python summarise_coconat_segments.py --input_file=results/YOUR_PROTEIN.tsv --output_file=results_annotated/YOUR_PROTEIN_annotated.tsv

# This script processes the raw CoCoNat output, detecting coiled-coil regions and annotating each segment with start/end positions, register sequence, and predicted oligomeric state.
