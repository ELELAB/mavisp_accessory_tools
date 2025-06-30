# CoCoNat Coiled-Coil Prediction & Annotation Pipeline

This tool provides scripts to predict coiled-coil regions in protein sequences and annotate the prediction output with summarized segment information. 
It processes raw CoCoNat TSV files to detect coiled-coil segments, summarizing each segment with residue ranges, predicted oligomeric states, and associated probabilities.

---

## Requirements

1. Ensure you have Python 3.8+ and required dependencies (e.g., pandas) installed. 
2. Make sure you have the CoCoNat PLM models downloaded and available at `/path/to/coconat-plms` or your preferred directory. 
3. Place your protein FASTA files in a convenient folder, e.g., `example/` or your working directory.

---
## Usage

Below is an example of running the pipeline step-by-step.

## Running the pipeline manually (protein by protein)

### Step 1: Run CoCoNat ab initio prediction

python run_coconat_abinitio_docker.py --fasta_file=example/YOUR_PROTEIN.fasta \
--output_file=results/YOUR_PROTEIN.tsv --plm_dir=/path/to/coconat-plms

# Replace `/path/to/coconat-plms` with the directory where you have the CoCoNat PLM models downloaded

# This step generates residue-level coiled-coil predictions in a TSV file.

### Step 2: Annotate coiled-coil segments

python summarise_coconat_segments.py --input_file=results/YOUR_PROTEIN.tsv --output_file=results_annotated/YOUR_PROTEIN_annotated.tsv


# The annotation script:

- Identifies coiled-coil segments based on CC_CLASS predictions.

- Summarizes each segment with start-end residue ranges (ResRange).

- Includes predicted oligomeric state (OligoState) and associated probability (POligo).

- Outputs a TSV file with these summarized segment annotations.

## Annotated Output Columns

| Column     | Description                                                                                                              |
|------------|--------------------------------------------------------------------------------------------------------------------------|
| ID         | Protein accession, as reported in the input FASTA file                                                                   |
| ResRange   | Residue range of the coiled-coil segment (e.g., 10-45)                                                                   |
| OligoState | Predicted oligomeric state for the segment, consistent across all residues in the helix. Possible values: A = antiparallel dimer, P = parallel dimer, 3 = trimer, 4 = tetramer |
| POligo     | Probability of the predicted oligomeric state (same for all residues in the segment)                                       |


## Notes

The script validates the input file for expected columns and reports errors if the file is missing or malformed.

If no coiled-coil segments are found, a warning is printed and no output file is created.

The output is formatted as a tab-separated file for easy integration with other analysis pipelines.


