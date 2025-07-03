# CoCoSum: CoCoNat Coiled-Coil Prediction & Summarisation 
This tool uses CoCoNat (https://doi.org/10.1093/bioinformatics/btad495) to predict coiled-coil regions for a given protein. 

The per-residue TSV output of CoCoNat is then summarised by the coco_sum.py script, providing segment-wise information;
residue ranges, predicted oligomeric states, and associated probabilities.

The script validates the input file for expected columns and reports errors if the file is missing or malformed.

If no coiled-coil segments are found, a warning is printed and no output file is created.

The output is formatted as a tab-separated file for easy integration with other analysis pipelines.

## Output 
| Column     | Description                                                                                                              |
|------------|--------------------------------------------------------------------------------------------------------------------------|
| ID         | Protein accession, as reported in the input FASTA file                                                                   |
| ResRange   | Residue range of the coiled-coil segment (e.g., 10-45)                                                                   |
| OligoState | Predicted oligomeric state for the segment, consistent across all residues in the helix. Possible values: A = antiparallel dimer, P = parallel dimer, 3 = trimer, 4 = tetramer |
| POligo     | Probability of the predicted oligomeric state (same for all residues in the segment)                                       |

## Requirements
- python >= 3.7
- pandas
- CoCoNat (see https://github.com/BolognaBiocomp/coconat)

## Usage
usage: coco_sum.py [-h] -i INPUT_FILE -o OUTPUT_FILE

Summarise CoCoNat coiled-coil segments.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input_file INPUT_FILE
                        Path to CoCoNat output
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Path to save the summary file

## Example
cd example
bash do.sh Q6ZNE5.fasta ATG14_Q6ZNE5