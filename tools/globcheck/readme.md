# GlobCheck: MAVISp Combined Globularity Metrics
Calculate globularity metrics for a trimmed protein domain PDB using MDAnalysis and CATH AlphaFlow.

The purpose of this script is to determined whether a structure/domain is suitable to use for 
the long_range module using AlloSigMA2. 

If yes, the structure will have 'Pass' in the 'Satisfied_thresholds' column of the output CSV. 

## Requirements

- python >= 3.7

## Usage 
usage: test.py [-h] -p PDB -r RANGE [-o OUTPUT]

Run globularity metrics using CATH AlphaFlow and MDAnalysis.

options:
  -h, --help            show this help message and exit
  -p PDB, --pdb PDB     Input trimmed PDB file
  -r RANGE, --range RANGE
                        Residue trimming range (e.g., 160-214)
  -o OUTPUT, --output OUTPUT
                        Output CSV file

## To run 
python globularity_metrics.py -p <input.pdb> -r <start-end> [-o <output.csv>]

## Output
CSV with asphericity, normalized radius of gyration, packing density, and pass/fail status.

## Pass criteria
Asphericity < 0.1
Normalized Rg < 0.356
Packing Density > 10.333
