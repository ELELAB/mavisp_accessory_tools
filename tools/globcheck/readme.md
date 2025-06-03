# GlobCheck: MAVISp Combined Globularity Metrics
Calculate globularity metrics for a trimmed protein domain PDB using MDAnalysis and CATH AlphaFlow.

In particular, this script calculates ...
  - asphericity (https://pubs.acs.org/doi/10.1021/jp037128y)
  - normalised radius of gyration (10.1016/j.jmb.2024.168551)
  - packing density (10.1016/j.jmb.2024.168551)

The final purpose of this script is to determine whether a structure/domain is suitable to use for 
the long_range module using AlloSigMA2. In order to do that, the three aforamentioned measures
are calculated on the input structure, and if the structure is deemed satisfactory considering all
three of them, it is deemed to be globular enough to be used in the long_range workflow.
If this is the case, the structure will have 'True' in the 'satisfies_thresholds' column of the output CSV.

The thresholds we have selected are:

Asphericity < 0.1
Normalized Rg < 0.356
Packing Density > 10.333

The threshold for asphericity was determined by observing the distribution of asphericity values for single-chain globular proteins in (https://doi.org/10.48550/arXiv.q-bio/0310023).
Whereas, the thresholds for normalised radius of gyration and packing density were determined in (DOI: 10.1126/science.adq4946), which established these thresholds based on an analysis of over 193 million TED-100 domains deemed to be of sufficient structural quality.

The script considers the whole protein structures by default, but a custom trimming range can
be specified instead with option -r (start-end). In this case, the structure will be processed
by keeping only the residues in the trimming range before calculating th emeasures.

Custom thresholds for our measures can also be specified on the command line.

## Requirements

- python >= 3.7
- pandas
- MDAnalysis
- BioPython
- cath_alphaflow (see https://github.com/UCLOrengoGroup/cath-alphaflow)

## Usage 
usage: globcheck.py [-h] -p PDB [-r RANGE] [-o OUTPUT] [-a ASPHERICITY] [-pd PACKING] [-rg RADIUS]

Run globularity metrics using CATH AlphaFlow and MDAnalysis.

options:
  -h, --help            show this help message and exit
  -p PDB, --pdb PDB     Input trimmed PDB file (default: None)
  -r RANGE, --range RANGE
                        Residue trimming range (e.g., 160-214) (default: None)
  -o OUTPUT, --output OUTPUT
                        Output CSV file (default: None)
  -a ASPHERICITY, --asphericity ASPHERICITY
                        Threshold for asphericity (default: 0.1)
  -pd PACKING, --packing PACKING
                        Threshold for packing density (default: 10.333)
  -rg RADIUS, --radius RADIUS
                        Threshold for normalised radius of gyration (default: 0.356)

## To run 
python globularity_metrics.py -p <input.pdb> -r <start-end> [-o <output.csv>]

## Output
CSV with asphericity, normalized radius of gyration, packing density, and pass/fail status.


