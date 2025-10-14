# MAVISp dot plot
*Last updated*: 22/04/2024
## Description

The script `piechart.py` takes in input the MAVISp table, the internal mavisp dictionary for clinvar records and returns a piechart plot of variants available the aggregated csv file, using the internal classification of MAVISp .

## Requirements

- python >= 3.7
- matplotlib
- pandas
- seaborn
- argparse
- argcomplete

## Inputs 

- MAVISp csv

## Usage

`python piechart.py [-h] -i INPUT -o OUTPUT NAME -d INTERNAL_MAVISP_DICT`

- `-i` : input MAVISp csv
- `-d` : input internal mavisp dictionary 
- `-o` : plot(s) output filename (.pdf/.png)

N.B
The internal mavisp dictionary can be found here in the official repository of Github of MAVISp : 
https://github.com/ELELAB/MAVISp/blob/main/mavisp/data/clinvar_interpretation_internal_dictionary.txt
