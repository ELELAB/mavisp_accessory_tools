# LolliPlot
## Description
The script `lolliplot.py` takes as input and plots the output csv from dotplot.py named `alphamissense_out.csv`, containing the mutations with an AlphaMissense 'pathogenic' classification and an identified MAVISp effect. 

The plot groups the identified effects in the following broad categories: stability, ptm, long range, local int, functional.

## Requirements
module load python/3.10/modulefile 

## Usage
usage: lolliplot.py [-h] -i INPUT_FILE [-x XLIM] [-s]

Options:
- i INPUT_FILE  AlphaMissense output csv of dotplot.py (**required**)
- x XLIM        Number of mutations to plot on the x axis. (default 15)
- s             Save individual plots as PNGs

## Output 
The script returns the following:
- `lolliplot.pdf` containing # number of individual plots, depending on the xlim used. Each lolli is colored according to the type of effect as described in the legend.
- if `-s` then individual lolliplots (`lolliplot_#`) are saved in `png` format, in addition to the PDF.

## Example
Example generated using dotplot.py for ARID3A-ensemble_mode.csv
