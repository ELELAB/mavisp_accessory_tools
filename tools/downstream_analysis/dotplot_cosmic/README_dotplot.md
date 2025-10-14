# MAVISp dot plot
*Last updated*: 26/09/24
## Description

The script `dot_plot.py` takes in input the MAVISp table and returns a dotplot and a summary of the effect of the mutations on the single modules of the MAVISp framework, dividing them in structural/functional effects.

Additionally, a csv file containing mutations with an identified module effect and a 'pathogenic' AlphaMissesnse classification.
The last column of the file summarises the MAVISp identified effects for the given mutation 'MAVISp Effects'.

## Requirements

- python >= 3.7
- matplotlib
- numpy
- pandas

## Inputs 

- MAVISp csv
- dictionary.csv (Internal dictionary for ClinVar Annotations)

## Usage

`python dot_plot.py [-h] -i INPUT [-o OUTPUT] [-m MUTATIONS [MUTATIONS ...]] [-r RESIDUES [RESIDUES ...]] [-R REVEL_THRESHOLD] [-D DEMASK_THRESHOLD] [-G GEMME_THRESHOLD] [-x X_LIM] [-f FIGSIZE FIGSIZE] [-pltR] [-pltD]
[-pltC {all,uncertain,benign,likely_benign,pathogenic,likely_pathogenic,conflicting} [{all,uncertain,benign,likely_benign,pathogenic,likely_pathogenic,conflicting} ...]][-colC]`

- `-i` : input MAVISp csv
- `-o` : plot(s) output filename (default: 'dot_plot_#.[png/pdf]') 
- `-m` : selected mutations to be plotted to be provided comma-separated (e.g., A4G,F55K) 
- `-r` : selected residues to be plotted to be provided comma-separated (e.g., 4,55) 
- `-R` : REVEL score threshold (default: 0.5)
- `-D` : Threshold to classify a mutation according to the DeMask score. (Default = 0.25)
- `-G` : Threshold to classify a mutation according to the GEMME score. (Default = 3)
- `-x` : x axis limit (i.e., number of mutations to be plotted on the x-axis) (default: 15); multiple plots will be output if needed
- `-f` : figure size (default: 14 4). It is suggested to use the default figsize with 40/50 mutations on the x-axis and 7/8 labels on the y-axis
- `-pltR` : include in the plot the  Revel scores. (Default = None)
- `-pltD` : include in the plot the DeMaSk predicted consequence (LOF/GOF) for mutations satisying the DeMaSk threshold. (Default = None)
- `-pltC` : select mutations to be plotted according to the ClinVar Interpretation (e.g pathogenic, uncertain); multiple ClinVar categories can be included. Required additional input: dictionary.csv (Default = None)
- `-colC` : Color xtick labels according to the ClinVar Interpretation of the given mutation. Required additional input: dictionary.csv (Default=None)

## Example

Please, see the `example` directory and the `do.sh` script within.

## Output

The script returns 3 types of output:
- dot plots (`dot_plot_#`) (the number depends on the number of mutations and the x-axis limit). Each dot is colored according to the type of effect as described in the legend. The plots are saved in `pdf` as default and `png` format if `-m`/`-r` flags are used.
- a summary file `log.txt` containing useful information to summarize the effects of each module (GitBook compatible)
- a csv file `alphamissense_out.csv` containing mutations with identified module effect and AlphaMissense 'pathogenic' classification

N.B. If `-m`/`-r` are used to delineate specific mutations/residues the plots will be output as `png` and `alphamissense_out.csv` will be filtered. 