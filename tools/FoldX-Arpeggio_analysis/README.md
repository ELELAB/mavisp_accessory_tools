# FoldX-Arpeggio analysis

## Introduction

The Snakefile implements a fully automated pipeline that takes as input the FoldX results organized into folders named by the wild-type (WT) residue chain and position (e.g., AA25, FA102, …).
For each residue, the workflow:

- extracts energetic contributions from FoldX, plotting them with bars plot for a better visualization

- Generates all mutations using the mutation list provided to FoldX

- Selects the corresponding FoldX models

- Runs Arpeggio on both WT and mutant structures

- Produces WT–mutant comparative interaction files highlighting interactions that are unique, lost, or altered.

- it generates “collapsed” summary files by aggregating contributions across all models/replicates for each mutation. 

## Affiliation

Cancer Structural Biology, Danish Cancer Society Research Center, 2100, Copenhagen, Denmark
Cancer Systems Biology, Department of Health and Technology, Section for Bioinformatics, 2800, Lyngby, Denmark

## Description

The Snakemake pipeline automates the analysis of mutation effects by combining energetic contributions from FoldX with geometrical distance analysis through Arpeggio. The Snakefile takes as input a configuration file specifying the FoldX output path, the mutation list, and the desired name for the output directories.
The workflow automatically parses the WT residue, position, and chain from the FoldX folder structure. It loads the list of mutant amino acid letters, generates all possible mutations in the format W25F, W25Y, etc., and sets up the output directories and parameters defined in the YAML configuration file.
The pipeline executes the script extract_foldx_energetic_contributions.py to parse the FoldX output, extracting for each mutation the energetic contributions of both the mutant and the WT across the five models generated during the FoldX run. It then computes, for each model, the difference between mutant and WT for each energetic term and aggregates these values into output file called delta.csv for inspection. Additionally, it identifies all energetic contributions where the mutant–WT subtraction results in a positive value (indicating a destabilizing effect) and produces a positives.csv file summarizing these contributions for all mutations.The pipeline also generates two additional files. The first called delta_mean.csv reports, for each energetic term, the difference between the mean value across the five mutant models and the mean value across the five WT models, together with the corresponding standard deviations. It then extracts those average values that remain positive (destabilizing after averaging across WT and mutant models) and outputs them together with the associated WT and mutant standard deviations in a file called delta_mean.csv. Finally, a file named statistics.csv summarizes information such as the number of residues involved in interfaces, clashes, and other relevant metrics.
Next, the pipeline processes the FoldX repair models for each mutant and the validated WT _repair.pdb structure, generating the necessary input files for running Arpeggio (see the mavisp-accessory_tools repository, https://github.com/ELELAB/mavisp_accessory_tools/). Arpeggio is executed on the WT model and on each of the five models for every mutant. After the analysis, the pipeline aggregates all contacts computed across the five models into a single file, distinguishing between contacts and clashes. It then extracts those interactions present in the mutants but absent in the WT, and those present in the WT but lost in the mutants, saving each category into separate files.

Finally, the pipeline collapses these files across models, retaining only those contacts and clashes consistently present in all five mutant models.

## Requirements

Required Software

- Snakemake ≥ 7
- Python ≥ 3.8, with:
- pandas
- pyyaml
- shutil
- subprocess
- re
- os

FoldX (already executed; pipeline starts from FoldX outputs)

Arpeggio installed in the active environment

## Input

The Snakefile takes as input a configuration file with the following structure:

**templates_path**: path to the folder containing the templates and scripts required by the pipeline

**foldx_mutations_input_directory**: path to the FoldX output directory, which must contain one subfolder per mutated position, including the WT and mutant FoldX models and the associated energetic contributions

**mutation_list**: list of mutations used for the FoldX run

**main_output_folder**: main directory where the pipeline will store outputs from both FoldX and Arpeggio analyses

**foldx_output_folder**: name of the subfolder to store the processed FoldX energetic readouts

**arpeggio_output_folder**: name of the subfolder collecting the Arpeggio contact analysis results

**scripts_parameters**: general parameters controlling Arpeggio execution and any additional scripts (for reference, see the ArpeggioInteraction pipeline)

### foldx_mutations_input_directory

The directory produced by FoldX and required by the pipeline must follow this structure:
(Example rooted at /data/raw_data/computational_data/mutatex_data/marinara/pole/complexes/DNA/9B8T_24-1198_retromut_pdbs/mutations/)

9B8T_Q07864_24-1198_filtered_model0_checked_Repair/
├── GA498
└── TA1090

Each subfolder must contain:

**pdb files**: PDB files for each mutation must be provided, including all five mutant structures generated by the five independent FoldX runs, as well as the corresponding WT structure (Repair.pdb).

```
    ├── 9B8T_Q07864_24-1198_filtered_model0_checked_Repair_1_0.pdb
    ├── 9B8T_Q07864_24-1198_filtered_model0_checked_Repair_1_1.pdb
    ├── 9B8T_Q07864_24-1198_filtered_model0_checked_Repair_1_2.pdb
    ├── 9B8T_Q07864_24-1198_filtered_model0_checked_Repair_1_3.pdb
    ├── 9B8T_Q07864_24-1198_filtered_model0_checked_Repair_1_4.pdb
    ├── 9B8T_Q07864_24-1198_filtered_model0_checked_Repair.pdb 
```

**Interfaces Residues files** Fxout-format files reporting the energetic contributions computed by FoldX for both the WT and the corresponding mutant.

```
    ├── Interface_Residues_9B8T_Q07864_24-1198_filtered_model0_checked_Repair_1_0_AC.fxout
    ├── Interface_Residues_9B8T_Q07864_24-1198_filtered_model0_checked_Repair_1_1_AC.fxout
    ├── Interface_Residues_9B8T_Q07864_24-1198_filtered_model0_checked_Repair_1_2_AC.fxout
    ├── Interface_Residues_9B8T_Q07864_24-1198_filtered_model0_checked_Repair_1_3_AC.fxout
    ├── Interface_Residues_9B8T_Q07864_24-1198_filtered_model0_checked_Repair_1_4_AC.fxout
    ├── Interface_Residues_WT_9B8T_Q07864_24-1198_filtered_model0_checked_Repair_1_0_AC.fxout
    ├── Interface_Residues_WT_9B8T_Q07864_24-1198_filtered_model0_checked_Repair_1_1_AC.fxout
    ├── Interface_Residues_WT_9B8T_Q07864_24-1198_filtered_model0_checked_Repair_1_2_AC.fxout
    ├── Interface_Residues_WT_9B8T_Q07864_24-1198_filtered_model0_checked_Repair_1_3_AC.fxout
    ├── Interface_Residues_WT_9B8T_Q07864_24-1198_filtered_model0_checked_Repair_1_4_AC.fxout
```

The mutation Interfaces Residues files and the WT ones are labeled with numbers: the first digit following the ‘Repair’ substring corresponds to the position of the mutation in the mutation list used by FoldX (i.e., ‘1’ refers to the first mutation in the list, ‘2’ to the second, and so on). The second digit identifies the index of the structural model (out of the five models generated by FoldX) corresponding to the five independent runs performed for that mutation.

YAML parameter file, example:

```
templates_path: templates
foldx_mutations_input_directory: "/data/raw_data/computational_data/mutatex_data/marinara/pole/complexes/DNA/9F6K_27-1175_pdbs/mutations/9F6K_Q07864_27-1175_model0_checked_Repair/"
mutation_list: "/data/raw_data/computational_data/mutatex_data/marinara/pole/complexes/DNA/9F6K_27-1175_pdbs/mutation_list.txt"
output:
   main_output_folder: "output"
   foldx_output_folder: "foldx_output"
   arpeggio_output_folder: "arpeggio_output"
scripts_parameters:
  ArpeggioInteractionPipeline:
    interface: False
    output_file: "arpeggio_output.csv"
    keep_proximal: False
```

### Template scripts

This folder contains the scripts used automatically by the pipeline and should not be modified.
For customization, refer to the config.yaml file provided to the Snakefile.
The ArpeggioInteractionPipeline.py script and its configuration file are already described in the corresponding repository ( https://github.com/ELELAB/mavisp_accessory_tools/).

#### extract_foldx_energetic_contributions.py

The script takes as input the foldx_mutations_input_directory specified in the configuration file.
It performs the following steps:

- Reads all FoldX .fxout files located inside the position-specific subfolders.
- Identifies and pairs WT and mutant structures based on their Repair_X_Y identifiers, using the mutation list as reference.
- Computes per-model energetic deltas (mut − wt) for all numeric FoldX columns.
- Removes structural count columns (e.g., interface residues, clashes) that are not relevant for energy analysis, saving them separately in the statistics.csv file.

The script generates the following output files:

- delta.csv,  per-model Δ values
- positives.csv, Δ values > 0 only (all others set to NA)
- delta_mean.csv, mean and standard deviation for WT and mutant groups, plus the mean delta
- statistics.csv, summary of structural statistics
- positive_mean.csv,  mean Δ > 0 (generated when --enable-positive-mean is enabled)


#### convert_dna_nomenclature.py

Used for Arpeggio analysis when the PDB file contains DNA nucleotides with a nomenclature not supported by Arpeggio.
Functionality:
- Reads a PDB file line by line.
- Processes only ATOM and HETATM records (or selected chains).
- Replaces residue names as follows:
   - A → DA
   - T → DT
   - C → DC
   - G → DG

Writes the modified PDB to a new file.

#### plot_foldX_deltas.py

This script generates per-mutation bar plots summarizing FoldX energy terms along with their propagated uncertainties.
For each mutation, the script takes as input 
a CSV file where each row corresponds to a mutation, and columns follow the FoldX naming scheme:

For each mutation, the script requires a CSV file where each row corresponds to a mutation, and columns follow the FoldX naming convention:

**Feature**, difference between Feature_mean_wt and Feature_mean_mut

**Feature_mean_wt**, mean value for the wild-type

**Feature_mean_mut** mean value for the mutant

**Feature_std**, propagated standard deviation from Feature_wt_std and Feature_mut_std

**Feature_wt_std**, standard deviation of the wild-type

**Feature_mut_std**, standard deviation of the mutant

**Position**, residue position

**Mutation**, mutation identifier

example of input file:

|Position|Mutation|Interaction_Energy|Interaction_Energy_wt|Interaction_Energy_wt|Interaction_Energy_wt_std|Interaction_Energy_mut_std|...|
|--------|--------|------------------|---------------------|---------------------|-------------------------|--------------------------|---|
|702	   |S       |	1.23	           |       -3.23         |   4.0               | 0.04	                   | 0.05                     |...|


The script extracts all relevant FoldX energy components.
It computes propagated errors from the WT and mutant standard deviations. Specifically, the standard deviation associated with the difference between WT_mean and Mutant_mean is calculated by propagating the WT and mutant errors using the formula:

$\sigma_X = \sqrt{\sigma_{\text{WT}}^2 + \sigma_{\text{MUT}}^2}$ 

Each component with a non-zero value is visualized as a bar with an error bar. The Interaction_Energy bar, representing the sum of all energy contributions, is placed as the first bar in the plot.

The script exports one .png file per mutation.

## Output

The pipeline produces an output folder named according to the string specified in the configuration file.
Inside this folder, there are two subfolders, each named based on strings from the configuration file: one contains the FoldX energy contribution analysis, and the other contains the Arpeggio outputs.

### FoldX energy contirbution output
The folder contains five files, including:

#### delta.csv
This file collects all differences for each energetic contribution reported as columns, comparing each FoldX-generated mutant model to the corresponding WT model. It has the following structure:

**Position**, the residue position mutated in the protein

**Mutation**, the mutation under investigation

**Model**, index of the model under investigation

**Feature**, difference between the mutant and WT values for a specific energetic contribution in that model

Here an exmaple of delta.csv:

```
|Position|Mutation|Model|IntraclashesGroup1|
|--------|--------|-----|------------------|
|GA702|S|1_0|0.054400000000001114|
|GA702|S|1_1|0.02790000000000248|

```

N.B In the Model column, the index refers to the model number as the second digit. The first digit refers to the mutation’s position in the mutation list; in the example, it is the first mutation in the list.

#### delta_mean.csv
This file collects, for each energetic contribution, the mean value across the five models for both mutant and WT, as well as the difference between these two means. It also reports the standard deviation of the WT and mutant values (calculated across the five models). The file has the following structure:

**Position**, residue position mutated in the protein

**Mutation**, mutation under investigation

**Model**, index of the model under investigation

**Feature_wt_mean**, mean value of the energetic contribution for the WT across the five models

**Feature_mut_mean**, mean value of the energetic contribution for the mutant across the five models

**Feature_wt_std**, standard deviation of the WT value across the five models

**Feature_mut_std**, standard deviation of the mutant value across the five models

**Feature**, difference between the mutant and WT mean values for that energetic contribution


```
|Position|Mutation|Backbone_Hbond_wt_mean|Backbone_Hbond_wt_std|Backbone_Hbond_mut_mean|Backbone_Hbond_mut_std|Backbone_Hbond|...|
|--------|--------|----------------------|---------------------|-----------------------|----------------------|--------------|...|
|GA702|S|-11.8202|0.0|-11.87902|0.0014855975228842203|-0.05882000000000076|...|

```
#### positives.csv
This file is a filtered version of delta.csv, containing only the positive contributions for each line. Negative values or zeroes are replaced with NaN automatically.

#### positives_mean.csv
This file is a filtered version of delta_mean.csv, containing only positive contributions for each line. Negative values or zeroes are replaced with NaN automatically.

#### statistics.csv
This file contains fields not directly linked to energetic contributions. It includes:

**Number_of_Residues_wt_mean**, number of residues in the WT

**Number_of_Residues_mut_mean**, number of residues in the mutant

**Number_of_Residues**, difference between Number_of_Residues_mut_mean and Number_of_Residues_wt_mean

**Interface_Residues_wt_mean**, number of residues at interfaces in the WT model

**Interface_Residues_mut_mean**, number of residues at interfaces in the mutant model

**Interface_Residues**, difference between mutant and WT interface residues

**Interface_Residues_Clashing_wt_mean**, number of interface residues clashing in the WT

**Interface_Residues_Clashing_mut_mean**, number of interface residues clashing in the mutant

**Interface_Residues_Clashing**, difference between mutant and WT clashing residues

**Interface_Residues_VdW_Clashing_wt_mean**, number of interface residues with VdW clashes in the WT

**Interface_Residues_VdW_Clashing_mut_mean**, number of interface residues with VdW clashes in the mutant

**Interface_Residues_VdW_Clashing**, difference between mutant and WT VdW clashing residues

**Interface_Residues_BB_Clashing_wt_mean**, number of interface residues with backbone clashes in the WT

**Interface_Residues_BB_Clashing_mut_mean**, number of interface residues with backbone clashes in the mutant

An example of statistics.csv

|Number_of_Residues_wt_mean|Number_of_Residues_mut_mean|Number_of_Residues|Interface_Residues_wt_mean|Interface_Residues_mut_mean|Interface_Residues|Interface_Residues_Clashing_wt_mean|Interface_Residues_Clashing_mut_mean|Interface_Residues_Clashing|Interface_Residues_VdW_Clashing_wt_mean|Interface_Residues_VdW_Clashing_mut_mean|Interface_Residues_VdW_Clashing|Interface_Residues_BB_Clashing_wt_mean|Interface_Residues_BB_Clashing_mut_mean|Interface_Residues_BB_Clashing|
|--------------------------|---------------------------|------------------|--------------------------|---------------------------|------------------|-----------------------------------|------------------------------------|---------------------------|---------------------------------------|----------------------------------------|-------------------------------|--------------------------------------|---------------------------------------|------------------------------|
|1080.0|1080.0|0.0|71.0|72.0|1.0|0.0|2.0|2.0|0.0|2.0|2.0|0.0|0.0|0.0|


Finally, the folder contains a subfolder called plots, which includes one bar plot per mutation. Each plot shows the difference between the mean WT and mutant values for every energetic contribution that is non-zero.

### Arpeggio output

The Arpeggio output folder contains subfolders for the WT and each mutation.
Each mutation subfolder contains five subfolders named model_0 to model_4, which include the Arpeggio output files for each model.
Each mutation folder also contains a folder aggregated_results with two subfolders: clashes and contacts.
The WT folder contains Arpeggio results for a single representative structure only (e.g., 9F6K_Q07864_27-1175_model0_checked_Repair.pdb). The arpeggio_output.csv in the WT folder is used for comparison against mutants.

The contents of the clashes and contacts folders are structured the same, except:

  - clashes contains contacts that introduce clashes.
  - contacts contains contacts that do not introduce clashes.

Inside each folder, there are three key files:

#### aggregated_{contacts/clashes}.csv
Aggregates the results from model_0 to model_4 for a specific mutation. Essentially, it concatenates the five Arpeggio CSV outputs, adding a model column to indicate the model index.

#### unique_{contacts/clashes}.csv
Reports only contacts or clashes that are unique to either the WT or the mutant. This file retains the same column structure as the original Arpeggio output, with two additional columns:

**model**  model index
**origin** specifies whether the unique contact/clash is from the mutant or WT

#### collapsed_unique_contacts.csv
Collapses the unique_{contacts/clashes}.csv across models, summing contacts or clashes that are identical across multiple models. The model column lists all models sharing the same contact or clash, separated by underscores (e.g., model_0_model_1_model_2_model_3_model_4).
Here an example of output for a run with the pipeline:

output/
├── arpeggio_output
│   ├── G702S
│   │   ├── aggregated_results
│   │   │   ├── clashes
│   │   │   │   ├── aggregated_clashes.csv
│   │   │   │   ├── collapsed_unique_clashes.csv
│   │   │   │   └── unique_clashes.csv
│   │   │   └── contacts
│   │   │       ├── aggregated_contacts.csv
│   │   │       ├── collapsed_unique_contacts.csv
│   │   │       └── unique_contacts.csv
│   │   ├── model_0
│   │   │   ├── 6TNYapt_78-1107_dna_model0_checked_Repair_1_0.cif
│   │   │   ├── 6TNYapt_78-1107_dna_model0_checked_Repair_1_0_converted.cif
│   │   │   ├── 6TNYapt_78-1107_dna_model0_checked_Repair_1_0_converted.pdb
│   │   │   ├── 6TNYapt_78-1107_dna_model0_checked_Repair_1_0.pdb
│   │   │   ├── ArpeggioInteractionPipeline.py
│   │   │   ├── arpeggio_output
│   │   │   │   └── 6TNYapt_78-1107_dna_model0_checked_Repair_1_0_converted_SER_702.json
│   │   │   ├── arpeggio_output.csv
│   │   │   ├── clashes_contacts.csv
│   │   │   ├── config.yaml
│   │   │   ├── interacting_entities_list.txt
│   │   │   ├── interacting_entities_saturation_mutlist.txt
│   │   │   └── residues.csv
│   │   ├── model_1
│   │   │   ├── 6TNYapt_78-1107_dna_model0_checked_Repair_1_1.cif
│   │   │   ├── 6TNYapt_78-1107_dna_model0_checked_Repair_1_1_converted.cif
│   │   │   ├── 6TNYapt_78-1107_dna_model0_checked_Repair_1_1_converted.pdb
│   │   │   ├── 6TNYapt_78-1107_dna_model0_checked_Repair_1_1.pdb
│   │   │   ├── ArpeggioInteractionPipeline.py
│   │   │   ├── arpeggio_output
│   │   │   │   └── 6TNYapt_78-1107_dna_model0_checked_Repair_1_1_converted_SER_702.json
│   │   │   ├── arpeggio_output.csv
│   │   │   ├── clashes_contacts.csv
│   │   │   ├── config.yaml
│   │   │   ├── interacting_entities_list.txt
│   │   │   ├── interacting_entities_saturation_mutlist.txt
│   │   │   └── residues.csv
│   │   ├── model_2
│   │   │   ├── 6TNYapt_78-1107_dna_model0_checked_Repair_1_2.cif
│   │   │   ├── 6TNYapt_78-1107_dna_model0_checked_Repair_1_2_converted.cif
│   │   │   ├── 6TNYapt_78-1107_dna_model0_checked_Repair_1_2_converted.pdb
│   │   │   ├── 6TNYapt_78-1107_dna_model0_checked_Repair_1_2.pdb
│   │   │   ├── ArpeggioInteractionPipeline.py
│   │   │   ├── arpeggio_output
│   │   │   │   └── 6TNYapt_78-1107_dna_model0_checked_Repair_1_2_converted_SER_702.json
│   │   │   ├── arpeggio_output.csv
│   │   │   ├── clashes_contacts.csv
│   │   │   ├── config.yaml
│   │   │   ├── interacting_entities_list.txt
│   │   │   ├── interacting_entities_saturation_mutlist.txt
│   │   │   └── residues.csv
│   │   ├── ...
│   └── WT
│       ├── 6TNYapt_78-1107_dna_model0_checked_Repair.cif
│       ├── 6TNYapt_78-1107_dna_model0_checked_Repair_converted.cif
│       ├── 6TNYapt_78-1107_dna_model0_checked_Repair_converted.pdb
│       ├── 6TNYapt_78-1107_dna_model0_checked_Repair.pdb
│       ├── ArpeggioInteractionPipeline.py
│       ├── arpeggio_output
│       │   └── 6TNYapt_78-1107_dna_model0_checked_Repair_converted_GLY_702.json
│       ├── arpeggio_output.csv
│       ├── clashes_contacts.csv
│       ├── config.yaml
│       ├── interacting_entities_list.txt
│       ├── interacting_entities_saturation_mutlist.txt
│       └── residues.csv
└── foldx_output
    ├── delta.csv
    ├── delta_mean.csv
    ├── plots
    │   └── GA702_S.png
    ├── positive_mean.csv
    ├── positives.csv
    └── statistics.csv

## How to Run

```
module load python/3.10/modulefile
snakemake -s Snakefile -c 1
```

to run the example:

```
bash run.sh
```


