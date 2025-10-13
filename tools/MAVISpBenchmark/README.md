# MAVISpBenchmark.py

MAVISpBenchmark.py provides a comprehensive framework for analyzing and comparing MAVISp classifications and variant effect predictors (currently GEMME and DeMaSK) against experimental data. The script can operate either by using MAVISp classification columns directly or by combining multiple columns to generate ensemble classifications. The combination logic is fully configurable through a YAML file, which supports three alternative approaches—voting, priority, and weighted schemes—allowing for a flexible and detailed classification system tailored to specific analytical needs.

As output, the pipeline produces informative plots and performance metrics for DeMaSK and GEMME individually, evaluates MAVISp’s classification performance, and enables direct comparison of DeMaSK and GEMME on MAVISp-filtered datasets. In addition, it analyzes mutations with divergent classifications, highlighting discrepancies between predictor outputs and experimental reference data.

## Description

The script accepts as input one or multiple mavisp.csv files, in either simple or ensemble mode, along with a configuration file (config.yaml) that specifies the analysis mode and required parameters. Optional flags can be used to:

include or exclude proline and glycine residues,

specify a colormap for histograms and confusion matrices,

exclude specific residues through a dedicated file when structural problems prevent reliable MAVISp calculations.

The script first verifies that all required inputs (CSV, configuration, and mode) are provided and validates the format of any optional exclusion file. The configuration file defines which columns should be analyzed and plotted. One column must contain the binary experimental classification (Damaging or Neutral), either derived directly from experimental scores or from a threshold-based classification. If continuous scores are provided, the corresponding threshold(s) must also be specified.

The configuration further defines normalization rules, mapping predictor-specific outputs into standardized categories (Damaging, Neutral, or Uncertain). It also specifies which MAVISp column(s) should be compared—either single columns or combined sets—along with the classification scheme to be applied. A validation step ensures consistency between the configuration and the input CSV file.

After error handling, the script merges and organizes CSV files according to the selected mode (simple or ensemble) and begins processing. GEMME and DeMaSK scores are normalized according to predefined thresholds:

**DeMaSK**

- Damaging: score < -0.25 or score > 0.25 (loss/gain of function)

- Neutral: -0.25 ≤ score ≤ 0.25

**GEMME**

- Damaging: score < -3 or score > 3

- Neutral: -3 ≤ score ≤ 3

Additional normalization is performed based on the rules defined in config.yaml. A new classification column is then generated, either from a single MAVISp column or from the combination of multiple columns according to the specified scheme.

Subsequently, the script compares GEMME, DeMaSK, and MAVISp classifications against experimental data, after removing classifications that are not labeled as Damaging or Neutral in both the experimental classification column and the predictor column being compared, generating several types of plots.

### Confusion Matrices

Three confusion matrices are generated:

- Comparing GEMME and DeMaSK classifications with experimental data (gold standard).

- Comparing MAVISp classification (single or combined column, as specified in the config) with experimental data.

- Comparing GEMME and DeMaSK classifications with experimental data on a subset of variants used in the MAVISp–experiment comparison (excluding variants labeled as “Uncertain” by MAVISp). This enables direct performance comparison across MAVISp, GEMME, and DeMaSK.

For each comparison, the script also saves the underlying data to confusion_matrix_data.csv.

Classification criteria:

**True Positive (TP)**: Damaging according to predictor/MAVISp and confirmed as Damaging experimentally.

**True Negative (TN)**: Neutral according to predictor/MAVISp and confirmed as Neutral experimentally.

**False Positive (FP)**: Damaging according to predictor/MAVISp but Neutral experimentally.

**False Negative (FN)**: Neutral according to predictor/MAVISp but Damaging experimentally.

From these values, the script computes standard performance metrics (saved in metrics.csv):

**Sensitivity**: TP / (TP + FN)

**Specificity**: TN / (TN + FP)

**Accuracy**: (TP + TN) / (TP + TN + FP + FN)

**Precision**: TP / (TP + FP)

**F1 Score**: 2 × (Precision × Sensitivity) / (Precision + Sensitivity)

### Scatter Plots

For each MAVISp column specified in the configuration, the script plots GEMME or DeMaSK scores (x-axis) against experimental scores (y-axis), with points colored by MAVISp classification. Spearman correlation coefficients are reported to quantify the relationship between predictor scores and experimental data. Plot data is saved in scatter_plot_data.csv. The plots are saved in the same location of the first the confusion matrix (see section above and output section below)

### Histograms of Performance

An Histograms are generated from the results of the second and third confusion matrices (see section above), allowing direct comparison of predictor performance across the same variant datasets.

### Distribution Plots and Hisotgrams plots of aminoacid composition

The script identifies mutations where MAVISp classifications diverge from experimental data, saving them in divergent_classifications.csv.
By leveraging this dataset, divergent mutations are visualized on distribution plots of the experimental data. The experimental threshold is shown as a vertical dashed line. Two additional CSV files are produced from divergent_classifications.csv:

discordant_classification_neutral.csv containing divergent Neutral mutations (experimentally Neutral, classified Damaging by MAVISp),

discordant_classification_damaging.csv divergent Damaging mutations (experimentally Damaging, classified Neutral by MAVISp).

For each subset, the script generates:

Histograms of amino acid residues: showing the percentage distribution of residues involved in the dataset of discordant mutations. Two perspectives are considered:
(i) wild-type residues, to identify which amino acids are most frequently affected by discordant classifications when mutated;
(ii) mutant residues, to highlight which amino acids, when introduced as mutations, are most often associated with discordant classifications.

Distribution plots of solvent accessibility values: assessing whether discordant classifications preferentially occur in mutations located on solvent-exposed or buried residues.

Histograms of amino acid class types: showing the percentage distribution of chemical classes in the dataset of discordant mutations. Again, both wild-type residues (to identify which classes are most frequently affected when mutated) and mutant residues (to determine which classes are most often associated with discordant classifications) are analyzed.

## Requirements

Python >= 3.8

Packages:

- pandas
- numpy
- matplotlib
- seaborn
- scipy
- scikit-learn
- pyyaml

## Input

The script accepts MAVISp CSV files in either simple or ensemble mode, along with a YAML configuration file. The configuration specifies the two columns to be plotted and includes a section defining normalization rules, which map predictor-specific outputs into standardized categories (Damaging, Neutral, or Uncertain). It also indicates which MAVISp column(s) should be compared—either individual columns or combined sets—together with the classification scheme to be applied. 

### Config file

The config file has the following structure:

```
modes:
  <mode_name>:               # e.g., "simple" or "ensemble
    experimental_column: <string>
    experimental_column_classification: <string>
    experimental_thresholds: <list of numbers> # not need in case of classification in experimental_column
    normalization:
      Damaging: <list of strings>   # all labels considered "Damaging"
      Neutral: <list of strings>    # all labels considered "Neutral"
      Uncertain: <list of strings>  # all labels considered "Uncertain"
    comparisons:
      - name: <comparison_name>
        columns: <list of strings>
        priority:
          class_order: ["Damaging", "Neutral", "Uncertain"]
          priority_for_classification:
            Damaging: <list of columns>   # columns used to classify Damaging
            Neutral: <list of columns>    # columns used to classify Neutral
            Uncertain: <list of columns>  # columns used to classify Uncertain
        weighted:                      # optional, if using weighted classification
          weights:
            <column_name>: <weight_value>
            ...
          threshold: <numeric_value>
        voting:                        # optional, if using voting logic
          target_class: <class_name>   # e.g., "Damaging"
          logic: <string>              # e.g., "majority", "at_least_one" or "all"
          fallback: <class_name>       # fallback class
          handle_uncertain: <string>   # e.g., "return_uncertain",  "as_fallback"
      - name: <comparison_name>
        ...
      # additional comparisons possible
```
Here an explanation of all the fields:

#### modes
This field can only take the values simple or ensemble, depending on the MAVISp CSV file of interest. If the analysis is to be performed on both simple and ensemble MAVISp CSV files, the configuration file must include two separate sections (one for simple and one for ensemble), each specifying the appropriate columns of interest for the corresponding file type.

#### experimental_column
The column name containing the raw scores from the experimental assay.

#### experimental_column_classification
The column name containing the categorical classification derived from the experimental assay.

#### experimental_thresholds
A list of thresholds used to map the values in the experimental_column into categories in the experimental_column_classification column. These thresholds are also displayed on the distribution plot to provide a clearer indication of where the Damaging and Neutral classifications are located.

#### normalization

This field defines the rules for mapping predictor-specific and experimental outputs into standardized categories (Damaging, Neutral, or Uncertain). This normalization ensures that all classifications provided in the MAVISp CSV file are harmonized before performing the analysis.

The field is structured as follows:

```
Damaging:
  - "Damaging"
  - "damaging"
  - "Destabilizing"
  - "destabilizing"
  - "pathogenic"
  - "Pathogenic"
Neutral:
  - "Neutral"
  - "neutral"
  - "Benign"
  - "benign"
Uncertain:
  - "Uncertain"
  - "uncertain"
```

The keys (Damaging, Neutral, Uncertain) represent the standardized categories that will replace the corresponding values found in the columns specified in the configuration file, as well as in the GEMME and DeMaSK classification columns, for both MAVISp module columns and the experimental classification column.
The values listed under each key must include all possible labels that should be converted into that standardized category for the column of interest.

#### comparisons

This section includes the column that needs to be comapred aginst the expeirmental data (a part form the GEMME and DeMaSK columns that will be peroformed automatically).
The comparison can include as much as comparisons desisred, so potentially it woule be possible ot perofrm several comaprison against experimental data in one run. In order to perofrm a comaprisn the follwoing information must be provided:

##### name 

The name of the output column that will be used to generate legend titles, axis labels, and output folder paths. It should be as informative as possible about the column’s content, while remaining concise to avoid excessively long labels in plots and folder names. The output CSV file generated from it will include an additional column named after this field, containing the classification used for the comparison.

##### columns
A list of one or more columns from the MAVISp CSV file to be used for the comparison.
It can be used in two ways:

 - **one column** Only the classification from that column will be used for the comparison. An additional column (named after the field 'name') will be created containing the same classification.


 - **multiple columns**  If multiple columns are specified, they will be combined into a new column (named after the field 'name') using the classification system defined below. This combined classification will then be compared against the experimental column

##### classification system

This section must be specified below the columns field. The script will apply one of the classification systems to the columns listed in columns, generating a combined classification under the new column specified in the name field. If only one column is specified, the new column will simply replicate its classification. Three classification schemes are available and are described below.

###### priority

The priority logic assigns a classification based on a defined order of precedence.

- class_order: a list specifying the order of classifications to consider. For each mutation, the script first checks whether any of the columns contain the first classification in the list.

- priority_for_classification  defines, for each classification, which columns have priority over others. The order is determined by the list. The script checks the first column in the priority list for the desired classification, then the next, and so on. If none of the specified columns contain the classification, the result is set to "Uncertain".

Example
```
priority:
        class_order: ["Damaging", "Neutral"]
        priority_for_classification:
          Damaging: ["GEMME classification", "Stability classification, alphafold, (RaSP, FoldX)"]
          Neutral:  ["GEMME classification", "Stability classification, alphafold, (RaSP, FoldX)"]
```

In this example, the script first looks for "Damaging" in "GEMME classification". If present, it assigns "Damaging". If not, it checks "Stability classification, alphafold, (RaSP, FoldX)". If not found, it then checks for "Neutral" in "GEMME classification" and subsequently in "Stability classification, alphafold, (RaSP, FoldX)".

###### weighted

The weighted logic works by assigning:

- weights: each column listed in columns must be assigned a weight. The sum of all weights must equal 1, otherwise an error is raised.A weighted score is calculated by summing the weights of all columns classified as "Damaging" and normalizing by the total weight.
It Returns
"Damaging" if normalized score ≥ threshold
"Neutral" if below threshold
"Uncertain" if total weight = 0

- threhsold: the cutoff value used to assign "Damaging".

###### voting

the voting logic works by assigning:

- target_class: the classification that this system will evaluate.

- logic: defines the rule for assigning the classification in target_class. Possible options are:
   - *majority*: if the majority of columns match the target_class, it is assigned.
   - *at_least_one*: if at least one column matches the target_class, it is assigned.
   - *all*: if all columns match the target_class, it is assigned.

 - fallback: the classification to assign if the logic defined above is not met.

 - handle_uncertain: specifies how to classify if all columns are "Uncertain". Options are:
 
   - *as_fallback*: return the classification defined in fallback
   - *ignore*: skip the assignment
   - *return_uncertain*: return "Uncertain".

N.B: 
- Multiple names (column combinations) can also be defined with different configurations. The script will process each accordingly, saving results in the corresponding folders.
- Multiple classification systems can be defined for the same column combination (name). For each system, the script will perform the appropriate analysis and save the results in dedicated folders.

 Here an example of config file:

```
modes:
  ensemble:
    experimental_column: "Experimental data (VAMP-seq, Stability (MAVEdb))"
    experimental_column_classification: "Experimental data classification (VAMP-seq, Stability (MAVEdb))"
    experimental_thresholds: [1, 2, 3, 4]
    normalization:
      Damaging:
        - "Damaging"
        - "damaging"
        - "Destabilizing"
        - "destabilizing"
        - "pathogenic"
        - "Pathogenic"
      Neutral:
        - "Neutral"
        - "neutral"
        - "Benign"
        - "benign"
      Uncertain:
        - "Uncertain"
        - "uncertain"
    comparisons:
      - name: "RaSP-FoldX_consensus_md"
        columns: ["Stability classification, (RaSP, FoldX) [md]"]
        priority:
          class_order: ["Damaging", "Neutral"]
          priority_for_classification:
            Damaging: ["Stability classification, (RaSP, FoldX) [md]"]
            Neutral: ["Stability classification, (RaSP, FoldX) [md]"]

  simple:
    experimental_column: "Experimental data (VAMP-seq, Stability (MAVEdb))"
    experimental_column_classification: "Experimental data classification (VAMP-seq, Stability (MAVEdb))"
    experimental_thresholds: [1, 2, 3, 4]
    normalization:
      Damaging:
        - "Damaging"
        - "damaging"
        - "Destabilizing"
        - "destabilizing"
        - "pathogenic"
        - "Pathogenic"
      Neutral:
        - "Neutral"
        - "neutral"
        - "Benign"
        - "benign"
      Uncertain:
        - "Uncertain"
        - "uncertain"
    comparisons:
      - name: "Rosetta-FoldX_consensus_with_GEMME"
        columns: ["Stability classification, alphafold, (Rosetta, FoldX)", "GEMME classification"]
        priority:
          class_order: ["Damaging", "Neutral"]
          priority_for_classification:
            Damaging: ["GEMME classification", "Stability classification, alphafold, (Rosetta, FoldX)"]
            Neutral: ["GEMME classification", "Stability classification, alphafold, (Rosetta, FoldX)"]
      - name: "Rosetta-FoldX_consensus_with_DeMaSK"
        columns: ["Stability classification, alphafold, (Rosetta, FoldX)", "DeMaSk classification"]
        priority:
          class_order: ["Damaging", "Neutral"]
          priority_for_classification:
            Damaging: ["DeMaSk classification", "Stability classification, alphafold, (Rosetta, FoldX)"]
            Neutral: ["DeMaSk classification", "Stability classification, alphafold, (Rosetta, FoldX)"]
      - name: "RaSP-FoldX_consensus_with_GEMME"
        columns: ["Stability classification, alphafold, (RaSP, FoldX)", "GEMME classification"]
        priority:
          class_order: ["Damaging", "Neutral"]
          priority_for_classification:
            Damaging: ["GEMME classification", "Stability classification, alphafold, (RaSP, FoldX)"]
            Neutral: ["GEMME classification", "Stability classification, alphafold, (RaSP, FoldX)"]
      - name: "RaSP-FoldX_consensus_with_DeMaSK"
        columns: ["Stability classification, alphafold, (RaSP, FoldX)", "DeMaSk classification"]
        priority:
          class_order: ["Damaging", "Neutral"]
          priority_for_classification:
            Damaging: ["DeMaSk classification", "Stability classification, alphafold, (RaSP, FoldX)"]
            Neutral: ["DeMaSk classification", "Stability classification, alphafold, (RaSP, FoldX)"]
```
The script must be run with the following flags:

**-S --simple** Selects the CSV files generated in MAVISp simple mode.

**-E --ensemble** Selects the CSV files generated in MAVISp ensemble mode.

**-y --config_yaml** Specifies the configuration file containing the parameters required to run the analysis (see below).

Optional flags can be provided to customize the analysis:

**-c --colormap** Specifies the colormap for histograms and confusion matrix plots. Default: viridis.

**-rg --remove_glycine** Excludes mutations to glycine when generating confusion matrices that compare predictions from MAVISp stability modules with experimental data. 

**-rp --remove_proline** Excludes mutations to proline when generating confusion matrices that compare predictions from MAVISp stability modules with experimental data.

**-sp --scatter_plot** Includes a scatter plot of the experimental data, highlighting the threshold used to define classifications.

**-rr --remove_residues** Excludes specific residues listed in a dedicated file, useful in cases where structural issues prevent reliable MAVISp calculations.


N.B The -sp flag can only be used if the experimental_column in the input file contains numerical scores from the experimental assay. Do not use this flag if the column contains categorical data (e.g., experimental classifications).

N.B If the experimental_column contains categorical (non-numerical) values, such as classifications derived from experimental results, scatter plots and distribution plots will not be generated, and the experimental_thresholds key will not be required. In this case, the -sp flag is unnecessary.

### File to exclude residues from the analysis
The -rr flag can be used to provide a csv file specifying residues that should be excluded from the entire analysis (MAVISp column comparisons). The file must contain the following columns:

**protein** The name of the protein, which must match the identifiers used in the MAVISp CSV file. This ensures that, when processing multiple files, the script can uniquely identify which mutations belong to each protein.

**position** The residue position to be excluded. All mutations affecting this position will be removed from the analysis.

Here is an example of a remove_residues file:

|protein|position|
|-------|--------|
|MLH1|E34|
|MLH1|C39|
|MLH1|N38|

## Output

The script generates an output directory called *output* containing one subfolder for each MAVISp mode (simple, ensemble, or both). Inside each subfolder, the following results are produced:

### Variant effect predicitors vs experimental data comparison

The output directory for variant effect predictors is structured as follow *classifier_comparison/confusion_matrix/* where: 
**classifier_comparison** is the classifier method under investigation (GEMME or DeMaSK)

This folder contains the 

- confusion_matrix_data.csv dataset used to make the confusion matrix

- confusion_matrix.pdf/png confusion matrix obtained by comparing the predictor of variant effects with the experimental data

- metrics.csv csv file containing the performances obtained from the comparison (see description section for a detailed explanation of the content)

Here an example of metrics.csv

Below is an example of the output files:

|sensitivity|specificity|accuracy|precision|F1 score| 
|-----------|-----------|--------|---------|--------|
|0.567|0.937|0.817|0.81|0.667|

If the -sp flag is activated, an additional output directory is generated with the following structure:
*classifier_comparison*/scatter_plots/mavisp_column/classification_system/ where:

**classifier_comparison** is the classifier method under investigation

**mavisp_column** is the MAVISp column specified in the config file under the name field

**classification_system** corresponds to the classification scheme specified in the config file

This folder contains:

- scatter_plot_data.csv. dataset used to make the scatter plot

- scatter_plot.pdf/png scatter plots comparing the predictor of variant effects with the experimental data


### MAVISp columns vs experimental data comparison

The output directory is structured as *mavisp_column/classification_system/*, where:

**mavisp_column** is the column specified in the config file under the *name* field (a folder is created for each entry in the config file).

**classification_system** is the classification scheme defined for that column (if multiple are specified, multiple folders with the appropriate analyses will be generated).

Each folder contains:

- confusion_matrix_data.csv dataset used to generate confusion matrices for MAVISp columns compared with experimental data

- confusion_matrix.pdf/png confusion matrix obtained by comparing the mavisp effects with the experimental data

- metrics.csv file containing the performances obtained from the comparison (see description section for a detailed explanation of the content)


### Variant effect predictiors on MAVISp dataset vs experimental data 

Inside each *mavisp_column/classification_system/* folder, there is a subfolder called:

*predictors_confusion_matrices_for_comparison/*
This contains two subfolders:
- GEMME
- DeMaSK
Each includes the confusion matrix and metrics of the corresponding variant effect predictor, applied to the dataset used for the MAVISp column analysis (*mavisp_column/classification_system/confusion_matrix_data.csv*).

This enables a direct comparison between MAVISp and the variant effect predictors on the same experimental dataset.
Additionally, a histogram plot is generated (histogram_plot.pdf and histogram_plot.png) summarizing the metrics from MAVISp and from GEMME/DeMaSK comparisons. This is stored in *mavisp_column/classification_system/*.

### analysis of divergent classifications

Inside each *mavisp_column/classification_system/* folder, there is also a subfolder (*divergent_classifications*) designed to collect analyses of mutations with divergent classifications between predictors and experimental data:

This directory contains two subfolders:

**discordant_neutral**  containing the analysis on mutations classified as Neutral by the predictor (MAVISp column) but Damaging by experimental data
**discordant_damaging** containing the analysis on mutations classified as Damaging by the predictor but Neutral by experimental data

Each folder contains the following:

- A CSV file listing the discordant mutations

- divergent_classification_accessibility.pdf/png  distribution of solvent accessibility values for discordant mutations

- divergent_classification_aminaocid_distribution.pdf/png percentage distribution of residues involved in the discordant dataset

- divergent_classification_aminoacid_class.pdf/png percentage distribution of amino acid chemical classes in the discordant dataset

Here an example of the output structure tree:

```
└── simple
    ├── DeMaSk_comparison
    │   ├── confusion_matrix
    │   │   ├── confusion_matrix_data.csv
    │   │   ├── confusion_matrix.pdf
    │   │   ├── confusion_matrix.png
    │   │   └── metrics.csv
    │   └── scatter_plots
    │       ├── RaSP-FoldX_consensus_with_DeMaSK_local
    │       │   ├── priority
    │       │   │   ├── scatter_plot_data.csv
    │       │   │   ├── scatter_plot.pdf
    │       │   │   └── scatter_plot.png
    │       │   ├── voting
    │       │   │   ├── scatter_plot_data.csv
    │       │   │   ├── scatter_plot.pdf
    │       │   │   └── scatter_plot.png
    │       │   └── weighted
    │       │       ├── scatter_plot_data.csv
    │       │       ├── scatter_plot.pdf
    │       │       └── scatter_plot.png
    │       └── RaSP-FoldX_consensus_with_GEMME_local
    │           └── priority
    │               ├── scatter_plot_data.csv
    │               ├── scatter_plot.pdf
    │               └── scatter_plot.png
    ├── GEMME_comparison
    │   ├── confusion_matrix
    │   │   ├── confusion_matrix_data.csv
    │   │   ├── confusion_matrix.pdf
    │   │   ├── confusion_matrix.png
    │   │   └── metrics.csv
    │   └── scatter_plots
    │       ├── RaSP-FoldX_consensus_with_DeMaSK_local
    │       │   ├── priority
    │       │   │   ├── scatter_plot_data.csv
    │       │   │   ├── scatter_plot.pdf
    │       │   │   └── scatter_plot.png
    │       │   ├── voting
    │       │   │   ├── scatter_plot_data.csv
    │       │   │   ├── scatter_plot.pdf
    │       │   │   └── scatter_plot.png
    │       │   └── weighted
    │       │       ├── scatter_plot_data.csv
    │       │       ├── scatter_plot.pdf
    │       │       └── scatter_plot.png
    │       └── RaSP-FoldX_consensus_with_GEMME_local
    │           └── priority
    │               ├── scatter_plot_data.csv
    │               ├── scatter_plot.pdf
    │               └── scatter_plot.png
    ├── RaSP-FoldX_consensus_with_DeMaSK_local_comparison
    │   ├── priority
    │   │   └── confusion_matrix
    │   │       ├── confusion_matrix_data.csv
    │   │       ├── confusion_matrix.pdf
    │   │       ├── confusion_matrix.png
    │   │       ├── distribution_plot.pdf
    │   │       ├── distribution_plot.png
    │   │       ├── divergent_classifications.csv
    │   │       ├── histogram_plot.pdf
    │   │       ├── histogram_plot.png
    │   │       ├── metrics.csv
    │   │       └── predictors_confusion_matrices_for_comparison
    │   │           ├── DeMaSk
    │   │           │   ├── confusion_matrix_data.csv
    │   │           │   ├── confusion_matrix.pdf
    │   │           │   ├── confusion_matrix.png
    │   │           │   └── metrics.csv
    │   │           └── GEMME
    │   │               ├── confusion_matrix_data.csv
    │   │               ├── confusion_matrix.pdf
    │   │               ├── confusion_matrix.png
    │   │               └── metrics.csv
    │   ├── voting
    │   │   └── confusion_matrix
    │   │       ├── confusion_matrix_data.csv
    │   │       ├── confusion_matrix.pdf
    │   │       ├── confusion_matrix.png
    │   │       ├── distribution_plot.pdf
    │   │       ├── distribution_plot.png
    │   │       ├── divergent_classifications.csv
    │   │       ├── histogram_plot.pdf
    │   │       ├── histogram_plot.png
    │   │       ├── metrics.csv
    │   │       └── predictors_confusion_matrices_for_comparison
    │   │           ├── DeMaSk
    │   │           │   ├── confusion_matrix_data.csv
    │   │           │   ├── confusion_matrix.pdf
    │   │           │   ├── confusion_matrix.png
    │   │           │   └── metrics.csv
    │   │           └── GEMME
    │   │               ├── confusion_matrix_data.csv
    │   │               ├── confusion_matrix.pdf
    │   │               ├── confusion_matrix.png
    │   │               └── metrics.csv
    │   └── weighted
    │       └── confusion_matrix
    │           ├── confusion_matrix_data.csv
    │           ├── confusion_matrix.pdf
    │           ├── confusion_matrix.png
    │           ├── distribution_plot.pdf
    │           ├── distribution_plot.png
    │           ├── divergent_classifications.csv
    │           ├── histogram_plot.pdf
    │           ├── histogram_plot.png
    │           ├── metrics.csv
    │           └── predictors_confusion_matrices_for_comparison
    │               ├── DeMaSk
    │               │   ├── confusion_matrix_data.csv
    │               │   ├── confusion_matrix.pdf
    │               │   ├── confusion_matrix.png
    │               │   └── metrics.csv
    │               └── GEMME
    │                   ├── confusion_matrix_data.csv
    │                   ├── confusion_matrix.pdf
    │                   ├── confusion_matrix.png
    │                   └── metrics.csv
    └── RaSP-FoldX_consensus_with_GEMME_local_comparison
        └── priority
            └── confusion_matrix
                ├── confusion_matrix_data.csv
                ├── confusion_matrix.pdf
                ├── confusion_matrix.png
                ├── distribution_plot.pdf
                ├── distribution_plot.png
                ├── divergent_classifications.csv
                ├── histogram_plot.pdf
                ├── histogram_plot.png
                ├── metrics.csv
                └── predictors_confusion_matrices_for_comparison
                    ├── DeMaSk
                    │   ├── confusion_matrix_data.csv
                    │   ├── confusion_matrix.pdf
                    │   ├── confusion_matrix.png
                    │   └── metrics.csv
                    └── GEMME
                        ├── confusion_matrix_data.csv
                        ├── confusion_matrix.pdf
                        ├── confusion_matrix.png
                        └── metrics.csv
```

## Usage

```
module load python/3.10/modulefile
python MAVISpBenchamrk.py -S -E -y config.yaml -sp -rr residue_to_exlude.csv
```

To run the analysis in the example folder:

```
bash run.sh
```