# MAVISpBenchmark.py
MAVISpBenchmark.py provides a comprehensive framework for analyzing and comparing MAVISp classifications and variant effect predictors against experimental data. The script evaluates one or more classification columns from MAVISp or external variant effect predictors (VEPs), individually or in combination. It supports both MAVISp simple-mode and ensemble-mode input files. The YAML configuration supports voting, priority and weighted combination strategies.

As output, the pipeline produces informative plots and performance metrics. In addition, it analyzes mutations with divergent classifications, highlighting discrepancies between predictor outputs and experimental reference data.

## Description
The script accepts one or more MAVISp CSV files in simple or ensemble mode and a YAML configuration specifying the experimental reference and classification strategies. Optional flags can be used to:

- Exclude mutations to glycine or proline.
- Choose a colormap for plots.
- Exclude specific residues through a separate CSV when structural calculations are unreliable.

The script verifies that a mode and configuration are provided and checks relevant input CSV columns and validates the optional exclusion file. The configuration file defines which columns should be analyzed and plotted. The input CSV must already contain a categorical experimental reference column that can be normalized to Damaging or Neutral. An experimental score column is also required by the current script; it may contain numerical assay scores or categorical assay readouts. `experimental_thresholds` only controls reference lines in relevant plots and does not generate experimental classes. The configuration further defines normalization rules, mapping predictor-specific outputs into standardized categories (Damaging, Neutral, or Uncertain). It also specifies which prediction columns should be compared—either single columns or combined sets—along with the classification scheme to be applied. A validation step ensures consistency between the configuration and the input CSV file.

After error handling, the script merges and organizes CSV files according to the selected mode (simple or ensemble) and begins processing. Columns with quantitative scores are classified using the thresholds specified in the YAML configuration under `score_classifications`.

The script normalizes configured prediction and experimental labels. Each strategy declared in `comparisons` then generates its own prediction column using the configured voting, priority or weighted logic, whether the strategy has one input column or several.

For each strategy explicitly listed under `comparisons`, the script evaluates variants with both a Damaging/Neutral experimental reference and a Damaging/Neutral prediction. It then generates performance metrics and, when enabled and applicable, plots and variant characterization outputs.

### Confusion Matrices
A confusion matrix is generated for **each classification strategy specified in `comparisons`**. This includes individual VEPs such as GEMME or DeMaSk when explicitly included in the YAML, MAVISp module classifications, and combinations of multiple prediction columns. The same procedure is used for one-column and multi-column comparisons; there is no special benchmark path for GEMME or DeMaSk. Each strategy is evaluated on the variants with both an experimental and a predicted Damaging/Neutral classification.

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

**Matthews Correlation Coefficient (MCC)**: a correlation-based measure of binary classification performance calculated from TP, TN, FP and FN.

**MCC_norm**: a coverage-normalized MCC, calculated as:

`MCC_norm = MCC × sqrt(N / N_ref)`

Here, `N` is the number of variants used in the individual comparison (after excluding non-binary predictions), and `N_ref` is the total number of variants with an evaluable experimental classification before filtering by the predictor. This metric makes the coverage of individual classifiers explicit.

**MCC_common**: an MCC calculated for each classification strategy on exactly the same set of variants. The script intersects the **final** `Damaging`/`Neutral` predictions from all strategies defined under `comparisons` in the YAML configuration, and retains only variants with a valid experimental classification. The MCC for each strategy is recalculated on this shared subset.

**N_common**: the number of variants in the common subset; it is identical for all the methods included in that particular run. The shared dataset is written to `output/<mode>/common_variants.csv`. Only classification strategies explicitly defined in `comparisons` determine this intersection. GEMME and DeMaSk are included in `MCC_common` only when their classification columns are also defined as comparisons.

The `metrics.csv` files include `MCC`, `MCC_norm`, `MCC_common`, `N_common`, and `dataset_size` in addition to sensitivity, specificity, accuracy, precision and F1 score. Changing the methods listed in `comparisons` may change both `N_common` and `MCC_common`.

### Scatter Plots
When plotting is enabled and the experimental scores are numerical, scatter plots are generated for each continuous score column declared under `score_classifications`, colored by the final classification of each strategy defined in `comparisons`. The x-axis is the raw numerical predictor score, the y-axis is the numerical experimental score, and colors indicate the final predicted class for variants included in that binary benchmark (Damaging or Neutral). For example, GEMME, DeMaSk, AlphaMissense and EVE can each be used as the x-axis if their numerical scores are declared in `score_classifications`. The plots include a linear regression and Pearson correlation. Missing/non-numeric scores are excluded **from scatter calculations, not from the underlying categorical benchmark**. Scatter plots are skipped when the experiment is categorical or too few valid numerical observations remain; confusion matrices and other benchmarking analyses are still computed. Each plot's underlying data are saved as `scatter_plot_data.csv`.

### Histograms of Performance
A performance comparison chart is generated once per mode (`performance_comparison.png`) from all methods configured in `comparisons`; `performance_summary.csv` contains the metrics of every strategy. The individual MCC metrics may be based on different numbers of variants. Use `MCC_common` when comparing performance on the exact same variants shared by all configured methods.

### Classification agreement and variant characterization
For each YAML classification strategy, the script separates variants into four groups according to predicted and experimental classifications:
- **True positive (TP):** predicted Damaging, experimentally Damaging.
- **True negative (TN):** predicted Neutral, experimentally Neutral.
- **False positive (FP):** predicted Damaging, experimentally Neutral.
- **False negative (FN):** predicted Neutral, experimentally Damaging.

TP and TN are stored under `classification_agreement/concordant/`, whereas FP and FN are stored under `classification_agreement/discordant/`. Each subgroup contains `variants.csv` listing the corresponding variants. When plots are enabled and there are sufficient data, the following visualizations are also generated:
- Histograms showing the wild-type and mutant amino acid composition.
- Histograms of wild-type and mutant amino acid chemical classes.
- Solvent accessibility distributions, when the appropriate structural feature is available and sufficient numerical data are present.

For the discordant variants, an additional distribution plot highlights incorrect predictions on the experimental score distribution, when numerical experimental values are available. These analyses are applied consistently to **all** strategies configured in the YAML, including one-column GEMME or DeMaSk comparisons.

## Requirements
Python >= 3.10

Packages:
- pandas
- numpy
- matplotlib
- seaborn
- scipy
- scikit-learn
- pyyaml

## Input
Input consists of one or more MAVISp CSVs for the selected mode and a YAML configuration. Input filenames must follow the `PROTEIN-simple_mode*.csv` or `PROTEIN-ensemble_mode*.csv` pattern, as appropriate; see **Config file** for how to specify the experimental reference, predictor columns and classification strategies. 

### Config file
The following schematic template shows the available fields (placeholders are not runnable values). Replace them with the exact CSV column names and chosen parameters:

```
modes:
  <mode_name>:               # e.g., "simple" or "ensemble"
    score_classifications:        # optional; convert any continuous score, e.g. GEMME, DeMaSk, EVE
      GEMME:
        column: <raw GEMME score column>
        output_column: "GEMME classification"
        damaging_threshold: -3.0
        direction: lower
        otherwise_class: Neutral
        missing_class: Uncertain
      DeMaSk:
        column: <raw DeMaSk score column>
        output_column: "DeMaSk classification"
        damaging_threshold: -0.25
        direction: lower
        otherwise_class: Neutral
        missing_class: Uncertain
    experimental_column: <string>
    experimental_column_classification: <string>
    experimental_thresholds: <list of numbers> # display thresholds only; use [] when none apply
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

The configuration fields are described below:

#### modes
This field can only take the values simple or ensemble, depending on the MAVISp CSV file of interest. If the analysis is to be performed on both simple and ensemble MAVISp CSV files, the configuration file must include two separate sections (one for simple and one for ensemble), each specifying the appropriate columns of interest for the corresponding file type.

#### score_classifications
This optional section converts raw numerical scores into categorical predictions before any strategy in `comparisons` is evaluated. It supports any named predictor (GEMME, DeMaSk, AlphaMissense, EVE, etc.) whose score column is available in the CSV. Required parameters are `column`, `output_column`, `direction` and the applicable threshold(s).

Supported directions and their Damaging conditions:

| `direction` | Damaging condition | Required cutoff(s) |

| --- | --- | --- |

| `lower` | score ≤ threshold | `damaging_threshold` |

| `higher` | score ≥ threshold | `damaging_threshold` |

| `inside_range` | lower ≤ score ≤ upper | `lower_threshold`, `upper_threshold` |

| `outside_range` | score ≤ lower or score ≥ upper | `lower_threshold`, `upper_threshold` |

`otherwise_class` (default `Neutral`) is used for valid numerical scores not satisfying the Damaging condition. `missing_class` (default `Uncertain`) is used for missing and non-numeric values. Missing or non-numeric scores are explicitly assigned `missing_class` rather than being left unclassified. Both fields can be set to a supported class (`Damaging`, `Neutral`, `Uncertain`); for binary experimental benchmarks the recommended settings are `otherwise_class: Neutral` and `missing_class: Uncertain`. Range boundaries are inclusive on the Damaging side for `outside_range` and inside the interval for `inside_range`.

Example:

```yaml
score_classifications:
  GEMME:
    column: "GEMME Score"
    output_column: "GEMME classification"
    direction: lower
    damaging_threshold: -3.0
    otherwise_class: Neutral
    missing_class: Uncertain
  DeMaSk:
    column: "DeMaSk delta fitness"
    output_column: "DeMaSk classification"
    direction: lower
    damaging_threshold: -0.25
    otherwise_class: Neutral
    missing_class: Uncertain
```

**Reproducibility:** The MSH2/MSH6 benchmark uses `direction: lower` and `damaging_threshold: -0.25` for DeMaSk, as shown above. An alternative `outside_range` rule can be configured with `lower_threshold` and `upper_threshold`, but it is **not equivalent** and may change the benchmark results.

Place `score_classifications` inside the selected `simple` and/or `ensemble` mode. Merely adding a score classifier does not benchmark it: to evaluate it, add its generated `output_column` to a comparison (with `name`, `columns` and `priority`, `voting` or `weighted`, even for one column). All methods follow this same rule.

#### experimental_column
The input CSV column containing raw experimental scores or readouts. It is required in the current configuration and used to select observations with a non-missing assay value; numerical values also enable scatter and experimental-score distribution plots.

#### experimental_column_classification
The input CSV column containing experimental classifications. These are normalized to Damaging/Neutral before benchmarking; the script does not derive them from `experimental_thresholds`.

#### experimental_thresholds
Optional numerical reference values displayed as vertical lines in the experimental-score distribution plot. **These thresholds do not classify experimental observations.** The `experimental_column_classification` column must already be present in the input CSV and contain categories that the normalization rules can map to Damaging or Neutral. Use `[]` when no thresholds apply.

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

The keys (Damaging, Neutral, Uncertain) represent standardized categories for the experimental classification and all prediction columns referenced under `comparisons`, including classes generated by `score_classifications`. No predictor receives special treatment.

The values listed under each key must include all possible labels that should be converted into that standardized category for the column of interest.

#### comparisons
This section lists every predictor or combined method to benchmark against the experimental data. No predictor is benchmarked automatically: GEMME and DeMaSk are included only if specified here.

Any number of comparisons can be configured in one run. Each comparison requires `name`, `columns` (one or more columns) and at least one strategy (`priority`, `voting` or `weighted`).

##### name 
A concise, descriptive identifier used for output paths and plot labels. Each configured strategy creates a prediction column named `<name>_<strategy>` (for example, `GEMME_voting`); no unqualified prediction column named only `<name>` is created.

##### columns
A list of one or more available classification columns, either already present in the MAVISp CSV or generated through `score_classifications`.

It can be used in two ways:
- **One column:** the declared `priority`, `voting` and/or `weighted` strategy is applied to that column. The result is saved in a prediction column named `<name>_<strategy>`; it is not automatically copied unchanged (for example, `Uncertain` inputs can be handled differently by voting or weighted strategies).
- **Multiple columns:** each configured strategy combines the selected columns into a prediction named `<name>_<strategy>`, then compares it with the experimental reference.

##### classification system
This section must be specified below the columns field. The script applies each classification strategy explicitly defined for a comparison to the columns listed under `columns`. This rule also applies when `columns` contains **only one column**: `priority`, `voting` and/or `weighted` must still be specified. There is **no implicit `direct` strategy**. Each strategy generates a separate prediction named `<name>_<strategy>` and an output folder under `output/<mode>/<name>/<strategy>/`. Multiple strategies may be defined on the same comparison.

###### priority
The priority logic assigns a classification based on a defined order of precedence.
- class_order: a list specifying the order of classifications to consider. For each mutation, the script first checks whether any of the columns contain the first classification in the list.
- priority_for_classification  defines, for each classification, which columns have priority over others. The order is determined by the list. The script checks the first column in the priority list for the desired classification, then the next, and so on. If none of the specified columns contain the classification, the result is set to "Uncertain".

Example

```
priority:
        class_order: ["Damaging", "Neutral"]
        priority_for_classification:
          Damaging: ["GEMME classification", "Stability classification, AFDB, (RaSP, FoldX)"]
          Neutral:  ["GEMME classification", "Stability classification, AFDB, (RaSP, FoldX)"]
```

In this example, the script first looks for "Damaging" in "GEMME classification". If present, it assigns "Damaging". If not, it checks "Stability classification, AFDB, (RaSP, FoldX)". If not found, it then checks for "Neutral" in "GEMME classification" and subsequently in "Stability classification, AFDB, (RaSP, FoldX)".

###### weighted
The weighted logic works by assigning:
- weights: each column listed in columns must be assigned a weight. Weights are relative values and do not need to sum to one. The weighted score is the sum of weights for columns classified as `Damaging`, divided by the sum of weights for **all** configured columns. `Neutral` and `Uncertain` contribute zero to the numerator but their weights remain in the denominator; missing predictions are not excluded from the denominator.

It returns:

"Damaging" if normalized score ≥ threshold

"Neutral" if below threshold

"Uncertain" if total weight = 0 (otherwise a binary class is returned)

- threshold: the cutoff value used to assign "Damaging".

###### voting
The voting logic assigns a class according to:
- target_class: the classification that this system will evaluate.
- logic: defines the rule for assigning the classification in target_class. Possible options are:

   - **majority**: if the majority of columns match the target_class, it is assigned.

   - **at_least_one**: if at least one column matches the target_class, it is assigned.

   - **all**: if all columns match the target_class, it is assigned.

- fallback: the classification to assign if the logic defined above is not met.

- handle_uncertain: specifies how to classify if all columns are "Uncertain". Options are:

   - **as_fallback**: return the configured `fallback` when all input columns are Uncertain

   - **ignore**: return the configured `fallback` when all input columns are Uncertain (identical to `as_fallback` for this case in the current implementation)

   - **return_uncertain**: return `Uncertain`.

Note: 
- Multiple names (column combinations) can also be defined with different configurations. The script will process each accordingly, saving results in the corresponding folders.
- Multiple classification systems can be defined for the same column combination (name). For each system, the script will perform the appropriate analysis and save the results in dedicated folders.

Example configuration:

```
modes:
  ensemble:
    score_classifications:
      GEMME:
        column: "GEMME Score"
        output_column: "GEMME classification"
        damaging_threshold: -3.0
        direction: lower
        otherwise_class: Neutral
        missing_class: Uncertain
      DeMaSk:
        column: "DeMaSk delta fitness"
        output_column: "DeMaSk classification"
        damaging_threshold: -0.25
        direction: lower
        otherwise_class: Neutral
        missing_class: Uncertain
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
        columns: ["Stability classification, (RaSP, FoldX) [md]_ensemble"]
        priority:
          class_order: ["Damaging", "Neutral"]
          priority_for_classification:
            Damaging: ["Stability classification, (RaSP, FoldX) [md]_ensemble"]
            Neutral: ["Stability classification, (RaSP, FoldX) [md]_ensemble"]
  simple:
    score_classifications:
      GEMME:
        column: "GEMME Score"
        output_column: "GEMME classification"
        damaging_threshold: -3.0
        direction: lower
        otherwise_class: Neutral
        missing_class: Uncertain
      DeMaSk:
        column: "DeMaSk delta fitness"
        output_column: "DeMaSk classification"
        damaging_threshold: -0.25
        direction: lower
        otherwise_class: Neutral
        missing_class: Uncertain
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
        columns: ["Stability classification, AFDB, (Rosetta, FoldX)", "GEMME classification"]
        priority:
          class_order: ["Damaging", "Neutral"]
          priority_for_classification:
            Damaging: ["GEMME classification", "Stability classification, AFDB, (Rosetta, FoldX)"]
            Neutral: ["GEMME classification", "Stability classification, AFDB, (Rosetta, FoldX)"]
      - name: "Rosetta-FoldX_consensus_with_DeMaSk"
        columns: ["Stability classification, AFDB, (Rosetta, FoldX)", "DeMaSk classification"]
        priority:
          class_order: ["Damaging", "Neutral"]
          priority_for_classification:
            Damaging: ["DeMaSk classification", "Stability classification, AFDB, (Rosetta, FoldX)"]
            Neutral: ["DeMaSk classification", "Stability classification, AFDB, (Rosetta, FoldX)"]
      - name: "RaSP-FoldX_consensus_with_GEMME"
        columns: ["Stability classification, AFDB, (RaSP, FoldX)", "GEMME classification"]
        priority:
          class_order: ["Damaging", "Neutral"]
          priority_for_classification:
            Damaging: ["GEMME classification", "Stability classification, AFDB, (RaSP, FoldX)"]
            Neutral: ["GEMME classification", "Stability classification, AFDB, (RaSP, FoldX)"]
      - name: "RaSP-FoldX_consensus_with_DeMaSk"
        columns: ["Stability classification, AFDB, (RaSP, FoldX)", "DeMaSk classification"]
        priority:
          class_order: ["Damaging", "Neutral"]
          priority_for_classification:
            Damaging: ["DeMaSk classification", "Stability classification, AFDB, (RaSP, FoldX)"]
            Neutral: ["DeMaSk classification", "Stability classification, AFDB, (RaSP, FoldX)"]
```

The script must be run with the following flags:

**-S --simple** Selects the CSV files generated in MAVISp simple mode.

**-E --ensemble** Selects the CSV files generated in MAVISp ensemble mode.

**-y --config_yaml** Specifies the configuration file containing the parameters required to run the analysis (see below).

Optional flags can be provided to customize the analysis:

**-c --colormap** Specifies the colormap for histograms and confusion matrix plots. Default: viridis.

**-rg --remove_glycine** Excludes mutations to glycine when generating confusion matrices that compare predictions from MAVISp stability modules with experimental data. 

**-rp --remove_proline** Excludes mutations to proline when generating confusion matrices that compare predictions from MAVISp stability modules with experimental data.

**-sp --scatter_plot** Retained for compatibility with older command lines. Scatter plots are now generated automatically for continuous scores listed under `score_classifications` when numerical experimental values are available and plotting is enabled.

**-rr --remove_residues** Excludes specific residues listed in a dedicated file, useful in cases where structural issues prevent reliable MAVISp calculations.

Note: Scatter plots need numeric predictor and experimental scores. If the experimental column is categorical, the scatter plot is skipped; other analyses still run.

Note: For categorical experimental data, numerical scatter/distribution plots are not produced. The experimental classification must still contain Damaging/Neutral labels for binary benchmarking.

### File to exclude residues from the analysis
The -rr flag can be used to provide a csv file specifying residues that should be excluded from the entire analysis (MAVISp column comparisons). The file must contain the following columns:

**protein** The name of the protein, which must match the identifiers used in the MAVISp CSV file. This ensures that, when processing multiple files, the script can uniquely identify which mutations belong to each protein.

**position** The residue position to be excluded. All mutations affecting this position will be removed from the analysis.

Here is an example of a remove_residues file:

| protein | position |
|---|---|
| MLH1 | E34 |
| MLH1 | C39 |
| MLH1 | N38 |

## Output
The script writes results under `output/<mode>/` where `<mode>` is `simple` or `ensemble`. It uses an identical folder layout for all predictors and comparison strategies; GEMME, DeMaSk and MAVISp combinations are **not** treated differently.

### Per-mode summary and common-variant data
- `common_variants.csv` contains the experimental class and the final predictions from **every configured classification strategy**, restricted to variants with Damaging/Neutral experimental and predicted classes across all methods. This is the dataset used to compute `MCC_common` and `N_common`.
- `performance_summary.csv` contains one row per comparison/strategy, with sensitivity, specificity, accuracy, precision, F1, MCC, `MCC_norm`, `MCC_common`, `N_common` and `dataset_size`.
- `performance_comparison.png` visually compares performance metrics across the configured strategies (omitted with `--no-plots`).

### Per-classifier and per-strategy outputs
For each comparison under the YAML `comparisons` list, and for each configured strategy (`priority`, `voting`, `weighted`), the script creates `output/<mode>/<name>/<strategy>/` with:
- `metrics.csv`: metrics for this particular classification strategy.
- `confusion_matrix/`: `confusion_matrix_data.csv`, `confusion_matrix.png` and `confusion_matrix.pdf` (figures omitted with `--no-plots`).
- `scatter_plot/<score_classifier>/` (when continuous scores and numeric experimental data are available): `scatter_plot_data.csv`, `scatter_plot.png`, `scatter_plot.pdf` for each score defined in `score_classifications` when both predictor and experimental score columns contain sufficient numerical data. Points are colored by the final class of the current strategy.
- `classification_agreement/concordant/true_positive/` and `.../true_negative/`: lists and molecular characterization of correctly classified variants.
- `classification_agreement/discordant/false_positive/` and `.../false_negative/`: lists and molecular characterization of incorrectly classified variants. The discordant directory also holds experimental-score distribution figures when possible.

Every TP/TN/FP/FN subfolder includes `variants.csv` (even if the group is empty). When sufficient data are available, molecular characterization figures are saved as PNG and PDF. The variant lists and metrics are always generated, including when plotting is disabled.

A typical output tree is:

```text
output/
└── simple/
    ├── common_variants.csv
    ├── performance_summary.csv
    ├── performance_comparison.png
    ├── GEMME/
    │   └── priority/
    │       ├── metrics.csv
    │       ├── confusion_matrix/
    │       ├── scatter_plot/
    │       │   ├── GEMME/
    │       │   └── DeMaSk/
    │       └── classification_agreement/
    │           ├── concordant/
    │           │   ├── true_positive/
    │           │   └── true_negative/
    │           └── discordant/
    │               ├── false_positive/
    │               └── false_negative/
    ├── DeMaSk/
    │   └── priority/
    │       └── ...
    └── MAVISp_combined/
        ├── priority/
        ├── voting/
        └── weighted/
```

`performance_comparison.png` and `scatter_plot/` are not generated when plots are disabled or their required numeric inputs are unavailable. The other result folders follow the YAML configuration rather than any hard-coded list of predictor names.

## Usage
```
module load python/3.10/modulefile
python MAVISpBenchmark.py -S -E -y config.yaml -rr residue_to_exclude.csv
```

Use `--no-plots` to generate all benchmark CSVs (including `performance_summary.csv`, `common_variants.csv`, individual metrics, confusion matrix data and TP/TN/FP/FN variant tables) without rendering images. This is useful for fast regression tests. Numerical scatter plots are generated automatically from `score_classifications`; categorical/insufficient experimental data cause the plots to be skipped.

To run the analysis in the example folder:

```
bash run.sh
```
