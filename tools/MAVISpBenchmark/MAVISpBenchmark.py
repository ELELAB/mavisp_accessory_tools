import pandas as pd
import os
import re
import argparse
import numpy as np
import glob
import matplotlib.pyplot as plt
import yaml
from scipy.stats import pearsonr
import seaborn as sn
from sklearn.metrics import confusion_matrix
import numpy as np
import warnings
import matplotlib.patches as mpatches


def plot_confusion_matrix(cm, 
                          predictor_column, 
                          experiment_column, 
                          output_path,
                          color_map,
                          labels=["Damaging", "Neutral"]
                          ):
    '''
     Plot the confusion matrix using Seaborn.
     
     Parameters:
         cm (np.ndarray): The confusion matrix (2x2 numpy array).
         title (str): The title of the confusion matrix.
         labels (list): Labels for the classes (default is Positive and Negative).
    '''
    plt.figure(figsize=(10,6))
    sn.heatmap(cm, annot=True, fmt="d", cmap=color_map, cbar=False,
                 xticklabels=[f"{label}" for label in labels],
                 yticklabels=[f"{label}" for label in labels])

    title=f"({experiment_column}_VS_{predictor_column})"
    plt.title(title)
    plt.xlabel(predictor_column)
    plt.ylabel(experiment_column)

    for extension in [".png",".pdf"]:
        plt.savefig(f"{output_path}/confusion_matrix{extension}", bbox_inches='tight', dpi=300)

    plt.close()
     

    

def generate_confusion_matrix(df, 
                              gold_col, 
                              pred_col):
     '''
     Generates a confusion matrix by comparing two columns of a DataFrame.
     
     Parameters:
         df: pd.DataFrame 
           The DataFrame containing the columns to compare.
         gold_col: str 
           The name of the column containing the golden standard.
         pred_col: str 
           The name of the column containing the predictions.

     Outputs:
         np.ndarray: The confusion matrix.
     '''

     # Ensure both columns exist
     if gold_col not in df.columns or pred_col not in df.columns:
         raise ValueError(f"Columns '{gold_col}' and/or '{pred_col}' not found in the DataFrame.")

     # Extract the columns for the confusion matrix
     df[gold_col] = df[gold_col].astype(str)
     df[pred_col] = df[pred_col].astype(str)
     y_true = df[gold_col]
     y_pred = df[pred_col]

     # Compute the confusion matrix
     cm = confusion_matrix(y_true, y_pred)

     return cm


import pandas as pd

def metrics(array, 
            predictor_column, 
            experiment_column, 
            output_path):
    '''
    Computes various performance metrics based on the provided confusion
    matrix, comparing a "golden standard" simulation time with another 
    simulation time. The function generates sensitivity, specificity, 
    precision, F1 score, accuracy, MCC, and dataset size, and outputs these
    metrics to a CSV file.

    Parameters
    ----------
    array: np array 
        A 2x2 confusion matrix where:
        array[0][0]: True Positives (short-destabilized, long-destabilized)
        array[1][1]: True Negatives (short-neutral, long-neutral)
        array[0][1]: False Positives (short-neutral, long-destabilized)
        array[1][0]: False Negatives (short-destabilized, long-neutral)
    predictor_column: str 
        The simulation time being evaluated.
    experiment_column: str 
        The "golden standard" simulation time for comparison.
    output_path: str 
        Directory where the output CSV file will be saved.

    Outputs
    ----------
    CSV file: saved to the output directory with performance metrics.
    '''

    pred_dest_exp_dest = array[0][0]  # TP
    pred_neut_exp_neut = array[1][1]  # TN
    pred_neut_exp_dest = array[0][1]  # FP
    pred_dest_exp_neut = array[1][0]  # FN

    # Dataset size
    dataset_size = pred_dest_exp_dest + pred_neut_exp_neut + pred_neut_exp_dest + pred_dest_exp_neut

    # Sensitivity (Recall)
    if (pred_neut_exp_dest + pred_dest_exp_dest) > 0:
        sensitivity = round(pred_dest_exp_dest / (pred_neut_exp_dest + pred_dest_exp_dest), 3)
    else:
        sensitivity = "nan"

    # Specificity
    if (pred_neut_exp_neut + pred_dest_exp_neut) > 0:
        specificity = round(pred_neut_exp_neut / (pred_neut_exp_neut + pred_dest_exp_neut), 3)
    else:
        specificity = "nan"

    # Precision
    if (pred_dest_exp_dest + pred_dest_exp_neut) > 0:
        precision = round(pred_dest_exp_dest / (pred_dest_exp_dest + pred_dest_exp_neut), 3)
    else:
        precision = "nan"

    # F1 Score
    if precision != "nan" and sensitivity != "nan":
        if (precision + sensitivity) != 0:
            F1_score = round(2 * (precision * sensitivity / (precision + sensitivity)), 3)
        else:
            F1_score = "nan"
    else:
        F1_score = "nan"

    # Accuracy
    if dataset_size > 0:
        accuracy = round((pred_dest_exp_dest + pred_neut_exp_neut) / dataset_size, 3)
    else:
        accuracy = "nan"

    # MCC
    numerator = (pred_dest_exp_dest * pred_neut_exp_neut) - (pred_neut_exp_dest * pred_dest_exp_neut)
    denominator = ((pred_dest_exp_dest + pred_neut_exp_dest) *
                   (pred_dest_exp_dest + pred_dest_exp_neut) *
                   (pred_neut_exp_neut + pred_neut_exp_dest) *
                   (pred_neut_exp_neut + pred_dest_exp_neut)) ** 0.5
    if denominator > 0:
        mcc = round(numerator / denominator, 3)
    else:
        mcc = "nan"

    N_ref = GLOBAL_DATASET_SIZE  

    if mcc != "nan" and dataset_size > 0 and N_ref > 0:
        mcc_norm = round(mcc * np.sqrt(dataset_size / N_ref), 3)
    else:
        mcc_norm = "nan"

    # Collect metrics
    metric = {
        "sensitivity": sensitivity,
        "specificity": specificity,
        "accuracy": accuracy,
        "precision": precision,
        "F1 score": F1_score,
        "MCC": mcc,
        "MCC_norm": mcc_norm,
        "dataset_size": dataset_size
    }


    # Save to CSV
    protein_performance_metric = pd.DataFrame([metric])
    protein_performance_metric.to_csv(f"{output_path}/metrics.csv", index=False)
    
    return metric



def plot_residue_mutation_distribution(df, 
                                       mut_col, 
                                       wt_col,
                                       output_path):
    """
    Plots side-by-side barplots: percentage of mutated vs WT residues,
    coloring bars according to type, and also plots amino acid class distribution.

    Parameters:
    ----------
    df : pd.DataFrame
        DataFrame with mutation and WT residue columns.
    mut_col : str
        Column indicating mutated residue (e.g., 'A', 'R', etc.).
    wt_col : str
        Column indicating WT residue with position (e.g., 'A23').
    output_path : str
        Base path to save the plots.
    """
    # --- Preprocessing ---
    df["WT_res"] = df[wt_col].str.extract(r"([A-Z])")

    mut_counts = df[mut_col].value_counts(normalize=True) * 100
    wt_counts = df["WT_res"].value_counts(normalize=True) * 100

    residues = sorted(set(mut_counts.index) | set(wt_counts.index))

    plot_df = pd.DataFrame({
        "Residue": residues,
        "Mutated": [mut_counts.get(r, 0) for r in residues],
        "WT": [wt_counts.get(r, 0) for r in residues]
    }).melt(id_vars="Residue", value_vars=["Mutated", "WT"], var_name="Type", value_name="Percentage")

    # --- Plot histogram ---
    plt.figure(figsize=(12,6))
    sn.barplot(data=plot_df, x="Residue", y="Percentage", hue="Type")
    plt.title("Residue Mutation Distribution")
    plt.ylabel("Percentage (%)")
    plt.xlabel("Residue")
    plt.legend(title="Type", loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.savefig(f"{output_path}/divergent_classification_aminaocid_distribution.png", dpi=300)
    plt.savefig(f"{output_path}/divergent_classification_aminaocid_distribution.pdf")
    plt.close()

    # --- Amino acid class plot ---
    aa_classes = {
        "A": "Nonpolar", "V": "Nonpolar", "L": "Nonpolar", "I": "Nonpolar", "M": "Nonpolar",
        "F": "Aromatic", "Y": "Aromatic", "W": "Aromatic",
        "S": "Polar", "T": "Polar", "N": "Polar", "Q": "Polar", "C": "Polar",
        "D": "Acidic", "E": "Acidic",
        "K": "Basic", "R": "Basic", "H": "Basic",
        "G": "Nonpolar", "P": "Nonpolar"
    }

    df["Mut_class"] = df[mut_col].map(aa_classes)
    df["WT_class"] = df["WT_res"].map(aa_classes)

    mut_class_counts = df["Mut_class"].value_counts(normalize=True) * 100
    wt_class_counts = df["WT_class"].value_counts(normalize=True) * 100
    classes = sorted(set(mut_class_counts.index) | set(wt_class_counts.index))

    plot_df_class = pd.DataFrame({
        "Class": classes,
        "Mutated": [mut_class_counts.get(c,0) for c in classes],
        "WT": [wt_class_counts.get(c,0) for c in classes]
    }).melt(id_vars="Class", value_vars=["Mutated","WT"], var_name="Type", value_name="Percentage")

    plt.figure(figsize=(12,6))
    sn.barplot(data=plot_df_class, x="Class", y="Percentage", hue="Type")
    plt.title("Amino Acid Class Distribution")
    plt.ylabel("Percentage (%)")
    plt.xlabel("Class")

    # Custom legend handles for residue classes
    class_legend = [
        mpatches.Patch(color="white", label="Residue Classes:"),
        mpatches.Patch(color="white", label="Nonpolar: A, V, L, I, M, G, P"),
        mpatches.Patch(color="white", label="Aromatic: F, Y, W"),
        mpatches.Patch(color="white", label="Polar: S, T, N, Q, C"),
        mpatches.Patch(color="white", label="Acidic: D, E"),
        mpatches.Patch(color="white", label="Basic: K, R, H")
    ]

    # Legend for WT/Mutated + Residue classes
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles=handles + class_legend,
               loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)

    plt.tight_layout()
    plt.savefig(f"{output_path}/divergent_classification_aminoacid_class.png", dpi=300)
    plt.savefig(f"{output_path}/divergent_classification_aminoacid_class.pdf")
    plt.close()

def score_classification(df,column_name, new_column_name):
    """
    Adds a new classification column to the DataFrame based on the values of an existing column.

    Parameters:
    ----------
    df : pandas.DataFrame
        The original DataFrame.
    column_name : str
        The name of the column to use for classification.
    new_column_name : str, optional
        The name of the new classification column (default: 'classification').

    Returns:
    --------
    pandas.DataFrame
        The DataFrame with the new classification column added.
    """
    
    
    # Apply conditional function for the values classiifcation
    if "GEMME" in column_name:
        df[new_column_name] = df[column_name].apply(lambda x: 'Damaging' if x <= -3 or x >= 3 else 'Neutral')
    if "DeMaSk" in column_name:
        df[new_column_name] = df[column_name].apply(lambda x: 'Damaging' if x <= -0.25 or x >= 0.25 else 'Neutral')
    
    return df


def plot_accessibility_kde(df, 
                           col_name,
                           threshold,
                           output_path):
    """
    Plots the KDE distribution of side-chain solvent accessibility with a vertical threshold line.

    Parameters:
    ----------
    df : pd.DataFrame
        DataFrame containing accessibility values.
    col_name : str
        Column name with accessibility values.
    threshold : float
        Threshold to display as vertical line.
    output_path : str
        Base path for saving the figure (PDF and PNG).
    """
    if "divergent" in output_path:
        output_file_name = "divergent_classification_accessibility"
    else:
        output_file_name = "concordant_classification_accessibility"
    plt.figure(figsize=(12,6))
    sn.kdeplot(df[col_name], fill=True, color="blue", alpha=0.5)
    plt.axvline(threshold, color="red", linestyle="--", label=f"Threshold = {threshold}")
    plt.xlabel("Relative Side Chain Solvent Accessibility")
    plt.ylabel("Density")
    plt.title("Distribution of Side Chain Solvent Accessibility")
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.savefig(f"{output_path}/{output_file_name}.png", dpi=300)
    plt.savefig(f"{output_path}/{output_file_name}.pdf")
    plt.close()






def scatter_plot(df, variant_predictor_column, experiment_column, mavisp_column, output_path, color_mapping=None):
    """
    Create a scatter plot and a regression line from the provided data.

    Parameters
    ----------
    df: pandas DataFrame
        DataFrame containing the data
    variant_predictor_column: string
        Name of the column for the x-axis (predictions)
    experiment_column: string
        Name of the column for the y-axis (experiments)
    stability_column: string
        Name of the column for the color of the points
    output_path: string
        Path to store the output files
    color_mapping: dict, optional
        Dictionary that maps stability_column values to specific colors

    Returns
    -------
    scatter_plot: plt plot
        png and pdf files of the plots
    """

    # Create the scatter plot
    plt.figure(figsize=(12, 6))

    #df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=[prediction_column, experiment_column])

    df[mavisp_column] = df[mavisp_column].astype(str)

    # Check if color mapping is provided
    if color_mapping is None:
        color_mapping = {
            'Damaging': 'red',
            'Neutral': 'blue',
            'Uncertain':'grey',
            # Add more mappings as needed
        }

    # Map the colors using the provided or default color mapping
    if df[mavisp_column].dtype == 'object':
        colors = df[mavisp_column].map(color_mapping)

        scatter = plt.scatter(df[variant_predictor_column], df[experiment_column], c=colors, edgecolor='k', s=25)

        # Add a legend
        handles = [plt.Line2D([0], [0], marker='o', color=color_mapping[label], label=label, markersize=10) for label in color_mapping]
        plt.legend(handles=handles, title=mavisp_column)
    else:
        scatter = plt.scatter(df[variant_predictor_column], df[experiment_column], c=df[mavisp_column], cmap='viridis', edgecolor='k', s=25)
        plt.colorbar(scatter, label=mavisp_column)


    # Calculate Pearson correlation coefficient
    pearson_corr, _ = pearsonr(df[variant_predictor_column], df[experiment_column])
    
    # Fit a regression line
    slope, intercept = np.polyfit(df[variant_predictor_column], df[experiment_column], 1)
    regression_line = slope * df[variant_predictor_column] + intercept

    # Plot the regression line
    plt.plot(df[variant_predictor_column], regression_line, color='red', label='Regression Line', linewidth=2)

    title = f"{variant_predictor_column} vs {experiment_column}"

    # Legend adjustment
    legend_title = mavisp_column
    plt.legend(handles=handles + [plt.Line2D([0], [0], color='black', lw=0, marker='o', label=f'Pearson Correlation: {pearson_corr:.2f}')], 
               loc='upper left', bbox_to_anchor=(1.05, 1), title=legend_title, title_fontsize='small',fontsize='smaller', markerscale=0.3, labelspacing=0.3)

    plt.title(title)
    plt.xlabel(variant_predictor_column)
    plt.ylabel(experiment_column)
    plt.grid(False)  # Disable grid
    output = output_path + "/" + "scatter_plot"
    for extension in [".png", ".pdf"]:
        plt.savefig(f"{output}{extension}", bbox_inches='tight', dpi=300)
    plt.close()


def plot_distribution_with_highlighted_mutations(df, data_column, mavisp_classification_column, mutations_list, output_path, experimental_thresholds=None):
    '''
    Creates a distribution plot with KDE (kernel density estimation) using the data in the specified column.
    Highlights the points corresponding to the mutations provided in the list and colors the points based on classification.
    
    Parameters:
    ----------
    df: pd.DataFrame
        The DataFrame containing 'Mutation', 'classification', and a column with the data to plot.
    
    data_column: str
        The name of the column in the DataFrame containing the data to be plotted.
        
    mutations_list: list of str
        A list of mutations (strings) to be highlighted in the plot.
    '''
    # Convert to numeric if possible
    numeric_data = pd.to_numeric(df[data_column], errors="coerce")

    if numeric_data.isna().all():
        print(f"Skipping KDE plot for {data_column}: categorical data")
        return
    fig, ax = plt.subplots(figsize=(12, 6))
    # Create the distribution plot for all data
    sn.kdeplot(df[data_column],fill = True)
    
    # Filter the DataFrame based on the list of mutations
    filtered_df = df[df['Mutation'].isin(mutations_list)]

    # Map old labels to new labels
    label_mapping = {
        'Damaging': 'Damaging',
        'Neutral': 'Neutral'
    }

    filtered_df = filtered_df.copy()
    filtered_df['legend'] = filtered_df[mavisp_classification_column].map(label_mapping)

    # Count Neutral and Damaging
    counts = filtered_df['legend'].value_counts().to_dict()
    
    # Add the highlighted points (dot plot), colored by classification
    sn.scatterplot(x=filtered_df[data_column], 
                   y=[0.01] * len(filtered_df), 
                   hue=filtered_df['legend'],  # Color by classification
                   s=10,
                   ax=ax)

    if experimental_thresholds:
        for threshold in set(experimental_thresholds):
            ax.axvline(x=threshold, color='red', linestyle='--')

    # Add title and labels
    ax.set_title(f"Distribution Plot ({data_column})")
    ax.set_xlabel("Experimental scores")
    ax.set_ylabel("Frequency")

    # Custom legend with counts
    handles, labels = ax.get_legend_handles_labels()
    new_labels = []
    for lab in labels:
        if lab in counts:
            new_labels.append(f"{lab} (n={counts[lab]})")
        else:
            new_labels.append(lab)

    legend_title = mavisp_classification_column
    ax.legend(
        handles,
        new_labels,
        title=f'Discordant classifications {legend_title}',
        fontsize=5,
        title_fontsize=6,
        bbox_to_anchor=(1.05, 1),   # fuori a destra
        loc='upper left'
    )
        
    # save the plot
    output = output_path + "/" + "distribution_plot"
    for extension in [".png", ".pdf"]:
        fig.savefig(f"{output}{extension}", bbox_inches='tight', dpi=300)
    plt.close(fig)


def performance_histogram_plot(df, output_path, color_palette):
    """
    This function creates a grouped bar plot from a DataFrame showing performance metrics 
    (sensitivity, specificity, accuracy, F1 score, and precision) for different methods.
    Each method is represented by a different color, with bars grouped by metrics.
    
    Parameters:
    df (pd.DataFrame): A DataFrame containing performance data. The first column should be 'method',
                       and the remaining columns should be the metrics (sensitivity, specificity, etc.).
    output_path (str): The path where the plot should be saved.
    """

    plot_df = df.drop(columns=["dataset_size"], errors="ignore").copy()
    # Step 1: Ensure 'method' is the first column (reorder columns)
    columns_order = ['method'] + [col for col in plot_df.columns if col != 'method']
    plot_df = plot_df[columns_order]

    # Step 2: Ensure all metric columns are of float type
    for col in plot_df.columns[1:]:
        plot_df[col] = plot_df[col].astype(float)

    # Step 3: Define metrics and number of methods
    metrics = plot_df.columns[1:]
    num_metrics = len(metrics)
    num_methods = len(plot_df['method'].unique())

    # Step 4: Create a color palette using seaborn
    palette = sn.color_palette(color_palette, num_methods)  # 'Set2' palette, can be changed
    color_map = dict(zip(df['method'], palette))  # Map each method to a color from the palette

    # Step 5: Create a figure for the plot
    fig, ax = plt.subplots(figsize=(12, 6))

    # Step 6: Define the width of each bar and the positions for the groups of metrics
    bar_width = 0.2
    indices = np.arange(num_metrics)  # Positions for each metric

    # Step 7: Plot the bars for each method side by side within each metric group
    for i, method in enumerate(plot_df['method']):
        # Offset the positions of the bars for each method
        offset = i * bar_width
        ax.bar(indices + offset, plot_df.iloc[i, 1:], width=bar_width, label=method, color=color_map[method], alpha=0.7)

    # Step 8: Set the x-axis labels and ticks
    ax.set_xticks(indices + bar_width * (num_methods - 1) / 2)
    ax.set_xticklabels(metrics)

    # Step 9: Add labels and title
    ax.set_xlabel('Metrics')
    ax.set_ylabel('Scores')
    ax.set_title('Performance Metrics')

    # Step 10: Add the legend outside the plot
    ax.legend(
        title='Method',
        bbox_to_anchor=(1.05, 1),
        loc='upper left'           
    )

    # Step 11: Save the plot
    output = output_path + "/" + "histogram_plot"
    for extension in [".png", ".pdf"]:
        plt.savefig(f"{output}{extension}", bbox_inches='tight', dpi=300)
    
    # Step 12: Close the plot to free memory
    plt.close()


def df_for_histogram_plot_creation(performance_dictionary):
    dfs = []
    for method,metrics in performance_dictionary.items():
        df = pd.DataFrame([metrics])
        df["method"] = method
        dfs.append(df)

    return pd.concat(dfs)


def remove_mutations(df,column,mutation):
    filtered_df = df.loc[df[column] != mutation]
    return filtered_df


def filter_different_rows(df, col1, col2):
    """
    Filters rows in the DataFrame where the values in the two specified columns are different.

    Parameters
    ----------
    df: pd.DataFrame
        Input DataFrame containing the data.
    col1: str
        Name of the first column to compare.
    col2: str
        Name of the second column to compare.

    Returns
    -------
    pd.DataFrame
        A new DataFrame with three columns: 'mutation', col1, and col2.
        Only includes rows where the values in col1 and col2 are different.
    """

    # Check if the specified columns exist in the DataFrame
    if col1 not in df.columns or col2 not in df.columns or 'Mutation' not in df.columns:
        raise ValueError("One or more specified columns are not present in the DataFrame.")

    # Create a boolean mask for rows where values in col1 and col2 are different
    mask = df[col1] != df[col2]

    # Create the new DataFrame with the specified columns and mutation
    result_df = df.loc[mask]

    return result_df

# ---------------------------
# Normalization
# ---------------------------
def normalize_classification(value, normalization_dict):
    """
    Maps a value to its main class according to the normalization_dict.
    Case-insensitive. Returns the original value if no match is found.
    """
    if pd.isna(value):
        return value  # keep NaN as is
    str_val = str(value).strip().lower()
    for class_name, variants in normalization_dict.items():
        if any(str(v).strip().lower() == str_val for v in variants):
            return class_name
    warnings.warn(f"Value '{value}' not found in normalization dictionary.")
    return value

# ---------------------------
# Voting logic
# ---------------------------


def apply_voting_logic(row, columns, config):
    """
    Apply voting logic for a specific target class.
    """
    values = [row[col] for col in columns]
    target = config["target_class"]
    logic = config["logic"]
    fallback = config.get("fallback", "Unknown")
    handle_uncertain = config.get("handle_uncertain", "ignore")

    # Case: all are Uncertain
    if all(v == "Uncertain" for v in values):
        if handle_uncertain in ("ignore", "as_fallback"):
            return fallback
        elif handle_uncertain == "return_uncertain":
            return "Uncertain"

    # Apply voting strategy
    if logic == "majority":
        count_target = sum(1 for v in values if v == target)
        if count_target > len(values) / 2:
            return target
        elif count_target == len(values) / 2:
            return fallback
        else:
            return fallback

    elif logic == "at_least_one":
        return target if target in values else fallback

    elif logic == "all":
        return target if all(v == target for v in values) else fallback

    return fallback

# ---------------------------
# Priority logic
# ---------------------------
def apply_priority_logic(row, config):
    """
    Apply priority-based classification logic.
    
    Parameters:
        row (pd.Series): single row of DataFrame
        config (dict): contains 'class_order' and 'priority_for_classification'
        
    Returns:
        str: combined classification
    """
    class_order = config["class_order"]
    priority_for_class = config["priority_for_classification"]
    
    for target_class in class_order:
        if target_class not in priority_for_class:
            continue
        for col in priority_for_class[target_class]:
            if row[col] == target_class:
                return target_class
    return "Uncertain"

# ---------------------------
# Weighted logic
# ---------------------------
def apply_weighted_logic(row, columns, config):
    """
    Apply weighted classification.
    """
    weights = config.get("weights", {})
    threshold = config.get("threshold", 0.5)

    # Calculate weighted score
    score = 0.0
    total_weight = 0.0
    for col in columns:
        weight = weights.get(col, 1.0)
        total_weight += weight
        if row[col] == "Damaging":
            score += weight

    # Normalize and compare
    if total_weight == 0:
        return "Uncertain"
    normalized_score = score / total_weight
    return "Damaging" if normalized_score >= threshold else "Neutral"

# ---------------------------
# YAML validation
# ---------------------------
def validate_yaml(config, input_flags):
    """
    Validate YAML configuration for required fields and correct structure.

    input_flags: list of str, can contain "S" for simple, "E" for ensemble
    Returns: dict of {mode_name: mode_config} for all validated modes
    """

    # ------------------ Check modes key ------------------ #
    if "modes" not in config or not isinstance(config["modes"], dict):
        raise ValueError("YAML must contain a 'modes' dictionary with available configurations.")

    available_modes = config["modes"].keys()

    # Map input flags to mode names
    flag_to_mode = {"S": "simple", "E": "ensemble"}

    # Determine which modes to validate
    modes_to_validate = []
    for flag in input_flags:
        if flag not in flag_to_mode:
            raise ValueError(f"Unknown input flag '{flag}'")
        mode_name = flag_to_mode[flag]
        if mode_name not in available_modes:
            raise ValueError(f"Mode '{mode_name}' indicated by flag '{flag}' not found in YAML.")
        modes_to_validate.append(mode_name)

    # ------------------ Validate each selected mode ------------------ #
    for mode_name in modes_to_validate:
        cfg = config["modes"][mode_name]

        # --- Experimental section --- #
        required_experimental_fields = [
            "experimental_column",
            "experimental_column_classification",
            "experimental_thresholds"
        ]
        for field in required_experimental_fields:
            if field not in cfg:
                raise ValueError(f"YAML '{mode_name}' mode must contain '{field}' in experimental section.")

        thresholds = cfg.get("experimental_thresholds", [])
        if not isinstance(thresholds, list):
            raise ValueError(f"'{mode_name}' mode: experimental_thresholds must be a list of numeric values.")

        # --- Comparisons section --- #
        if "comparisons" not in cfg:
            raise ValueError(f"YAML '{mode_name}' mode must contain 'comparisons' section.")

        for comp in cfg["comparisons"]:
            if "name" not in comp:
                raise ValueError(f"Each comparison in '{mode_name}' must have a 'name'.")
            if "columns" not in comp or not isinstance(comp["columns"], list):
                raise ValueError(f"Comparison {comp.get('name')} must have a list of 'columns'.")

            # Validate voting logic
            if "voting" in comp:
                voting = comp["voting"]
                if "target_class" not in voting:
                    raise ValueError(f"Voting in {comp['name']} must have 'target_class'.")
                if "logic" not in voting or voting["logic"] not in ["majority", "at_least_one", "all"]:
                    raise ValueError(f"Voting logic in {comp['name']} must be one of 'majority', 'at_least_one', 'all'.")
                # Optional: check that fallback class exists in normalization
                fallback_class = voting.get("fallback")
                if fallback_class and fallback_class not in cfg.get("normalization", {}):
                    raise ValueError(f"Voting fallback '{fallback_class}' in {comp['name']} not defined in normalization.")

            # Validate priority logic
            if "priority" in comp:
                priority = comp["priority"]
                if "priority_for_classification" not in priority:
                    raise ValueError(f"Priority in {comp['name']} must have 'priority_for_classification'.")
                if not isinstance(priority["priority_for_classification"], dict):
                    raise ValueError(f"'priority_for_classification' in {comp['name']} must be a dictionary.")
                # Optional: check that all listed columns exist
                for cls, cols in priority["priority_for_classification"].items():
                    for col in cols:
                        if col not in comp["columns"]:
                            raise ValueError(f"Priority column '{col}' for class '{cls}' in {comp['name']} not listed in 'columns'.")

            # Validate weighted logic
            if "weighted" in comp:
                weighted = comp["weighted"]
                if "weights" not in weighted or not isinstance(weighted["weights"], dict):
                    raise ValueError(f"Weighted logic in {comp['name']} must have a 'weights' dictionary.")

                # Check che tutte le colonne nei pesi esistano nei columns
                for col in weighted["weights"]:
                    if col not in comp["columns"]:
                        raise ValueError(f"Weighted column '{col}' in {comp['name']} not listed in 'columns'.")

                # Check threshold range
                if "threshold" in weighted and not (0 <= weighted["threshold"] <= 1):
                    raise ValueError(f"Threshold in {comp['name']} weighted logic must be between 0 and 1.")

                # Check somma pesi = 1 (entro tolleranza numerica)
               # weight_sum = sum(weighted["weights"].values())
                #if abs(weight_sum - 1.0) > 1e-6:
                    #raise ValueError(
                       # f"In '{comp['name']}' weighted logic, the sum of weights must be 1. "
                       # f"Currently: {weight_sum}"
                    #)
# ---------------------------
# CSV columns validation
# ---------------------------

def validate_csv_columns(df, config):
    """
    Ensure all columns referenced in YAML exist in the dataframe.
    """
    missing_cols = set()

    # --- Check comparison columns --- #
    for comp in config["comparisons"]:
        for col in comp["columns"]:
            if col not in df.columns:
                missing_cols.add(col)

        # Check priority columns if present
        if "priority" in comp:
            for col_list in comp["priority"].get("priority_for_classification", {}).values():
                for col in col_list:
                    if col not in df.columns:
                        missing_cols.add(col)

        # Check weighted columns if present
        if "weighted" in comp:
            for col in comp["weighted"].get("weights", {}).keys():
                if col not in df.columns:
                    missing_cols.add(col)

    # --- Check experimental columns --- #
    if "experimental_column" in config:
        if config["experimental_column"] not in df.columns:
            missing_cols.add(config["experimental_column"])
    if "experimental_column_classification" in config:
        if config["experimental_column_classification"] not in df.columns:
            missing_cols.add(config["experimental_column_classification"])

    # Raise error if any columns are missing
    if missing_cols:
        raise ValueError(f"The following columns are missing in CSV: {', '.join(missing_cols)}")
    
    # --- Extra check: scatter flag requires numeric experimental column --- #
    if args and getattr(args, "scatter", False):
        exp_col = config.get("experimental_column")
        if exp_col and exp_col in df.columns:
            try:
                pd.to_numeric(df[exp_col], errors="raise")
            except Exception:
                raise ValueError(
                    f"--scatter flag specified but the experimental column '{exp_col}' "
                    f"contains non-numeric values (e.g., strings)."
                )

# ---------------------------
# Modified main function with validation
# ---------------------------

def process_input_files(config, selected_mode):
    """
    Load, validate, normalize, and classify all input CSV files according to YAML config.
    Returns a list of processed dataframes ready for downstream analysis.
    """

    processed_dfs = []

    
    available_csvs = glob.glob(f"*{selected_mode}_mode*.csv")
    if not available_csvs:
        print(f"No CSV files found for {selected_mode} mode. Exiting...")
        exit(1)

    for csv_file in available_csvs:
        fname = os.path.basename(csv_file)
        m = re.match(f"(.*)-({selected_mode}_mode.*)\\.csv", fname)
        if not m:
            print(f"File {csv_file} name is not in the expected format. Skipping...")
            continue

        protein_name = m.groups()[0]
        df = pd.read_csv(csv_file)
        df.insert(0, "protein", protein_name)

        # Extract mutation info
        df["Mut_res"] = df["Mutation"].str[-1]
        df["wt_res_position"] = df["Mutation"].str[:-1]
        df["Position"] = df["Mutation"].str[1:-1].astype(int)
        # Apply score classification for GEMME and DeMaSk if columns exist
        for column_name, new_column_name in zip(
            ["GEMME Score", "DeMaSk delta fitness"],
            ["GEMME classification", "DeMaSk classification"]
        ):
            if column_name in df.columns:
                df = score_classification(df, column_name, new_column_name)

        # Validate that all required columns exist
        try:
            validate_csv_columns(df, config)
        except ValueError as e:
            print(f"Skipping {csv_file}: {e}")
            continue

        # Apply normalization of classification columns
        def apply_normalization(df, cfg):
            """
            Applies classification normalization to columns listed in comparisons
            plus the experimental classification column.
            """
            normalization_dict = cfg.get("normalization", {})
            
            # Columns to normalize: all comparison columns + experimental classification
            cols_to_normalize = []
            for comp in cfg.get("comparisons", []):
                cols_to_normalize.extend(comp.get("columns", []))
            
            exp_col_class = cfg.get("experimental_column_classification")
            if exp_col_class:
                cols_to_normalize.append(exp_col_class)

            # Apply normalization only to columns present in df
            for col in cols_to_normalize:
                if col in df.columns:
                    df[col] = df[col].apply(lambda v: normalize_classification(v, normalization_dict))
            
            return df
        df = apply_normalization(df,config)

        # Apply comparison logics
        for comp in config["comparisons"]:
            name = comp["name"]
            cols = comp["columns"]

            if "voting" in comp:
                df[f"{name}_voting"] = df.apply(
                    lambda row: apply_voting_logic(row, cols, comp["voting"]), axis=1
                )
            if "priority" in comp:
                df[f"{name}_priority"] = df.apply(
                    lambda row: apply_priority_logic(row, comp["priority"]), axis=1
                )
            if "weighted" in comp:
                df[f"{name}_weighted"] = df.apply(
                    lambda row: apply_weighted_logic(row, cols, comp["weighted"]), axis=1
                )

        processed_dfs.append(df)

    # Merge all processed dataframes
    if processed_dfs:
        df_merged = pd.concat(processed_dfs, ignore_index=True)
    else:
        print("No valid CSVs processed. Exiting...")
        exit(1)


    return df_merged

def get_comparison_map(df, config):
    """
    Build a mapping between comparison names and the classification
    logics applied to them, along with the corresponding dataframe columns.

    Example return structure:
    {
        "Comparison1": {"voting": "Comparison1_voting", "priority": "Comparison1_priority"},
        "Comparison2": {"weighted": "Comparison2_weighted"}
    }
    """
    comparison_map = {}
    logic_suffixes = {
        "voting": "voting",
        "priority": "priority",
        "weighted": "weighted"
    }

    for comp in config["comparisons"]:
        comp_name = comp["name"]
        comparison_map[comp_name] = {}

        for logic, suffix in logic_suffixes.items():
            if logic in comp:  # logic specified in YAML
                col_name = f"{comp_name}_{suffix}"
                if col_name in df.columns:
                    comparison_map[comp_name][logic] = col_name
                else:
                    print(f"⚠️ Warning: expected column {col_name} not found in dataframe.")

    return comparison_map


def run_full_analysis(df_confusion_matrix, config, args, selected_mode,remove_positions=None,mutations_to_be_excluded = None):
    """
    Run the full analysis pipeline:
    1. Independent GEMME/DeMaSK analysis
    2. Filtered dataset analysis for YAML-specified columns
    3. Generate confusion matrices, scatter plots, and histogram plots
    """

    experiment_column = config["experimental_column"]
    experiment_column_classification = config["experimental_column_classification"]
    experimental_thresholds = config["experimental_thresholds"]

    comparison_map = get_comparison_map(df_confusion_matrix,config)


    df_confusion_matrix = df_confusion_matrix.dropna(subset=[experiment_column])
    df_confusion_matrix = df_confusion_matrix.loc[df_confusion_matrix[experiment_column_classification].isin(['Damaging', 'Neutral'])]

    # ------------------------------
    # GEMME and DeMaSK independent analysis
    # ------------------------------
    for predictior_column_score, predictor_column_classification in zip(
        ["GEMME Score", "DeMaSk delta fitness"],
        ["GEMME classification", "DeMaSk classification"]
    ):
        predictor = predictior_column_score.split(" ")[0]
        output_path = f"output/{selected_mode}/{predictor}_comparison/confusion_matrix"
        os.makedirs(output_path, exist_ok=True)

        df_pred = df_confusion_matrix.loc[
            df_confusion_matrix[predictor_column_classification].isin(["Damaging", "Neutral"])
        ]

        labels = ["Damaging", "Neutral"]
        cm = generate_confusion_matrix(df_pred, experiment_column_classification, predictor_column_classification)
        plot_confusion_matrix(cm, predictor_column_classification, experiment_column_classification, output_path, args.colormap, labels=labels)
        metrics(cm, predictor_column_classification, experiment_column, output_path)
        df_pred = df_pred[["Mutation",predictior_column_score,predictor_column_classification,experiment_column,experiment_column_classification]]
        df_pred.to_csv(f"{output_path}/confusion_matrix_data.csv",sep=",",index=False)

        if args.scatter:
            for comp, logics in comparison_map.items():
                for logic, comp_col in logics.items():
                    scatter_path = f"output/{selected_mode}/{predictor}_comparison/scatter_plots/{comp}/{logic}"
                    os.makedirs(scatter_path, exist_ok=True)
                    # ------------------- Confusion matrices and scatter plots -------------------
                    df_comp_col = df_confusion_matrix.loc[df_confusion_matrix[comp_col].isin(['Damaging','Neutral'])].copy()
                    # Confusion matrix vs experiment
                    df_comp_col.to_csv(f"{scatter_path}/scatter_plot_data.csv", sep=",", index=False)
                    scatter_plot(df_comp_col, predictior_column_score, experiment_column, comp_col, scatter_path)


    # ------------------------------
    # YAML-specified columns analysis
    # ------------------------------
    for comp, logics in comparison_map.items():
        for logic, comp_col in logics.items():
            comp_output = f"output/{selected_mode}/{comp}_comparison/{logic}"
            os.makedirs(comp_output, exist_ok=True)

            # ------------------------------
            # Filter dataset according to user specifications
            # ------------------------------

            # Remove positions specified in prot_pos_to_remove
            if args.residues_to_exclude:
                for p, pos in remove_positions:
                    if not ((df_confusion_matrix["protein"] == p) & (df_confusion_matrix["wt_res_position"] == pos)).any():
                        print(f"Position {pos} in {p} not found in DataFrame")
                mask = df_confusion_matrix.apply(lambda row: (row["protein"], row["wt_res_position"]) in remove_positions, axis=1)
                df_confusion_matrix = df_confusion_matrix[~mask]

            # ------------------- Confusion matrices and scatter plots -------------------
            data_for_histogram_plot = {}
            df_confusion_matrix_plot = df_confusion_matrix.copy()
            # Remove Glycine/Proline if requested
            for mut in mutations_to_be_excluded:
                df_confusion_matrix_plot = remove_mutations(df_confusion_matrix_plot, "Mut_res", mut)
            df_confusion_matrix_plot = df_confusion_matrix_plot.loc[df_confusion_matrix_plot[comp_col].isin(['Damaging','Neutral'])].copy()
            # Confusion matrix vs experiment
            cm_path = f"{comp_output}/confusion_matrix/"
            os.makedirs(cm_path, exist_ok=True)
            cm = generate_confusion_matrix(df_confusion_matrix_plot, experiment_column_classification, comp_col)
            plot_confusion_matrix(cm, comp_col, experiment_column_classification, cm_path, args.colormap)
            performances = metrics(cm, comp_col, experiment_column, cm_path)
            data_for_histogram_plot[comp_col] = performances
            divergent_classification = filter_different_rows(df_confusion_matrix_plot, 
                                                            comp_col, 
                                                            experiment_column_classification)
            # divergent classification analysis
            divergent_classification_analysis_path = f"{cm_path}/divergent_classification_analysis"
            divergent_classification_neutral_path = f"{divergent_classification_analysis_path}/discordant_neutral"
            divergent_classification_damaging_path = f"{divergent_classification_analysis_path}/discordant_damaging"
            os.makedirs(divergent_classification_neutral_path, exist_ok=True)
            os.makedirs(divergent_classification_damaging_path, exist_ok=True)
            divergent_classification.to_csv(f"{divergent_classification_analysis_path}/divergent_classification.csv",sep=",")
            plot_distribution_with_highlighted_mutations(df_confusion_matrix_plot, 
                                                        experiment_column, 
                                                        comp_col, 
                                                        divergent_classification["Mutation"].tolist(), 
                                                        divergent_classification_analysis_path, 
                                                        experimental_thresholds = experimental_thresholds)
            divergent_classification_neutral = divergent_classification.loc[divergent_classification[comp_col] == "Neutral"]
            divergent_classification_damging = divergent_classification.loc[divergent_classification[comp_col] == "Damaging"]
            # =============================
            # Concordant classification
            # =============================
            concordant_classification = df_confusion_matrix_plot.loc[
                df_confusion_matrix_plot[comp_col] ==
                df_confusion_matrix_plot[experiment_column_classification]
            ].copy()

            concordant_path = f"{cm_path}/concordant_classification_analysis"
            neutral_conc_path = f"{concordant_path}/concordant_neutral"
            damaging_conc_path = f"{concordant_path}/concordant_damaging"

            os.makedirs(neutral_conc_path, exist_ok=True)
            os.makedirs(damaging_conc_path, exist_ok=True)

            concordant_classification.to_csv(
                f"{concordant_path}/concordant_classification.csv", index=False
            )

            concordant_neutral = concordant_classification.loc[
                concordant_classification[comp_col] == "Neutral"
            ]

            concordant_damaging = concordant_classification.loc[
                concordant_classification[comp_col] == "Damaging"
            ]
            # =============================
            # Plot settings
            # =============================
            if selected_mode == "simple":
                accessibility_column = "Relative Side Chain Solvent Accessibility in wild-type"
            else:
                accessibility_column = "Relative Side Chain Solvent Accessibility in wild-type (average) [md]"
            # =============================
            # Plot ALL groups
            # =============================
            plot_groups = [
                (divergent_classification_neutral, divergent_classification_neutral_path, "discordant_neutral"),
                (divergent_classification_damging, divergent_classification_damaging_path, "discordant_damaging"),
                (concordant_neutral, neutral_conc_path, "concordant_neutral"),
                (concordant_damaging, damaging_conc_path, "concordant_damaging"),
            ]

            for df_group, output_path, label in plot_groups:

                if df_group.empty:
                    print(f"⚠️ Skipping {label}: no data")
                    continue

                plot_accessibility_kde(
                    df_group,
                    accessibility_column,
                    20,
                    output_path
                )

                plot_residue_mutation_distribution(
                    df_group,
                    "Mut_res",
                    "wt_res_position",
                    output_path
                )

                df_group.to_csv(
                    f"{output_path}/{label}.csv",
                    index=False
                )
            # save in csv format the data used for the confusion matrix
            df_confusion_matrix_data = df_confusion_matrix_plot[["Mutation",comp_col,experiment_column,experiment_column_classification]]
            df_confusion_matrix_data.to_csv(f"{cm_path}/confusion_matrix_data.csv",sep=",",index=False)

            # -- confusion matrixes of impact variant predictors for comparison with stability data -- #

            for predictor_column_score, predictor_column_classification in zip(["GEMME Score", "DeMaSk delta fitness"],["GEMME classification", "DeMaSk classification"]):
                predictor = predictor_column_score.split(" ")[0]
                comp_predictor_path = f"{cm_path}/predictors_confusion_matrices_for_comparison/{predictor}"
                os.makedirs(comp_predictor_path, exist_ok=True)
                cm = generate_confusion_matrix(df_confusion_matrix_plot,experiment_column_classification,predictor_column_classification)
                plot_confusion_matrix(cm,predictor_column_classification,experiment_column_classification,comp_predictor_path,args.colormap)
                performances = metrics(cm,predictor_column_classification,experiment_column,comp_predictor_path)
                data_for_histogram_plot[predictor] = performances

                df_confusion_matrix_data = df_confusion_matrix_plot[["Mutation",predictor_column_score,predictor_column_classification,experiment_column,experiment_column_classification]]
                df_confusion_matrix_data.to_csv(f"{comp_predictor_path}/confusion_matrix_data.csv",sep=",",index=False)

            df_for_histogram_plot = df_for_histogram_plot_creation(data_for_histogram_plot)
            performance_histogram_plot(df_for_histogram_plot,cm_path,args.colormap)


def main():
    global GLOBAL_DATASET_SIZE
    prot_pos_to_remove = []
    if args.residues_to_exclude:
        position_to_remove = pd.read_csv(args.residues_to_exclude)
        expected_cols = ["protein", "position"]
        file_cols = list(position_to_remove.columns)

        if file_cols != expected_cols:
            print(f"The file {args.residues_to_exclude} has columns {file_cols} "\
                f"instead of {expected_cols}. Please rename the columns to "\
                f" 'protein' and 'position'.")
            exit(1)
        
        # regex: lettera maiuscola seguita da numeri
        pattern = r"^[A-Z][0-9]+$"

        # crea colonna booleana valid/invalid
        position_to_remove["valid"] = position_to_remove["position"].str.match(pattern)

        # segnala errori
        errors = position_to_remove[~position_to_remove["valid"]]

        if not errors.empty:
            print(f"The file {args.residues_to_exclude} contain position\s in the wrong" \
                f" format:")
            print(errors[["protein", "position"]].to_string(index=False))
            exit(1)
        else:
            print(f"The file {args.residues_to_exclude} passed the quality check... "
                f"Proceding with the rest of the analysis")
        prot_pos_to_remove = list(zip(position_to_remove["protein"], position_to_remove["position"]))

    mutations_to_be_excluded = []
    for boolean,mutation in zip([args.no_glycine,args.no_proline],["G","P"]):
            if boolean:
                mutations_to_be_excluded.append(mutation)
    
    if not args.simple and not args.ensemble:
        print("Please specify for which mode (simple or ensemble) the analyis should be performed...")
        exit(1)

    # ---------------------------- Input files parsing ------------------------ #
    selected_modes = []
    input_flags = []
    if args.simple:
        selected_modes.append('simple')
        input_flags.append("S")
    if args.ensemble:
        selected_modes.append('ensemble')
        input_flags.append("E")

    for mode in selected_modes:
         # Load YAML configuration
        with open(args.yaml, 'r') as f:
            config = yaml.safe_load(f)
        validate_yaml(config,input_flags)
        config = config["modes"][mode]        

        # ---------------------- process input file --------------------------- #
        df_merged = process_input_files(config,mode)
        exp_col_class = config["experimental_column_classification"]
        df_filtered = df_merged.dropna(subset=[exp_col_class])
        df_filtered = df_filtered[df_filtered[exp_col_class].isin(['Damaging','Neutral'])]

        # -----------------------
        # Set global dataset size here
        # -----------------------
        GLOBAL_DATASET_SIZE = len(df_filtered)
        # ------------------------ plots the result --------------------------- #
        run_full_analysis(df_merged,config,args,mode,prot_pos_to_remove,mutations_to_be_excluded)

parser=argparse.ArgumentParser(description = 'exp_vs_pred_plot.py is a script designed for generating\
                                              scatter plots and confusion_matrixes. It visualizes prediction data \
                                              obtained from pathogenicity classifiers against \
                                              experimental data to assess the effects of variants and create confusion \
                                              matrixes comparing the Stability classification of MAVISp and \
                                              pathogenicty predictors (GEMME and DeMaSk) classifications with the experimental\
                                              data (considered as golden standard)')

parser.add_argument("-E","--ensemble",
                       dest="ensemble", 
                       default=False, 
                       action="store_true",
                        help="choose to plot only the data"\
                             " from the ensemble mode")

parser.add_argument("-S","--simple",
                       dest="simple", 
                       default=False, 
                       action="store_true",
                        help="choose to plot only the data"\
                             " from simple mode")

parser.add_argument("-y","--config_yaml",
                       dest="yaml", 
                       type=str,
                       required = True, 
                       help="config file with the columns to plot")


parser.add_argument("-rg","--remove_glycine",
                       dest="no_glycine", 
                       action="store_true",
                       required = False, 
                       default = False,
                       help="Choose whether to exclude mutations to "\
                            "glycine when generating confusion "\
                            "matrices for comparing predictions "\
                            "from MAVISp stability modules with "\
                            "experimental data.")

parser.add_argument("-rp","--remove_prolines",
                       dest="no_proline", 
                       action="store_true",
                       required = False, 
                       default = False,
                       help="Choose whether to exclude mutations to "\
                            "proline when generating confusion "\
                            "matrices for comparing predictions "\
                            "from MAVISp stability modules with "\
                            "experimental data.")

parser.add_argument("-rr","--remove_residues",
                       dest="residues_to_exclude", 
                       required = False, 
                       help="CSV file containing a list of positions "\
                            " (in P63 format) to be excluded; mutations " \
                            " mapped to these positions will be filtered out.")

parser.add_argument("-sp","--scatter_plot",
                       dest = "scatter", 
                       action = "store_true",
                       required = False, 
                       default = False,
                       help = "Include a scatter plot of the experimental data,"\
                              " in which the threshold that defines the "\
                              "classifications is highlighted.")


parser.add_argument("-c", "--colormap", 
                        dest="colormap", 
                        default="viridis",
                        help = f"choose the colormap for the histograms and confusion matrices plots: \
                             {', '.join(plt.colormaps())}. Default: viridis")


args=parser.parse_args()

if __name__ == "__main__":
    main()














