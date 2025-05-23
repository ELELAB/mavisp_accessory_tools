import pandas as pd
import csv
import os
import numpy as np
import seaborn as sn
import matplotlib.pyplot as plt
import argparse
import re
import glob
from sklearn.metrics import confusion_matrix
import sys

def plot_confusion_matrix(cm, 
                          sim_time, 
                          gold_sim_time, 
                          output_path, 
                          consensus,
                          gene = "whole_datasets", 
                          labels=["Destabilizing", "Neutral"]):
    '''
    Plot the confusion matrix using Seaborn.
    
    Parameters:
        cm (np.ndarray): The confusion matrix (2x2 numpy array).
        title (str): The title of the confusion matrix.
        labels (list): Labels for the classes (default is Positive and Negative).
    '''

    plt.figure(figsize=(6,6))
    sn.heatmap(cm, annot=True, fmt="d", cmap="Blues", cbar=False,
                xticklabels=[f"{label}" for label in labels],
                yticklabels=[f"{label}" for label in labels])
    
    title=f"{consensus} for {gene} ({gold_sim_time}_VS_{sim_time})"
    plt.title(title)
    plt.xlabel(sim_time)
    plt.ylabel(gold_sim_time)
    if gene:
        filename = f"{gene}_{gold_sim_time}_VS_{sim_time}"
    else:
        filename = f"{gold_sim_time}_VS_{sim_time}"
    plt.savefig(f"{output_path}/{filename}.png")
    plt.savefig(f"{output_path}/{filename}.pdf")
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


def metrics(array, 
            sim_time, 
            gold_sim_time, 
            output_path, 
            gene=False):
    '''
    computes various performance metrics based on the provided confusion
    matrix, comparing a "golden standard" simulation time with another 
    simulation time. The function generates sensitivity, specificity, 
    precision, F1 score, and accuracy and outputs these metrics to a 
    CSV file.

    The function extracts the following values from the confusion matrix: 
    true positives, true negatives, false positives, and false negatives. 
    It then calculates the following performance metrics:

        Sensitivity (Recall): The proportion of actual positives correctly 
                              identified.
        Specificity: The proportion of actual negatives correctly identified.
        Precision: The proportion of predicted positives that are true 
                   positives.
        F1 Score: The harmonic mean of precision and sensitivity.
        Accuracy: The overall proportion of correctly classified instances.

    Finally, it generates a CSV file containing the calculated metrics, 
    saved to the specified output_path. The output filename is 
    dynamically generated based on the provided simulation times and 
    the optional gene identifier.

    Parameters
    ----------
    array: np array: 
         A 2x2 array that represents the confusion matrix where:
         array[0][0]: True Positives (short-destabilized, long-destabilized)
         array[1][1]: True Negatives (short-neutral, long-neutral)
         array[0][1]: False Positives (short-neutral, long-destabilized)
         array[1][0]: False Negatives (short-destabilized, long-neutral)
    sim_time:str 
         The simulation time being evaluated.
    gold_sim_time:str 
         The "golden standard" simulation time for comparison.
    output_path:str 
         The directory where the output CSV file will be saved.
    gene(optional): str 
         A gene name to include in the output filename.

    Outputs:
    ----------
        CSV file: file to the output directory with the performance metrics.
    '''

    short_dest_long_dest = array[0][0]
    short_neut_long_neut = array [1][1]
    short_neut_long_dest = array[0][1]
    short_dest_long_neut = array [1][0]


    if (short_neut_long_dest+short_dest_long_dest)>0:
        sensitivity = round(short_dest_long_dest/(short_neut_long_dest+short_dest_long_dest),3)
    else:
        sensitivity="nan"
    if (short_neut_long_neut+short_dest_long_neut)>0:
        specificity = round(short_neut_long_neut/(short_neut_long_neut+short_dest_long_neut),3) #recall_score
    else:
        specificity="nan"
    if (short_dest_long_dest+short_dest_long_neut)>0:
        precision = round(short_dest_long_dest/(short_dest_long_dest+short_dest_long_neut),3)
    else:
        precision="nan"
    if precision !="nan" and sensitivity != "nan":
        if (precision+sensitivity)!=0:
            F1_score = round(2*(precision*sensitivity/(precision+sensitivity)),3)
        else:
            F1_score="nan"
    else:
        F1_score="nan"

    if short_dest_long_dest>0 or short_neut_long_dest>0 or short_neut_long_neut>0 or short_dest_long_neut>0:
        accuracy= round((short_dest_long_dest+short_neut_long_neut)/(short_dest_long_dest+short_neut_long_dest+short_neut_long_neut+short_dest_long_neut),3)
    else:
        accuracy="nan"

#Calculate MCC score (TP*TN-FP*FN)/(sqrt(TP+FP)*(TP+FN)*(TN+FP)*(TN*FN))
    if (short_dest_long_dest+short_dest_long_neut)>0 and (short_dest_long_dest+short_neut_long_dest)>0 and (short_neut_long_neut+short_dest_long_neut)>0 and (short_neut_long_neut+short_neut_long_dest)>0:
        MCC_score = round((short_dest_long_dest*short_neut_long_neut-short_dest_long_neut*short_neut_long_dest)/(np.sqrt((short_dest_long_dest+short_dest_long_neut)*(short_dest_long_dest+short_neut_long_dest)*(short_neut_long_neut+short_dest_long_neut)*(short_neut_long_neut+short_neut_long_dest))),3)
    else:
        MCC_score = "nan"

    metric={"sensitivity":sensitivity,
             "specificity":specificity,
             "accuracy":accuracy,
             "precision":precision,
             "F1 score":F1_score,
             "MCC score":MCC_score}
    if gene:
        filename = f"performance_{gene}_{gold_sim_time}_VS_{sim_time}"
    else:
        filename = f"performance_{gold_sim_time}_VS_{sim_time}"

    protein_performance_metric = pd.DataFrame([metric])
    protein_performance_metric.to_csv(f"{output_path}/{filename}.csv", index=False)


def file_consistency_check(csvs,sim_time):
    '''
    Filters CSV files based on their names and the specified simulation time.
    
    Parameters:
    csvs (list): List of CSV files to process.
    sim_time (str): Simulation time to look for in the file names.

    Returns:
    tuple: A list of filtered CSV files and a list of files to remove.
    '''


    genes = {}
    csvs_sets_to_remove =[]

    for csv in csvs:
        fname = os.path.basename(csv)
        # Match the file name using the regex
        m = re.match(r'(.*)_(\d+ns|clusters)-(ensemble_mode)\.csv', os.path.basename(csv))

        # Extract the protein name and simulation length
        protein = m.group(1)

        # Append extracted data to lists
        if protein not in genes.keys():
            genes[protein] = [csv]
        else:
            genes[protein].append(csv)

    for gene, csv_list in genes.items():
        # If no file in the list contains the simulation time, mark it for removal
        if not any(sim_time in csv for csv in csv_list):
            csvs_sets_to_remove.append(csv_list)
    
    # Count occurrences of each protein
    filtered_csv = [csv for csv in genes.values() if csv not in csvs_sets_to_remove]

    return filtered_csv,csvs_sets_to_remove

def process_dfs(df1, 
                df2, 
                column_name, 
                gold_sim_time, 
                sim_time):
    '''
    Processes two DataFrames by renaming a specified column, merging the 
    DataFrames on a common column, and filtering the merged DataFrame based 
    on specific criteria.

    This function performs the following steps:
    1. Renames the specified column in each DataFrame by appending 
       the provided simulation times.
    2. Merges the two DataFrames on the 'Mutation' column.
    3. Drops rows with missing values in the renamed columns.
    4. Filters out rows where the values in either of the renamed 
       columns are 'Uncertain' or 'Stabilizing'.

    Parameters
    ----------
        df1: pd.DataFrame 
             The first DataFrame to process.
        df2: pd.DataFrame
             The second DataFrame to process.
        column_name: str 
             The name of the column to rename in both DataFrames.
        gold_sim_time: str 
             The simulation time to append to the column name for df1.
        sim_time: str 
             The simulation time to append to the column name for df2.

    Returns
    -------
        tuple: A tuple containing:
            - merged_df (pd.DataFrame): The processed and 
              filtered DataFrame after merging.
            - column_mapping (dict): A dictionary mapping 
              'golden_standard' to the new column name in df1,
              and 'to_test' to the new column name in df2.
    '''

    new_column_name_df1 = column_name+"_"+gold_sim_time
    new_column_name_df2 = column_name+"_"+sim_time
    renamed_df1 = df1.rename(columns={column_name:new_column_name_df1})
    renamed_df2 = df2.rename(columns={column_name:new_column_name_df2})
    merged_df = pd.merge(renamed_df1, renamed_df2, on=['Mutation','protein'])
    merged_df = merged_df.dropna(subset=[new_column_name_df1, new_column_name_df2])
    merged_df = merged_df[~((merged_df[new_column_name_df1].isin(['Uncertain', 'Stabilizing'])) | 
                           (merged_df[new_column_name_df2].isin(['Uncertain', 'Stabilizing'])))]
    
    return merged_df, {'golden_standard':new_column_name_df1,'to_test':new_column_name_df2}



#############################################################################
#                                                                           #
#                                                                           #
#                                                                           #
#############################################################################
#                                 SCRIPT                                    #
#############################################################################
#                                                                           #
#                                                                           #
#                                                                           #
#############################################################################


parser=argparse.ArgumentParser(description = 'mavisp_sim_lenghts_comparison.py: \
                                              Processes MAVISP CSV data in ensemble mode, \
                                              containing stability calculations (FoldX/RaSP) \
                                              performed on ensembles extracted at different \
                                              time points. It evaluates the classification  \
                                              based on FoldX/RaSP consensus across datasets \
                                              generated from different ensembles obtained  \
                                              at varying simulation times, with one dataset\
                                              defined as the golden standard based on the \
                                              simulation time of interest.')


parser.add_argument("-t","--gold_sim_time",
                       dest = "gold_sim_time", 
                       required = True,
                       type = str,
                       help = "simulation time of the ensemble mode dataset to "
                              "use as golden standard expressed as {number}{ns}."
                              " For example 1000ns")
parser.add_argument("-i", "--input_dir",
                    dest="input_dir",
                    required=True,
                    type=str,
                    help="Directory containing the input CSV files. ")

args=parser.parse_args()



#############################################################################
#                                                                           #
#                             ERRROR HANDLING                               #
#                                                                           #
#############################################################################

# --------------------- Csv name and input flags check -------------------- #

# Check the argument

gold_sim_time = args.gold_sim_time

if not re.search(r'(\d+ns)', gold_sim_time):
    print("ERROR: the specified simulation time is in the wrong format."\
          " Please provide a valid simulation time:"\
          "{number}{ns} (for example 1000ns)")
    exit(1)

# Check the consistency between parameter and input files

available_csvs = glob.glob(os.path.join(args.input_dir, "*ensemble_mode.csv"))


if not any(gold_sim_time in csv for csv in available_csvs):
    print(f"None of the input files contain a simulation length compatible "\
          f"with the provided golden standard: {gold_sim_time}. Please "\
          f"provide a valid simulation time.")
    exit(1)


# ------------------------ Stability column check ------------------------- #

rosetta_stab_checked_csvs = []        # list containing the csv files with the 
                                      # correct (FoldX/RaSP) column
rasp_stab_checked_csvs = []

stab_missing_csvs = []        # list containing the csv files with the 
                              # wrong (FoldX/RaSP) column

rasp_pattern = r'(Stability classification, \(RaSP, FoldX\) \[md\])'
rosetta_pattern = r'(Stability classification, \(Rosetta, FoldX\) \[md\])'

for csv in available_csvs:
    fname = os.path.basename(csv)
    # Match the file name using the regex
    #m = re.match(r'(.*)_(\d+ns)-(ensemble_mode)\.csv', fname)
    m = re.match(r'(.*)_(\d+ns|clusters)-(ensemble_mode)\.csv', fname)
    if not m or len(m.groups()) != 3:
        print(f"The {csv} file name is not in the right format. Exiting...")
        exit(1)
    df = pd.read_csv(csv)

    rasp_columns = [col for col in df.columns if re.search(rasp_pattern, col)]
    rosetta_columns = [col for col in df.columns if re.search(rosetta_pattern, col)]

    if not rasp_columns and not rosetta_columns:
        stab_missing_csvs.append(csv)

    else:
        if rosetta_columns:
                rosetta_stab_checked_csvs.append(csv)
        if rasp_columns:
                rasp_stab_checked_csvs.append(csv)

#--------------------- Sorting the csvs file to ignore -------------------- #

single_csvs_to_remove = []
csvs_sets_to_remove = []

for df in stab_missing_csvs:
    if gold_sim_time in df:
        gene = df.split("_")[0]
        file_to_remove = [name for name in right_format_csvs if name.startswith(gene)]
        file_to_remove.append(df)
        csvs_sets_to_remove.append(set(file_to_remove))
    else:
        single_csvs_to_remove.append(df)


#--------------- Provide to the user the discarded files ------------------ #

# Remove those csv files missing both Rosetta/FoldX and RaSP/FoldX consenus

if single_csvs_to_remove:
    print(f"The following files will not be considered because missing "
          f"Stability classification (RaSP/FoldX) [md] and Stability classification (Rosetta/FoldX) [md] columns:'")
    for csv in single_csvs_to_remove:
        print(csv)

if csvs_sets_to_remove:
    print(f"The following sets of files will not be considered because the file"
          f" containing the golden standard simulation time misses"
          f" Stability classification (RaSP/FoldX) [md] and Stability classification (Rosetta/FoldX) [md] columns:'")

    for csv_set in csvs_sets_to_remove:
        print(csv_set)

# Remove those csvs files for which the golden standard csv is missing

filtered_rosetta_stab_checked_csvs = [csv for csv in rosetta_stab_checked_csvs if csv not in single_csvs_to_remove and csv not in [item for subset in csvs_sets_to_remove for item in list(subset)]]
filtered_rosetta_stab_checked_csvs,rosetta_unusable_csvs = file_consistency_check(filtered_rosetta_stab_checked_csvs,args.gold_sim_time)


filtered_rasp_stab_checked_csvs = [csv for csv in rasp_stab_checked_csvs if csv not in single_csvs_to_remove and csv not in [item for subset in csvs_sets_to_remove for item in list(subset)]]
filtered_rasp_stab_checked_csvs,rasp_unusable_csvs = file_consistency_check(filtered_rasp_stab_checked_csvs,args.gold_sim_time)


# Inform the user about the csvs files that will be used for the analysis
print("The following set of csvs files will be considered for the confusion matrix with (Rosetta/FoldX) consensus: ")

for csv in filtered_rosetta_stab_checked_csvs:
    print(set(csv))

print(3*"\n","The following set of csvs files will be considered for the confusion matrix with (RaSP/FoldX) consensus: ")

for csv in filtered_rasp_stab_checked_csvs:
    print(set(csv))


# Flat the csvs files list
filtered_rosetta_stab_checked_csvs = [item for sublist in filtered_rosetta_stab_checked_csvs for item in sublist]
filtered_rasp_stab_checked_csvs = [item for sublist in filtered_rasp_stab_checked_csvs for item in sublist]



consensus_filtered_stab_checked_csvs = {"Rosetta-FoldX_consensus":filtered_rosetta_stab_checked_csvs,
                                        "RaSP-FoldX_consensus":filtered_rasp_stab_checked_csvs}


#############################################################################
#                                                                           #
#                               INPUT PARSING                               #
#                                                                           #
#############################################################################
for consensus,filtered_stab_checked_csvs in consensus_filtered_stab_checked_csvs.items():

    whole_comparison_dfs = {}
    gene_specific_comparison_dfs = {}
    stability_classification = ""

    if str(consensus) == "RaSP-FoldX_consensus":
        stability_classification = "Stability classification, (RaSP, FoldX) [md]"
    elif str(consensus) == "Rosetta-FoldX_consensus":
        stability_classification =  "Stability classification, (Rosetta, FoldX) [md]"

    for csv in filtered_stab_checked_csvs:
        fname = os.path.basename(csv)
        m = re.match(r'(.*)_(\d+ns|clusters)-(ensemble_mode)\.csv', fname)
        df = pd.read_csv(csv)   

        protein = m.groups()[0]
        sim_time = m.groups()[1]

        df.insert(loc=0, column="protein", value=protein)

        df["Mut_res"] = df["Mutation"].str[-1]
        df["Position"] = df["Mutation"].str[1:-1]
        df["Position"] = df["Position"].astype(int)

        if sim_time not in whole_comparison_dfs.keys():
            whole_comparison_dfs[sim_time] = [df]
        else:
            whole_comparison_dfs[sim_time].append(df)

        if protein not in gene_specific_comparison_dfs.keys():
            gene_specific_comparison_dfs[protein] = {sim_time:df}
        else:
            gene_specific_comparison_dfs[protein].update({sim_time:df})


    #############################################################################
    #                                                                           #
    #                             OUTPUT PRODUCTION                             #
    #                                                                           #
    #############################################################################

    # ----------------------- gene specific comparison ------------------------ #

    for gene in gene_specific_comparison_dfs:
        for sim_time in gene_specific_comparison_dfs[gene]:
            if sim_time != gold_sim_time:
                output_path = f'output/{consensus}/gene_specific_comparison_{gold_sim_time}/{gene}/{sim_time}'
                golden_standard = gene_specific_comparison_dfs[gene][gold_sim_time]
                testing_df = gene_specific_comparison_dfs[gene][sim_time]


                merged_df,columns = process_dfs(golden_standard,
                                                testing_df,
                                                stability_classification,
                                                gold_sim_time,
                                                sim_time)

                os.makedirs(output_path, exist_ok=True)

                cm = generate_confusion_matrix(merged_df,
                                               columns['golden_standard'],
                                               columns['to_test'])
                plot_confusion_matrix(cm, 
                                      sim_time, 
                                      gold_sim_time, 
                                      output_path,
                                      consensus,
                                      gene=gene)
                metrics(cm, 
                        sim_time, 
                        gold_sim_time, 
                        output_path, 
                        gene=gene)

                merged_df.to_csv(f"{output_path}/{gene}_{gold_sim_time}_VS_{sim_time}_dataset.csv",sep=",",index=False)

    # -------------------------- whole comparison ----------------------------- #

    if len(gene_specific_comparison_dfs) > 1:
        golden_standard_dfs = pd.concat(whole_comparison_dfs[gold_sim_time])
        for sim_time in whole_comparison_dfs:
            if sim_time != gold_sim_time:
                output_path = f'output/{consensus}/overall_comparison_{gold_sim_time}/{sim_time}'
                golden_standard_dfs_copy = golden_standard_dfs.copy()
                testing_concat_dfs = pd.concat(whole_comparison_dfs[sim_time])
                genes_to_test = testing_concat_dfs['protein']
                golden_standard_dfs_copy = golden_standard_dfs_copy.loc[golden_standard_dfs_copy['protein'].isin(genes_to_test)]

                merged_df,columns = process_dfs(golden_standard_dfs_copy,
                                                testing_concat_dfs,
                                                stability_classification,
                                                gold_sim_time,
                                                sim_time)

                os.makedirs(output_path, exist_ok=True)
                cm = generate_confusion_matrix(merged_df,
                                               columns['golden_standard'],
                                               columns['to_test'])
                plot_confusion_matrix(cm, 
                                      sim_time, 
                                      gold_sim_time, 
                                      output_path,
                                      consensus)
                metrics(cm, 
                        sim_time, 
                        gold_sim_time, 
                        output_path)
                merged_df.to_csv(f"{output_path}/{gold_sim_time}_VS_{sim_time}_dataset.csv",sep=",",index=False)
    else:
        print(f"For the following {consensus} consensus all the input files are refferring to one gene. No overall comparison will be performed")






