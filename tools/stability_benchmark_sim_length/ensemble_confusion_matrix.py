# This script is designed to calculate and visualise the performance of different md simulation times compared to a gold standard simulation time.

# Initialise
import pandas as pd
import csv
import os
import numpy as np
import seaborn as sn
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import argparse
from collections import defaultdict


# Compare simulations across simulation times
def comparison_data(path_components, gold_standard_simulation_time):

    gene_dict = defaultdict(list)
    for components in path_components:
        gene_name, structure_length, simulation_time = components

        gene_dict[gene_name].append((structure_length, simulation_time))


    results=[]

    for gene_name, components in gene_dict.items():
        golden_standard_components = None
        other_simulation_time=[]

        for structure_length, simulation_time in components:
            if simulation_time == gold_standard_simulation_time:
                golden_standard_components = (gene_name, structure_length, simulation_time)
            else:
                other_simulation_time.append((gene_name, structure_length, simulation_time))


        if golden_standard_components:
                
            comparison_results = {
                                       'gold_gene_name': gene_name,
                                       'structure_length': golden_standard_components[1],
                                       'gold_standard_simulation_time': gold_standard_simulation_time,
                                       'other_simulation_time': other_simulation_time
                                 }
            results.append(comparison_results)
        else:

            print(f"Golden standard time {gold_standard_simulation_time} not found for {gene_name}")

    return results




#Classify the ddG calculations (as classified in MAVISp protocol)
def classification(value):
    if value < -3:
        return "Stabilizing"
    elif value > 3:
        return "Destabilizing"
    elif -2 <= value <= 2:
        return "Neutral"
    else:
        return "Uncertain"
# Calculate components to use in confusion matrix
def confusion_matrix(standard_value, gold_standard_value, simulation_time, gold_standard_simulation_time, gene_name,gene_length):
    short_dest_long_dest = 0
    short_dest_long_neut = 0
    short_neut_long_dest = 0
    short_neut_long_neut = 0

# Find the number of true/false positives/negatives (positives = destabilizing)
    for index, rows in standard_value.iterrows():
        for cols in range(3, len(standard_value.columns)):
            standard = standard_value.iloc[index, cols]
            gold_standard = gold_standard_value.iloc[index, cols]

            if standard == "Destabilizing" and gold_standard == "Destabilizing":
                short_dest_long_dest += 1
            elif standard == "Destabilizing" and gold_standard == "Neutral":
                short_dest_long_neut += 1
            elif standard == "Neutral" and gold_standard == "Destabilizing":
                short_neut_long_dest += 1
            elif standard == "Neutral" and gold_standard == "Neutral":
                short_neut_long_neut += 1

    array = np.array([[short_neut_long_neut, short_neut_long_dest],
                      [short_dest_long_neut, short_dest_long_dest]])
    plot(array, simulation_time, gold_standard_simulation_time, gene_name,gene_length)
    metrics(array,simulation_time,gold_standard_simulation_time,gene_name,gene_length)
    return array

#Function to calculate performance metrics
def metrics(array,simulation_time,gold_standard_simulation_time,gene_name,gene_length):
    short_neut_long_neut = array[0][0]#true negatives
    short_neut_long_dest = array[0][1]#false negatives
    short_dest_long_neut = array [1][0]#false positives
    short_dest_long_dest = array [1][1]# true positives

# Calculate sensitivity TP/(FN+TP)
    if (short_neut_long_dest+short_dest_long_dest)>0:
        sensitivity = round(short_dest_long_dest/(short_neut_long_dest+short_dest_long_dest),3)
    else:
        sensitivity="nan"

#Calculate specificity TN/(TN+FP)
    if (short_neut_long_neut+short_dest_long_neut)>0:
        specificity = round(short_neut_long_neut/(short_neut_long_neut+short_dest_long_neut),3) #recall_score
    else:
        specificity="nan"

#Calculate precision TP/(TP+FP)
    if (short_dest_long_dest+short_dest_long_neut)>0:
        precision = round(short_dest_long_dest/(short_dest_long_dest+short_dest_long_neut),3)
    else:
        precision="nan"
#Calculate F1 score 2* TP/(2*TP+FP+FN)
    if (short_dest_long_dest+short_dest_long_neut+short_neut_long_dest)>0:
        if precision !="nan" and sensitivity != "nan":
            F1_score = round(2*short_dest_long_dest/(2*short_dest_long_dest+short_dest_long_neut+short_neut_long_dest),3)
        else:
            F1_score="nan"
    else:
        F1_score="nan"

#Calculate accuracy (TP+TN)/(TP+TN+FN+FP) 
    if short_dest_long_dest>0 or short_neut_long_dest>0 or short_neut_long_neut>0 or short_dest_long_neut>0:
        accuracy= round((short_dest_long_dest+short_neut_long_neut)/(short_dest_long_dest+short_neut_long_dest+short_neut_long_neut+short_dest_long_neut),3)
    else:
        accuracy="nan"
#Calculate MCC score (TP*TN-FP*FN)/(sqrt(TP+FP)*(TP+FN)*(TN+FP)*(TN*FN))
    if (short_dest_long_dest+short_dest_long_neut)>0 and (short_dest_long_dest+short_neut_long_dest)>0 and (short_neut_long_neut+short_dest_long_neut)>0 and (short_neut_long_neut+short_neut_long_dest)>0:
        MCC_score = round((short_dest_long_dest*short_neut_long_neut-short_dest_long_neut*short_neut_long_dest)/(np.sqrt((short_dest_long_dest+short_dest_long_neut)*(short_dest_long_dest+short_neut_long_dest)*(short_neut_long_neut+short_dest_long_neut)*(short_neut_long_neut+short_neut_long_dest))),3)
    else:
        MCC_score = "nan"

#Combine into a file

    metric={"sensitivity":sensitivity,
             "specificity":specificity,
             "accuracy":accuracy,
             "precision":precision,
             "F1 score":F1_score,
             "MCC score":MCC_score,
             "gene length":gene_length}
    file_title = f"{gene_name}_{simulation_time}_vs_{gold_standard_simulation_time}"
    plot_title = f"Stability comparison between {simulation_time}ns vs {gold_standard_simulation_time}ns of {gene_name}"

    protein_performance_metric = pd.DataFrame([metric])
    protein_performance_metric.to_csv(f"{output_dest}/performance_metrics_{file_title}.csv", index=False)

#Plot the confusion matrix and save as pdf and png file
def plot(array, simulation_time, gold_standard_simulation_time, gene_name,gene_length):
    Index = ["Neutral", "Destabilizing"]

    plot_title = f"Stability comparison between {simulation_time}ns vs {gold_standard_simulation_time}ns of {gene_name} of size {gene_length}"
    file_title = f"{gold_standard_simulation_time}ns_{simulation_time}ns_stability_comparison_{gene_name}"
    plt.figure(figsize=(12,8))
    sn.heatmap(array, annot=True, annot_kws={"size": 16}, fmt='g', xticklabels=Index, yticklabels=Index,cmap="Greens")
    plt.title(plot_title)
    plt.xlabel(f'{simulation_time}ns simulation', fontsize=15)
    plt.ylabel(f'{gold_standard_simulation_time}ns simulation', fontsize=15)

    plt.savefig(f"{output_dest}/{file_title}.png")
    plt.savefig(f"{output_dest}/{file_title}.pdf")
    plt.close()

input_directory=sys.argv[1]
gold_standard_simulation_time=sys.argv[2]
output_directory_name = sys.argv[3] if len(sys.argv) > 3 else "result"


current_directory = os.getcwd()
working_directory=os.path.join(current_directory,input_directory)
output_dest = os.path.join(current_directory,output_directory_name)
os.makedirs(output_dest, exist_ok=True)

data_dict = {}
path_components = []

for dirpath, dirnames, filenames in os.walk(working_directory):
    for filename in filenames:
        if filename.endswith('.csv'):
            path_folders = dirpath.split(os.sep)[-3:]
            if len(path_folders) == 3:
                 gene_name, structure_length, simulation_time = path_folders
            path_components.append((gene_name, structure_length, simulation_time))

            file_path = os.path.join(dirpath, filename)

            try:
                df = pd.read_csv(file_path)
                data_dict[file_path] = df
            except Exception as e:
                print(f"Could not read {file_path}: {e}")            

keys = list(data_dict.keys())
gold_standard_counter = 0
print("Classification started")
#Loop for classification
for path, df in data_dict.items():
    for index, row in df.iterrows():
        for col in df.columns[3:]:
            value = row[col]
            if type(value) == float or type(value) == int:
                df.at[index,col]=classification(value)


print("Classification done")
print("Confusion matrix calculation started")

#process the data 

compare_data = comparison_data(path_components, gold_standard_simulation_time)
confusion_matrix_data=[]
for i in range(len(data_dict)):
    for j in range(len(compare_data)):
        simulation_time = path_components[i][2]
        gene_name = path_components[i][0]
        gene_length = path_components[i][1]

        gold_standard_gene_name = compare_data[j]['gold_gene_name']
        gold_standard_sim_time = compare_data[j]['gold_standard_simulation_time']
        gold_gene_length = compare_data[j]['structure_length']


        if simulation_time != gold_standard_sim_time:
            if gene_name == gold_standard_gene_name :
                standard_value=data_dict[keys[i]]
                for k in keys:
                    keys_components=k.split('/')
                    if(keys_components[-2]==gold_standard_sim_time and keys_components[-4] == gold_standard_gene_name):
                        gold_standard_value = data_dict[k]
                        result=confusion_matrix(standard_value, gold_standard_value, simulation_time, gold_standard_simulation_time, gene_name,gene_length)

                        confusion_matrix_data.append((simulation_time,result))


summed_matrices = {}

for sim_time, matrix in confusion_matrix_data:
    if sim_time not in summed_matrices :
        summed_matrices[sim_time]=matrix
    else :
        summed_matrices[sim_time] += matrix

for sim_time, summed_matrix in summed_matrices.items():
    plot(summed_matrix,f"combined_{sim_time}", gold_standard_simulation_time , "all_genes","_all")
    metrics(summed_matrix, f"combined_{sim_time}" , gold_standard_simulation_time, "all_genes","_all")

print("Finished!")
