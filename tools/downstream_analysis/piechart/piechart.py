import pandas as pd
import re
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
import argcomplete
sns.set()

def piechart(df,out):

    # Select the "ClinVar Interpretation" column
    clinvar_interpretation = df['ClinVar Interpretation']
    # Count the occurrences of each value
    value_counts = clinvar_interpretation.value_counts()

    # Create a pie chart
    plt.figure(figsize=(12, 8))
    #colors = ['lightcoral', 'lightblue', 'lightgreen']  # Define your desired colors
    plt.pie(value_counts, startangle=140,autopct='')

    # Create a legend with percentages
    legend_labels = [f'{label} ({value_counts[label] / sum(value_counts) * 100:.1f}%)' for label in value_counts.index]
    plt.legend(legend_labels, title="Mavisp classification", loc="best",fontsize=18)
    plt.axis('equal')
    plt.tight_layout()
    # Show the pie chart
    plt.savefig(out,dpi=300)


if __name__ == "__main__":

    # Add arguments required to the script to a argparse.ArgumentParser instance.
    description = "Plot of the pie chart on the ClinVar intepretation " \
                  "column on mavisp aggregate csv. " 
    parser = argparse.ArgumentParser(description = description)

    i_helpstr = "Input: MAVISp aggregated csv."
    parser.add_argument("-i", "--input",
                        action = "store",
                        type = str,
                        help = i_helpstr,
                        required = True)

    d_helpstr = "Input: MAVISp dictionary file"
    parser.add_argument("-d", "--dictionary",
                        action = "store",
                        type = str,
                        help = d_helpstr,
                        required = True)    

    o_helpstr = "Input: MAVISp aggregated csv."
    parser.add_argument("-o", "--out_name",
                        action = "store",
                        type = str,
                        help = i_helpstr,
                        required = True) 
    
    argcomplete.autocomplete(parser)
    args = parser.parse_args()   
    
    # Read csv
    df = pd.read_csv(args.input)
    df["ClinVar Interpretation"] = df["ClinVar Interpretation"].str.replace(r",\s+", ",", regex=True)

    # Read  dictionary file and convert it to dictionary
    clinvar_dict = pd.read_csv(args.dictionary,sep="\t")
    clinvar_dict = clinvar_dict.set_index('#ClinVar').T.to_dict('records')
    clinvar_dict[0] = {k: (v.strip() if isinstance(v, str) else v) for k, v in clinvar_dict[0].items()}

    # Replace with Mavisp internal dictionary
    df["ClinVar Interpretation"] = df["ClinVar Interpretation"].replace(clinvar_dict[0])
    
    # Plot
    piechart(df,args.out_name)
