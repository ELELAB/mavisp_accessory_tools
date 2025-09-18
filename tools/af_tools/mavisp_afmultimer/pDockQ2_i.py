'''
@karolinakrzesinska

This is an adaptation of the script sourced from: https://gitlab.com/ElofssonLab/afm-benchmark/
Changes made include parsing through all models and respective pkl files. 
Generating two files, one containing the pDockQ scores per chain and one with more detailed metrics.

Addition of function parse_and_score, process_models and changes to main function have been made. 

'''

# Imports 
from Bio.PDB import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Selection import unfold_entities
import numpy as np
import sys,os
import itertools
from scipy.optimize import curve_fit
from Bio.PDB import PDBParser
import numpy as np
import argparse
import pickle
from scipy.optimize import curve_fit
import pandas as pd
import json
import os
import logging as log


def retrieve_IFplddt(structure, chain1, chain2_lst, max_dist):
    ## generate a dict to save IF_res_id
    chain_lst = list(chain1) + chain2_lst

    ifplddt = []
    contact_chain_lst = []
    for res1 in structure[0][chain1]:
        for chain2 in chain2_lst:
            count = 0
            for res2 in structure[0][chain2]:
                if res1.has_id('CA') and res2.has_id('CA'):
                   dis = abs(res1['CA']-res2['CA'])
                   ## add criteria to filter out disorder res
                   if dis <= max_dist:
                      ifplddt.append(res1['CA'].get_bfactor())
                      count += 1

                elif res1.has_id('CB') and res2.has_id('CB'):
                   dis = abs(res1['CB']-res2['CB'])
                   if dis <= max_dist:
                      ifplddt.append(res1['CB'].get_bfactor())
                      count += 1
            if count > 0:
              contact_chain_lst.append(chain2)
    contact_chain_lst = sorted(list(set(contact_chain_lst)))   


    if len(ifplddt)>0:
       IF_plddt_avg = np.mean(ifplddt)
    else:
       IF_plddt_avg = 0

    return IF_plddt_avg, contact_chain_lst


def retrieve_IFPAEinter(structure, paeMat, contact_lst, max_dist):
    ## contact_lst:the chain list that have an interface with each chain. For eg, a tetramer with A,B,C,D chains and A/B A/C B/D C/D interfaces,
    ##             contact_lst would be [['B','C'],['A','D'],['A','D'],['B','C']]

 
    chain_lst = [x.id for x in structure[0]]
    seqlen = [len(x) for x in structure[0]]
    ifch1_col=[]
    ifch2_col=[]
    ch1_lst=[]
    ch2_lst=[]
    ifpae_avg = []
    d=10
    for ch1_idx in range(len(chain_lst)):
      ## extract x axis range from the PAE matrix
      idx = chain_lst.index(chain_lst[ch1_idx])
      ch1_sta=sum(seqlen[:idx])
      ch1_end=ch1_sta+seqlen[idx]
      ifpae_col = []   
      ## for each chain that shares an interface with chain1, retrieve the PAE matrix for the specific part.
      for contact_ch in contact_lst[ch1_idx]:
        index = chain_lst.index(contact_ch)
        ch_sta = sum(seqlen[:index])
        ch_end = ch_sta+seqlen[index]
        remain_paeMatrix = paeMat[ch1_sta:ch1_end,ch_sta:ch_end]
        #print(contact_ch, ch1_sta, ch1_end, ch_sta, ch_end)        

        ## get avg PAE values for the interfaces for chain 1
        mat_x = -1
        for res1 in structure[0][chain_lst[ch1_idx]]:
          mat_x += 1
          mat_y = -1
          for res2 in structure[0][contact_ch]:
              mat_y+=1
              if res1['CA'] - res2['CA'] <=max_dist:
                 ifpae_col.append(remain_paeMatrix[mat_x,mat_y])
      ## normalize by d(10A) first and then get the average
      if not ifpae_col:
        ifpae_avg.append(0)
      else:
        norm_if_interpae=np.mean(1/(1+(np.array(ifpae_col)/d)**2))
        ifpae_avg.append(norm_if_interpae)

    return ifpae_avg
    

def calc_pmidockq(ifpae_norm, ifplddt):
    df = pd.DataFrame()
    df['ifpae_norm'] = ifpae_norm
    df['ifplddt'] = ifplddt
    df['prot'] = df.ifpae_norm*df.ifplddt
    fitpopt = [1.31034849e+00, 8.47326239e+01, 7.47157696e-02, 5.01886443e-03] ## from orignal fit function  
    df['pmidockq'] = sigmoid(df.prot.values, *fitpopt)

    return df

def sigmoid(x, L ,x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0)))+b
    return (y)



def fit_newscore(df, column):

    testdf = df[df[column]>0]

    colval = testdf[column].values
    dockq = testdf.DockQ.values
    xdata =colval[np.argsort(colval)]
    ydata = dockq[np.argsort(dockq)]

    p0 = [max(ydata), np.median(xdata),1,min(ydata)] # this is an mandatory initial guess
    popt, pcov = curve_fit(sigmoid, xdata, ydata,p0)# method='dogbox', maxfev=50000)
    
#    tiny=1.e-20
#    print('L=',np.round(popt[0],3),'x0=',np.round(popt[1],3), 'k=',np.round(popt[2],3), 'b=',np.round(popt[3],3))

     ## plotting
#    x_pmiDockQ = testdf[column].values
#    x_pmiDockQ = x_pmiDockQ[np.argsort(x_pmiDockQ)]
#    y_pmiDockQ = sigmoid(x_pmiDockQ, *popt)
#    print("Average error for sigmoid fit is ", np.average(np.absolute(y_pmiDockQ-ydata)))

    
    #sns.kdeplot(data=df,x=column,y='DockQ',kde=True,levels=5,fill=True, alpha=0.8, cut=0)
#    sns.scatterplot(data=df,x=column,y='DockQ', hue='class')
#    plt.legend([],[], frameon=False)
    
#    plt.plot(x_pmiDockQ, y_pmiDockQ,label='fit',color='k',linewidth=2)
    return popt

def parse_and_score(pdb_path, pkl_path, dist):
    """Retrieve metrics from PKL and PDB files and calculate pDockQ score"""
    pdbp = PDBParser(QUIET=True)
    iopdb = PDBIO()

    structure = pdbp.get_structure('', pdb_path)
    chains = [chain.id for chain in structure[0]]
    # Define dictionaries 
    remain_contact_lst = []
    plddt_lst = []
    ptm_lst = []
    iptm_lst = []
    ranking_confidence_lst = []
    # retrieve pLDDT from pdb file
    for idx in range(len(chains)):
        chain2_lst = list(set(chains) - set(chains[idx]))
        IF_plddt, contact_lst = retrieve_IFplddt(structure, chains[idx], chain2_lst, dist)
        plddt_lst.append(IF_plddt)
        remain_contact_lst.append(contact_lst)
   # retrieve metrics and PAE from pkl 
    with open(pkl_path, 'rb') as f:
        data = pickle.load(f)

    avgif_pae = retrieve_IFPAEinter(structure, data['predicted_aligned_error'], remain_contact_lst, dist)
    ptm_lst.append(data.get('ptm', 0))  
    iptm_lst.append(data.get('iptm', 0))  
    ranking_confidence_lst.append(data.get('ranking_confidence', 0)) 
    # calculate pDockQ score
    res = calc_pmidockq(avgif_pae, plddt_lst)

    return chains, res['pmidockq'].tolist(), ptm_lst, iptm_lst, ranking_confidence_lst, avgif_pae, plddt_lst

def process_models(order, af_dir, suffix, models, dist):
    """For each model retrieve metrics and store in dict"""
    data_list = []
    # match pkl files to ranked pdbs
    for i, fname in enumerate(order[:models]):
        pdb_path = os.path.join(af_dir, f"ranked_{i}.pdb")
        pkl_path = os.path.join(af_dir, f"{suffix}_{fname}.pkl")

        chains, scores, ptm_scores, iptm_scores, ranking_confidence_scores, avgif_pae, plddt_lst = parse_and_score(pdb_path, pkl_path, dist)
        # construct dictionary per model
        model_data = {
            'ranked_i': i,
            'chains': chains,
            'scores': scores,
            'ptm_scores': ptm_scores,
            'iptm_scores': iptm_scores,
            'ranking_confidence_scores': ranking_confidence_scores,
            'interface_PAE' : avgif_pae,
            'interface_pLDDT' : plddt_lst
        }
        #Append to overall dictionary
        data_list.append(model_data)

    return data_list

def main():
    parser = argparse.ArgumentParser(description="Calculate chain-level pDockQ_i.")
    parser.add_argument('af_dir', help="Directory containing AlphaFold results")
    parser.add_argument('-r', '--ranking', help="JSON file containing the names of the ranked model", default="ranking_debug.json")
    parser.add_argument('-s', '--suffix', help="Suffix preceding the file model name", default="result")
    parser.add_argument('-n', '--models', default=None, type=int, help="Number of models to be considered in ranking order. Default is all of them")
    parser.add_argument("-dist", help="Maximum distance of a contact", nargs='?', type=int, default=8)

    args = parser.parse_args()

    # Parse the JSON file which includes ranking
    log.info(f"Parsing {args.ranking} in {args.af_dir}")
    try:
        with open(os.path.join(args.af_dir, args.ranking)) as fh:
            order = json.load(fh)['order']
    except IOError:
        log.error(f"Couldn't open {args.ranking} in the specified input directory; exiting...")
        exit(1)
    except json.decoder.JSONDecodeError:
        log.error(f"Couldn't parse {args.ranking}; is it in the right format? Exiting...")
        exit(1)
    except KeyError:
        log.error(f"{args.ranking} doesn't contain order information; Exiting...")
        exit(1)

    
    # Process models and store the relevant information
    model_data_list = process_models(order, args.af_dir, args.suffix, args.models, args.dist)

    # Write pDockQ scores to a file
    with open("pDockQ_scores.txt", "w") as output_file:
        # for each model processed 
        for model_data in model_data_list:
            i = model_data['ranked_i']
            scores = model_data['scores']
            output_file.write(f"Model ranked_{i} Scores:\n")
            # for each chain report score in output file
            for chain, score in zip(model_data['chains'], scores):
                output_file.write(f"{chain}: {score:.3f}\n")
            output_file.write("\n")

   # Write detailed scores to a CSV file
    with open("detailed_scores.csv", "w") as output_file:
        output_file.write("Ranked_Model,Chain,pdockq2,Interface_PAE,Interface_plDDT,PTM,iPTM,Ranking_Confidence\n")
        # for each model processed 
        for model_data in model_data_list:
            #define variables of interest 
            i = model_data['ranked_i']
            chains = model_data['chains']
            scores = model_data['scores']
            ptm_score = model_data['ptm_scores'][0]
            iptm_score = model_data['iptm_scores'][0] 
            ranking_confidence = model_data['ranking_confidence_scores'][0]  
            if_pae = model_data['interface_PAE']
            if_plddt = model_data['interface_pLDDT']
            # for each chain report all values in outputfile 
            for chain, score, ifpae, if_plddt in zip(chains, scores, if_pae, if_plddt):
                output_file.write(f"ranked_{i},{chain},{score:.3f},{ifpae:.3f},{if_plddt:.3f},{ptm_score:.3f},{iptm_score:.3f},{ranking_confidence:.3f}\n")

if __name__ == '__main__':
    main()
