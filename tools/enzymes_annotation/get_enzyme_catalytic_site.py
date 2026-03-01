#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (C) 2026, Matteo Arnaudi  <mata@cancer.dk>,

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import requests
import pandas
import json
import argparse

def from_uniprot_ac_to_gene_name(uniprot_ac):

    '''
    Retrieve the gene name associated to the input uniprot ac

    The function takes as input an uniprot_ac and queries the
    Uniprot database to return the correpsonding Gene name asso
    ciated to the uniprot_ac provided as input.
    If the uniprot_ac does not contain any gene_name information
    it will be returned as no_gene

    Parameters
    ----------
    uniprot_ac: string
        Uniprot_ac of the protein of interest

    Returns
    -------
    genes: string
        genes name associated to the uniprot_ac
    no_genes: string
        uniprot_ac provided as input in case there are no gene names
        associated to it
    '''

    uniprot_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_ac}.json"
    uniprot_response = requests.get(uniprot_url)
    genes_name = []
    no_genes_name = ""
    if uniprot_response.status_code == 200:
        uniprot_entry_information = uniprot_response.json()
        if "genes" in uniprot_entry_information.keys():
            try:
                for gene in uniprot_entry_information['genes']:
                    if not 'geneName' in gene.keys():
                        no_genes_name = uniprot_ac
                        continue
                    genes_name.append(gene['geneName']['value'])
            except KeyError:
                raise KeyError("The structure of Uniprot Database changed. Please update the code")
        else:
            no_genes_name = uniprot_ac
        genes = "_".join(genes_name)
        return genes,no_genes_name
    else:
        raise requests.exceptions.HTTPError("Error during the uniprot query. Exiting...")


parser = argparse.ArgumentParser(description='get information about enzimes and catalytic sites')

parser.add_argument("-ac","--uniprot_ac",
                       dest="uniprot_ac",
                       required = True,
                       type=str,
                       help="Uniprot_ac to provide to check if the protein"\
                            " is an enzime and in case retirve the catalyitic site")

parser.add_argument("-m","--homologus_prediction",
                       dest="homologus_prediction",
                       required = False,
                       action = "store_true",
                       help="Choose to include the homologus_predictions for the "
                            " identification of catalytic site for the protein "
                            "of interest")

parser.add_argument("-o","--output",
                       dest="output_file",
                       required = False,
                       type=str,
                       help="Output file name. Default: "\
                            "uniprot_ac_catalytic_residues.csv")
args = parser.parse_args()

# URL JSON
url_all_entries = "https://www.ebi.ac.uk/thornton-srv/m-csa/api/residues/?format=json"
uniprot_url = "https://www.ebi.ac.uk/thornton-srv/m-csa/api/entries/?format=json&entries.proteins.sequences.uniprot_ids="+args.uniprot_ac
#url2= "https://www.ebi.ac.uk/thornton-srv/m-csa/api/residues/?format=json&entries.proteins.sequences.uniprot_ids=P07948" # in case we are interested only in catalytic sites
mcsa_id_url = "https://www.ebi.ac.uk/thornton-srv/m-csa/api/entries/?format=json&entries.mcsa_ids="
homologus_prediction_url = "https://www.ebi.ac.uk/thornton-srv/m-csa/api/homologues_residues.json" # for homologus catalitic residues
#file = "homologues_residues.json"

avaiable_proteins = []
entry_without_uniprot = []

# M-CSA database queries
user_request_response = requests.get(uniprot_url)
all_entries_response = requests.get(url_all_entries)

#############################################################################
#                                                                           #
#                       Homologus prediction retrival                       #
#                                                                           #
#############################################################################

if args.homologus_prediction:
    print("Accessing homologus prediction file")
    homologus_prediction_response = requests.get(homologus_prediction_url)

##########  keep only the entry associated to the target uniprot_ac #########

    if not homologus_prediction_response.status_code == 200:
        print("Error during the c-msa database query. Exiting...")
        exit(1)
    else:
        homologus_file = homologus_prediction_response.json()
        counter = 0
        mcsa_uniprot_ids = {}
        for entry in homologus_file:
            uniprot_ids = []
            for homologus in entry['residue_sequences']:
                uniprot_ids.append(homologus['uniprot_id'])
            mcsa_uniprot_ids[entry['mcsa_id']] = uniprot_ids
        filtered_mcsa = []
        for mcsa_id,uniprot_acs in mcsa_uniprot_ids.items():
            if args.uniprot_ac in uniprot_acs:
                filtered_mcsa.append(mcsa_id)

        selected_entries = []
                
        for mcsa_id in filtered_mcsa:
            for entry in homologus_file:
                if entry['mcsa_id'] == int(mcsa_id):
                    selected_entries.append(entry)

######### Collect the catalytic residues for the selected entries ###########

        homologus_df = {}
        for entry in selected_entries:
            reference_residues = ""
            target_residues = ""
            residues = entry['residue_sequences']
            mcsa_identifier = entry['mcsa_id']
            for i in residues:
                if i['is_reference'] == True:
                    reference_residues = i['code']+str(i['resid'])
                if i['uniprot_id'] ==  args.uniprot_ac:
                    target_residues = i['code']+str(i['resid'])
            information = {"reference_residue":reference_residues,
                           "target_residue":target_residues}
            if mcsa_identifier not in homologus_df.keys():
                homologus_df[mcsa_identifier] = {"residues":[information]}
            else:
                homologus_df[mcsa_identifier]["residues"].append(information)

############ Collect the information for the reference enzyme ###############

        for mcsa_id in homologus_df:
            mcsa_url = mcsa_id_url + str(mcsa_id)
            mcsa_request_response = requests.get(mcsa_url)
            if mcsa_request_response.status_code == 200:
                mcsa_entry = mcsa_request_response.json()
                for catalytic_site_information in mcsa_entry['results']:
                    ec_number = catalytic_site_information['all_ecs']
                    catalytic_residues = []
                    uniprot_reference = catalytic_site_information['reference_uniprot_id']
                    references_ac = uniprot_reference.split(", ")
                    for index, ac in enumerate(references_ac):
                        references_ac[index] = from_uniprot_ac_to_gene_name(ac)[0]
                    reference_gene_name = "_".join(references_ac)
                    reference_name = catalytic_site_information['enzyme_name']
                    homologus_df[mcsa_id].update({"ec_number":"_".join(ec_number)})
                    homologus_df[mcsa_id].update({"reference_uniprot_id":uniprot_reference})
                    homologus_df[mcsa_id].update({"reference_gene_name":reference_gene_name})

############ Convert the dictionary into a dataframe ###############

        data = []
        for mcsa_id, details in homologus_df.items():
            for residue in details['residues']:
                entry = {'Uniprot_ac': args.uniprot_ac,
                         'target_gene_name':from_uniprot_ac_to_gene_name(args.uniprot_ac)[0],
                         'target_residue': residue['target_residue'],
                         'mcsa_id_reference_entry': mcsa_id,
                         'reference_residue': residue['reference_residue'],
                         'reference_uniprot_ac': details.get('reference_uniprot_id'),
                         'reference_gene_name': details.get('reference_gene_name'),
                         'ec_number':details.get('ec_number')
                }
                data.append(entry)

####################### Write Output file ########################### 

        
        homologus_prediction_df = pandas.DataFrame(data)
        if not homologus_prediction_df.empty:
            homologus_prediction_df.to_csv(f"{args.uniprot_ac}_catalytic_sites_homologus_prediction.csv",index = False,sep=";")

# ---------- Collect the list with all the uniprot ids contained in M-CSA database ----------- #

# check successful request (status code 200)
if all_entries_response.status_code == 200:
    # Convert json content in python dictionary
    all_entries_json = all_entries_response.json()

    # extract all the available uniprot id of the M-CSA annotated entries
    for entry in all_entries_json:
            uniprot_id = entry['residue_sequences'][0]['uniprot_id']
            if uniprot_id:
                avaiable_proteins.append(uniprot_id)
            else:
                entry_without_uniprot.append(entry['mcsa_id'])

else:
    print("Error in the request to M-CSA database:", all_entries_response.status_code)


#############################################################################
#                                                                           #
#                  Mannually curated annotations retrival                   #
#                                                                           #
#############################################################################

# ----------------------- Process the user entry  ------------------------- #

entries_not_in_MCSA = []
entries_not_enzimes = []
output_data = {}
# check successful request (status code 200)
if user_request_response.status_code == 200:
    # Convert json content in python dictionary
    user_request_json = user_request_response.json()

# --- Retrive enzyme class form uniprot in case of lack of annotations --- #

    if not user_request_json['results'] and args.homologus_prediction and homologus_prediction_df.empty or not user_request_json['results'] and not args.homologus_prediction:
        entry_uniprot_information = f"https://rest.uniprot.org/uniprotkb/{args.uniprot_ac}.json"
        uniprot_response = requests.get(entry_uniprot_information)
        if uniprot_response.status_code == 200:
            uniprot_entry_information = uniprot_response.json()
            try: 
                description_keys = list(uniprot_entry_information['proteinDescription']['recommendedName'].keys())
            except KeyError:
                print("The Uniprot database changed the entry structure in json format. "
                      "The script needs to be updated")
            if 'ecNumbers' in description_keys:
                ec_numbers = []
                for i in uniprot_entry_information['proteinDescription']['recommendedName']['ecNumbers']:
                    ec_numbers.append(i['value'])
                if len(ec_numbers) == 1:
                    ec_number = ec_numbers[0]
                else:
                    ec_number = "_".join(ec_numbers)

                output_data['ec_number'] = ec_number
                output_data['catalytic residues'] = pandas.NA
                output_data["is_manually_curated"] = pandas.NA
                output_data["reference_uniprot"] = pandas.NA
                output_data['reference_enzyme_name'] = pandas.NA
                output_data['reference_gene_name'] = pandas.NA

                entries_not_in_MCSA.append(args.uniprot_ac)
            else:
                entries_not_enzimes.append(args.uniprot_ac)
        else:
            print(f"The uniprot_ac provided as input is not contained in the Uniprot database or"
                  f" it's in the wrong format. Exiting...")
            exit(1)

# --- Extract enzyme class and catalytic site annotations from the user query --- #

    else:
        for catalytic_site_information in user_request_json['results']:
            output_data['ec_number'] = catalytic_site_information['all_ecs']
            catalytic_residues = []
            uniprot_reference = catalytic_site_information['reference_uniprot_id']
            gene_name_reference = from_uniprot_ac_to_gene_name(uniprot_reference)[0]
            reference_name = catalytic_site_information['enzyme_name']
            for residue in catalytic_site_information['residues']:
                if residue['residue_sequences'][0]['uniprot_id'] == args.uniprot_ac:
                    position = str(residue['residue_sequences'][0]['resid'])
                    residue = residue['residue_sequences'][0]['code']
                    catalytic_residues.append(residue+position)
                else:
                    mcsa_id = residue['mcsa_id']
                    position = residue['residue_sequences'][0]['resid']
                    residue = residue['residue_sequences'][0]['code']
                    print(f"something strange with {residue}{position} residue. "\
                          f"It does not belong to the {args.uniprot_ac} uniprot ac ")
            output_data['catalytic residues'] = catalytic_residues
            if args.uniprot_ac in avaiable_proteins:
                isin_uniprot = "Yes"
            else:
                isin_uniprot = "No"
            output_data["is_manually_curated"] = [isin_uniprot]
            output_data["reference_uniprot"] = [uniprot_reference]
            output_data['reference_enzyme_name'] = [reference_name]
            output_data['reference_gene_name'] = [gene_name_reference]


else:
    print(f"Error in the request for {args.uniprot_ac} entry ", 
          user_request_response.status_code)

#------------------------------ Output section -----------------------------#

if entries_not_in_MCSA:
    print(f"{args.uniprot_ac} does not have any matches in M-CSA database,")
if entries_not_enzimes:
    print(f"{args.uniprot_ac} is not annotated as enzyme in Uniprot database")


if output_data:
    output_data['target_gene_name'] = from_uniprot_ac_to_gene_name(args.uniprot_ac)[0]
    columns_name = output_data.keys()
    columns_values = output_data.values()
    output = {args.uniprot_ac:columns_values}

    # DataFrame Creation 
    df = pandas.DataFrame.from_dict(output, orient = 'index', columns=columns_name)
    df = df.rename_axis('Uniprot_ac')
    df = df.applymap(lambda x: ', '.join(x) if isinstance(x, list) else x)

    columns_order = ["target_gene_name",
                     "ec_number",
                     "catalytic residues",
                     "is_manually_curated",
                     "reference_uniprot",
                     "reference_enzyme_name",
                     "reference_gene_name"]

    df = df[columns_order]
    if not args.output_file:
        if not entries_not_in_MCSA:
            out_name = args.uniprot_ac+"_catalytic_residues.csv"
        else:
            out_name = args.uniprot_ac+"_ec_number_annotation.csv"
    else:
        out_name= args.output_file
    df.to_csv(f"{out_name}",sep=";")



