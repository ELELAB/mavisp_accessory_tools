import pandas as pd
import requests
import time
import argparse

# Function to make an API call
def call_api(base_url, hg, distance, chro, coordinate, wt_nt, mut_nt):

    """
    Builds and sends a GET request to the specified API and handles the response.

    Parameters:
    -----------
    - base_url (str): The base URL of the API to query.
    - hg (str): The genome version, either '37' or '38'.
    - distance (int): The distance parameter for the SpliceAI model.
    - chro (str): The chromosome number (e.g., "1", "X").
    - coordinate (int): The genomic coordinate (e.g., "140300616").
    - wt_nt (str): The wild-type nucleotide (reference).
    - mut_nt (str): The mutant nucleotide (altered).

    Returns:
    --------
    - dict: The JSON response from the API if the request is successful (status code 200).
            Returns a dictionary with an error message if the request fails.
    """
    url = f"{base_url}/?hg={hg}&distance={distance}&variant=chr{chro}-{coordinate}-{wt_nt}-{mut_nt}"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    else:
        return {"error": f"Failed to fetch data for {variant}. Status code: {response.status_code}"}


parser=argparse.ArgumentParser(description='Plot the ΔΔG of Mutatex'\
                                           ' and Rosetta along with std')

parser.add_argument("-i","--input_file",
                       dest = "input_file",
                       required = True,
                       help = "Mavisp csv file with the Mutation and "\
                                "the genomic coordinates under the HGSv column")
parser.add_argument("-d","--distance",
                       dest = "distance",
                       required = True,
                       help = "Defines the genomic region analyzed for "\
                                 " splicing alterations caused by the mutation.")

parser.add_argument("-t","--time_delay",
                       dest = "time_delay",
                       required = False,
                       type = int,
                       default = 20,
                       help = "Specifies the delay between two consecutive API queries (Default: 20 seconds).")



args=parser.parse_args()

# Read the input CSV
df = pd.read_csv(args.input_file)

#############################################################################
#                                                                           #
#                               Error handling                              # 
#                                                                           #
#############################################################################

# Check that required columns are present
if "Mutation" not in df.columns or "HGVSg" not in df.columns:
    raise ValueError("The CSV file must contain the columns 'Mutation' and 'HGVS'.")


#############################################################################
#                                                                           #
#                          Input file Parsing                               # 
#                                                                           #
#############################################################################

# put the genomic coordinates in a list and put it in a dictionary with mutation
# as keys

mutations_genomic_coordinates = {}

for mutation,hgvsg in zip(df['Mutation'].to_list(),df['HGVSg']):
    hgvsgs = hgvsg.split(", ")
    mutations_genomic_coordinates[mutation] = hgvsgs

spliceai_output = []
pangolin_output = []
not_found = []

#############################################################################
#                                                                           #
#                               API query                                   # 
#                                                                           #
#############################################################################

for mutation, genomic_coordinates in mutations_genomic_coordinates.items():
    for genomic_coordinate in genomic_coordinates:

        # ------------------------ API querying --------------------------- #

        # Parameters for API query
        hg = genomic_coordinate.split(":")[0].split(",")[0][2:] # genome version
        if str(hg) == "19":
            hg = "37"
        chro = genomic_coordinate.split(":")[0].split(",")[1] # chromosome number
        coordinate = genomic_coordinate.split(".")[1][:-3] # mutation genomic coordinate
        wt_nt = genomic_coordinate.split(">")[0][-1] # WT nucletotide 
        mut_nt = genomic_coordinate.split(">")[1] # mutated nucleotide
        distance = args.distance # distance on the genome in which look for splicing functional element
        time_delay = args.time_delay # Seconds to wait between requests

        # first structure of the urls 
        base_urls = [f"http://127.0.0.1:8080/spliceai", 
                     f"http://127.0.0.1:8080/pangolin"] 

        for base_url in base_urls:
            result = call_api(base_url, hg, distance, chro, coordinate, wt_nt, mut_nt)

        # --------------------- API result process ------------------------ #

            # handle not found variants with that particular tool and with a 
            # specific genomic version
            if "error" in result.keys():
                row = {}
                row['Mutation'] = mutation
                row['HGVSg'] = genomic_coordinate
                row['note'] = result['error']
                not_found.append(row)

            for parameter,values in result.items():

                # SpliceAI results processs
                if "spliceai" in base_url:
                    if parameter == 'scores':
                        for transcript in values:
                            for delta_type,\
                                delta_score,\
                                position,\
                                REF_score,\
                                ALT_score in zip(['Acceptor_Loss','Donor_Loss','Acceptor_Gain','Donor_Gain'],\
                                                 ['DS_AL','DS_DL','DS_AG','DS_DG'],
                                                 ['DP_AL','DP_DL','DP_AG','DP_DG'],
                                                 ['DS_AL_REF','DS_DL_REF','DS_AG_REF','DS_DG_REF'],
                                                 ['DS_AL_ALT','DS_DL_ALT','DS_AG_ALT','DS_DG_ALT']):
                                row = {}
                                row["gene_id"] = transcript["g_id"]
                                row["transcript_id"] = transcript["t_id"]
                                row["ref_seq_id"] = transcript['t_refseq_ids']
                                row['Mutation'] = mutation
                                row["variant_coordinate"] = genomic_coordinate
                                row['affected_nucleotide_position'] = str(transcript[position])+ "bp"
                                row['Δ_tpye'] = delta_type
                                row['Δ_score'] = transcript[delta_score]
                                row['REF_score'] = transcript[REF_score]
                                row['ALT_score'] = transcript[ALT_score]
                                spliceai_output.append(row)

                # Pangolin results processs
                if "pangolin" in base_url:
                    if parameter == 'scores':
                        for transcript in values:
                            for delta_type,\
                                delta_score,\
                                position,\
                                REF_score,\
                                ALT_score in zip(['Splice_Loss','Splice_Gain'],\
                                                 ['DS_SL','DS_SG'],
                                                 ['DP_SL','DP_SG'],
                                                 ['SL_REF','SG_REF'],
                                                 ['SL_ALT','SG_ALT']):
                                row = {}
                                row["gene_id"] = transcript["g_id"]
                                row["transcript_id"] = transcript["t_id"]
                                row['ref_seq_id'] = transcript['t_refseq_ids']
                                row['Mutation'] = mutation
                                row["variant_coordinate"] = genomic_coordinate
                                row['affected_nucleotide_position'] = str(transcript[position])+ "bp"
                                row['Δ_tpye'] = delta_type
                                row['Δ_score'] = transcript[delta_score]
                                row['REF_score'] = transcript[REF_score]
                                row['ALT_score'] = transcript[ALT_score]
                                
                                pangolin_output.append(row)

        # Wait before making the next request
        time.sleep(time_delay)  

#############################################################################
#                                                                           #
#                            Output_production                              #  
#                                                                           #
#############################################################################

spliceai_df = pd.DataFrame(spliceai_output)
pangolin_df = pd.DataFrame(pangolin_output)

# ----------------------- Process the output df --------------------------- #

for tool,df in zip(['spliceai',
                    'pangolin'],[spliceai_df,pangolin_df]):

    df = df.sort_values(by=['Mutation', 'gene_id'], ascending=[False, True])
    df.to_csv(f"{tool}_output.csv",index=False)


# --------------- Process the not found entry output df ------------------ #

not_found_df = pd.DataFrame(not_found)

if not not_found_df.empty:
    not_found_df.to_csv("not_found_variants.csv",index=False)
else:
    print("All the variants have been analyzed.")


