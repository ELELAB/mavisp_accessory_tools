import pandas as pd
import requests
import time
import argparse
import yaml
import os
from pyfaidx import Fasta
import sys
# Aggiungi la cartella dove pip ha installato mygene
sys.path.append("/home/marnaudi/.local/lib/python3.8/site-packages")
import mygene
import re

mg = mygene.MyGeneInfo()

def construct_ref_alt(row, fasta, genome_build,maf_cols):
    """
    Build HGVSg string compatible with your API for SNV, insertion, deletion, and delins.
    """
    chrom = str(row[maf_cols["Chromosome"]])
    start = int(row[maf_cols["Start_Position"]])
    end = int(row[maf_cols["End_Position"]])
    ref = str(row[maf_cols["wt"]])
    alt = str(row[maf_cols["mutant"]])
    # --- Inversion check ---
    if alt.upper() == "INV" or alt.upper() == "inversion":
        ref_seq = fasta[chrom][start-1:end].seq.upper()
        # reverse complement
        complement = str.maketrans("ACGT", "TGCA")
        alt_seq = ref_seq[::-1].translate(complement)
    # Handle different variant types
    elif ref == "-" and alt != "-":  # insertion
        # anchor = previous base
        anchor = fasta[chrom][start-1].seq.upper()
        ref_seq = anchor
        alt_seq = anchor + alt  # full inserted sequence
    elif alt == "-" and ref != "-":  # deletion
        ref_seq = fasta[chrom][start-1:end].seq.upper()  # full deleted sequence
        alt_seq = ref_seq[0]  # keep first base as anchor for HGVS compatibility
    elif len(ref) > 1 and len(alt) > 1 and ref != alt:  # delins
        ref_seq = fasta[chrom][start-1:end].seq.upper()
        alt_seq = alt  # full inserted sequence
    else:  # SNV
        ref_seq = ref
        alt_seq = alt

    # Build HGVSg string: genome_build, chromosome:g.positionREF>ALT
    hgvsg = f"{genome_build},{chrom}:g.{start}{ref_seq}>{alt_seq}"
    return hgvsg


# Function to make an API call
def call_api(base_url, hg, distance, chro, coordinate, wt_nt, mut_nt):
    """
    Builds and sends a GET request to the specified API and handles the response.
    """
    url = f"{base_url}/?hg={hg}&distance={distance}&variant=chr{chro}-{coordinate}-{wt_nt}-{mut_nt}"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()
    else:
        return {"error": f"Failed to fetch data for chr{chro}-{coordinate}. Status code: {response.status_code}"}

def parse_hgvs_for_api(genomic_coordinate, fasta_19, fasta_38):
    """
    Convert HGVSg (e.g. hg19,2:g.48018217_48018218delinsTT)
    into 4 fields for API: chr, coordinate, ref, alt
    Produces variants like:
      chr8-140300616-T-G
      chr1-1042601-A-AGAGAG
      chr1-1042466-GGGC-G
    """
    build = genomic_coordinate.split(",")[0]
    hg = build.replace("hg", "")
    chro = genomic_coordinate.split(",")[1].split(":")[0]
    g_notation = genomic_coordinate.split(":")[1][2:]  # remove 'g.'

    # Choose FASTA
    fasta = fasta_19 if hg == "19" else fasta_38

    # Normalize chromosome name (adds "chr" prefix if missing)
    fasta_chro = chro if chro in fasta.keys() else f"chr{chro}"

    # --- SNV ---
    if ">" in g_notation:
        m = re.match(r"(\d+)([ACGT]+)>([ACGT]+)", g_notation)
        if not m:
            raise ValueError(f"Unrecognized SNV format: {g_notation}")
        coord, ref, alt = m.groups()
        return chro, coord, ref, alt

    # --- DELINS ---
    elif "delins" in g_notation:
        m = re.match(r"(\d+)_(\d+)delins([ACGT]+)", g_notation)
        if not m:
            raise ValueError(f"Unrecognized delins format: {g_notation}")
        start, end, alt = m.groups()
        ref = fasta[fasta_chro][int(start)-1:int(end)].seq.upper()
        return chro, start, ref, alt

    # --- INSERTION ---
    elif "ins" in g_notation:
        m = re.match(r"(\d+)_(\d+)ins([ACGT]+)", g_notation)
        if not m:
            raise ValueError(f"Unrecognized insertion format: {g_notation}")
        start, end, alt = m.groups()
        anchor = fasta[fasta_chro][int(start)-1].seq.upper()
        ref = anchor
        alt = anchor + alt
        return chro, start, ref, alt

    # --- DELETION ---
    elif "del" in g_notation:
        m = re.match(r"(\d+)(?:_(\d+))?del", g_notation)
        if not m:
            raise ValueError(f"Unrecognized deletion format: {g_notation}")
        start, end = m.groups()
        end = end or start
        ref = fasta[fasta_chro][int(start)-1:int(end)].seq.upper()
        anchor = fasta[fasta_chro][int(start)-2].seq.upper()
        ref = anchor + ref
        alt = anchor
        return chro, str(int(start)-1), ref, alt
    
    # --- INVERSION ---
    elif "inv" in g_notation:
        m = re.match(r"(\d+)_(\d+)inv", g_notation)
        if not m:
            raise ValueError(f"Unrecognized inversion format: {g_notation}")
        start, end = m.groups()
        ref_seq = fasta[fasta_chro][int(start)-1:int(end)].seq.upper()
        alt_seq = ref_seq[::-1]  # inverti la sequenza come rappresentazione semplice
        return chro, start, ref_seq, alt_seq
        
    else:
        raise ValueError(f"Unrecognized HGVSg format: {g_notation}")



parser = argparse.ArgumentParser(description="Run SpliceAI and Pangolin API queries and annotate results with gene symbols.")
parser.add_argument("-i", "--input_file", required=True, help="CSV with 'Mutation' and 'HGVSg' columns")
parser.add_argument("-c","--config_file", type=str, help="YAML config specifying column names (required if using -f/--maf)")
parser.add_argument("-d", "--distance", required=True, help="Distance parameter for splicing models")
parser.add_argument("-t", "--time_delay", type=int, default=20, help="Delay between API queries (default: 20s)")
parser.add_argument("-g", "--genome_build_path", type=str, help="Path containing the sequenced genome in .fa format for the 38 and 19 (37) genome build")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-m","--mavisp", action="store_true", help="Use this if input has Mutation | HGVSg format")
group.add_argument("-f","--maf", action="store_true", help="Use this if input is a full MAF file | HGVSc format")
args = parser.parse_args()
if args.maf and not args.config_file:
    raise ValueError("A config file (-c/--config) must be provided when using -f/--maf")

# Read the input CSV
df = pd.read_csv(args.input_file)
mutations_genomic_coordinates = {}

print(f"[+] Loading genome reference from hg38 ...")
fasta_38 = Fasta(f"{args.genome_build_path}/hg38.fa")
print(f"[+] Loading genome reference from hg19 ...")
fasta_19 = Fasta(f"{args.genome_build_path}/hg19.fa")
if args.mavisp:
    if "Mutation" not in df.columns or "HGVSg" not in df.columns:
        raise ValueError("The CSV file must contain the columns 'Mutation' and 'HGVSg'.")
    # Map mutations to genomic coordinates
    mutations_genomic_coordinates = {
            mutation: hgvsg.split(", ") for mutation, hgvsg in zip(df["Mutation"], df["HGVSg"])
            }
if args.maf:
    with open(args.config_file, "r") as f:
        cfg = yaml.safe_load(f)

    maf_cols = cfg.get("maf_file_columns", {})

    required_cols = list(maf_cols.values())
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Column '{col}' missing in {args.input_file} MAF file")
    for _, row in df.iterrows():
        genome_build = row[maf_cols["NCBI_Build"]]
        num = re.search(r'\d+', genome_build)
        if num:
            genome_build_digit = str(num.group())
        else:
            print("Error: The 'NCBI_Build' column contains genome build annotations "\
                  "in an invalid format. Accepted values are 'GRCh37' or 'GRCh38'. Exiting...")
            sys.exit(1)
        if genome_build == "GRCh38":
            fasta_path = fasta_38
        if genome_build == "GRCh37":
            fasta_path = fasta_19
        hgvsg = construct_ref_alt(row,fasta_path,genome_build_digit,maf_cols)
        mutation = row[maf_cols["HGVSp_Short"]]
        # Map mutations to genomic coordinates
        if mutation in mutations_genomic_coordinates:
            mutations_genomic_coordinates[mutation].append(hgvsg)
        else:
            mutations_genomic_coordinates[mutation] = [hgvsg]   

spliceai_output = []
pangolin_output = []
not_found = []
# ---------------------- API Query Loop ---------------------- #
for mutation, genomic_coordinates in mutations_genomic_coordinates.items():
    for genomic_coordinate in genomic_coordinates:
        hg = genomic_coordinate.split(":")[0].split(",")[0][2:]
        try:
            chro, coordinate, wt_nt, mut_nt = parse_hgvs_for_api(genomic_coordinate, fasta_19, fasta_38)
            hg = "37" if "hg19" in genomic_coordinate else "38"
        except Exception as e:
            print(f"Error parsing {genomic_coordinate}: {e}")
            not_found.append({
                "Mutation": mutation,
                "HGVSg": genomic_coordinate,
                "note": f"Parsing error: {e}"
            })
            continue
        distance = args.distance
        time_delay = args.time_delay

        base_urls = [
            "http://127.0.0.1:8080/spliceai",
            "http://127.0.0.1:8080/pangolin"
        ]

        for base_url in base_urls:
            result = call_api(base_url, hg, distance, chro, coordinate, wt_nt, mut_nt)

            if "error" in result.keys():
                not_found.append({
                    "Mutation": mutation,
                    "HGVSg": genomic_coordinate,
                    "note": result["error"]
                })
                continue

            for parameter, values in result.items():
                # SpliceAI results
                if "spliceai" in base_url and parameter == "scores":
                    for transcript in values:
                        for delta_type, delta_score, position, REF_score, ALT_score in zip(
                            ["Acceptor_Loss", "Donor_Loss", "Acceptor_Gain", "Donor_Gain"],
                            ["DS_AL", "DS_DL", "DS_AG", "DS_DG"],
                            ["DP_AL", "DP_DL", "DP_AG", "DP_DG"],
                            ["DS_AL_REF", "DS_DL_REF", "DS_AG_REF", "DS_DG_REF"],
                            ["DS_AL_ALT", "DS_DL_ALT", "DS_AG_ALT", "DS_DG_ALT"]
                        ):
                            row = {
                                "gene_id": transcript["g_id"],
                                "transcript_id": transcript["t_id"],
                                "ref_seq_id": transcript["t_refseq_ids"],
                                "Mutation": mutation,
                                "variant_coordinate": genomic_coordinate,
                                "affected_nucleotide_position": str(transcript[position]) + "bp",
                                "Δ_type": delta_type,
                                "Δ_score": transcript[delta_score],
                                "REF_score": transcript[REF_score],
                                "ALT_score": transcript[ALT_score]
                            }
                            spliceai_output.append(row)

                # Pangolin results
                if "pangolin" in base_url and parameter == "scores":
                    for transcript in values:
                        for delta_type, delta_score, position, REF_score, ALT_score in zip(
                            ["Splice_Loss", "Splice_Gain"],
                            ["DS_SL", "DS_SG"],
                            ["DP_SL", "DP_SG"],
                            ["SL_REF", "SG_REF"],
                            ["SL_ALT", "SG_ALT"]
                        ):
                            row = {
                                "gene_id": transcript["g_id"],
                                "transcript_id": transcript["t_id"],
                                "ref_seq_id": transcript["t_refseq_ids"],
                                "Mutation": mutation,
                                "variant_coordinate": genomic_coordinate,
                                "affected_nucleotide_position": str(transcript[position]) + "bp",
                                "Δ_type": delta_type,
                                "Δ_score": transcript[delta_score],
                                "REF_score": transcript[REF_score],
                                "ALT_score": transcript[ALT_score]
                            }
                            pangolin_output.append(row)

        time.sleep(time_delay)

# ---------------------- Create DataFrames ---------------------- #
spliceai_df = pd.DataFrame(spliceai_output)
pangolin_df = pd.DataFrame(pangolin_output)
not_found_df = pd.DataFrame(not_found)

# ---------------------- Annotate with Gene Symbols ---------------------- #
for name, df_tool in [("spliceai", spliceai_df), ("pangolin", pangolin_df)]:
    if not df_tool.empty:
        print(f" Annotating {name} results with gene symbols...")
        gene_ids = [g.split(".")[0] for g in df_tool["gene_id"].dropna().unique()]
        mapping = mg.querymany(gene_ids, scopes="ensembl.gene", fields="symbol", species="human", as_dataframe=False)
        gene_map = {r["query"]: r.get("symbol", "NA") for r in mapping}
        df_tool["gene_symbol"] = df_tool["gene_id"].apply(lambda g: gene_map.get(g.split(".")[0], "NA"))
        df_tool = df_tool.sort_values(by=["Mutation", "gene_id"], ascending=[False, True])
        if args.maf:
            patterns = df['Transcript_ID'].dropna().astype(str)
            regex = r'^(' + '|'.join(map(re.escape, patterns)) + r')(\.|$)'

            df_tool = df_tool[
                df_tool['transcript_id']
                    .astype(str)
                    .str.match(regex, na=False)
            ]
        df_tool.to_csv(f"{name}_output.csv", index=False)
        print(f" {name}_output.csv saved with gene symbols.")
    else:
        print(f"No results found for {name}.")

if not not_found_df.empty:
    not_found_df.to_csv("not_found_variants.csv", index=False)
    print(" Some variants were not found. Saved to not_found_variants.csv.")
else:
    print(" All variants analyzed successfully.")

