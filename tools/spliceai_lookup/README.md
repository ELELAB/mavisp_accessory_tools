# SpliceAI_lookup local API

splice_lookup.py is a python script designed to query the SpliceAI_lookup API (https://spliceailookup.broadinstitute.org) on a local server retriving the predictions about splicing alterations upon mutations (SNV and SV) from SpliceAI and Pangolin tools. 

## Requirements

The script requires the following Python version:

Python 3.10
The script also requires the following Python libraries:

pandas
request
time
argparse
yaml
os
pyfaidx
sys
mygene
re

In order to run the Splice_lookup locally install the requirements specified in the README.md
for the API_local_server installation

## Description

splice_lookup.py supports two types of files:

One file (MAVISp format) requires:

- A "Mutation" column with protein sequence mutations in the format S2P.

- An "HGVSg" column with the corresponding genomic coordinates in the format hg19,12:g.133263898A>G, hg19,5:g.79970811_79970812delinsTT ecc.

One file (MAF format) requires:

- A "NCBI_Build" column with the genome build information (GRCh38,GRCh19 or GRCh37).

- A "Chromosome" column with the chromosome number in which the variant occurs.

- A "Start_Position" column with genomic coordinates indicating the starting position of the mutation

- A "End_Position" column with genomic coordinates indicating the ending position of the mutation

- A "Reference_Allele" column with the nucleotide presents in the wild-type (reference) sequence at the mutated position.

- A "Tumor_Seq_Allele2" column with nucleotide observed in the altered (mutant) sequence at the same position.

- A "Transcript_ID" column with the information about the trancript id to which the variant under investigation belongs.

- A "HGVSp_Short" column with the mutation expressed in short HGVSp format (i.e p.R2236C. In case of deletion or insertion with an unclear effect also "?" is accepted) 


N.B. When using a MAF file, a configuration file specifying the required column names must be provided. This is necessary because MAF files report the mutational state of a given genomic position in two separate columns, each corresponding to one of the two alleles. Therefore, the user must indicate which column contains the mutation of interest to be analyzed.


The script performs the following steps:

- Check depending on the type of file provided that the input file contains the required columns

- For every line in the csv, it creates a mutation object in the HGVSg format ready to be converted in the input format required by the API

- Identifies mutations with valid genomic coordinates.

- For each valid mutation, it queries the local API, submitting the variant to SpliceAI and Pangolin for splicing alteration analysis.

- Collects the results, including splicing alteration scores, and saves them to a new CSV file.

- Generates an additional CSV file listing mutations that resulted in an error message from the API, where SpliceAI or Pangolin failed to process the input.

The following details are retrieved during the API queries to SpliceAI and Pangolin:

### genomic_coordinates

The genomic coordinate of the mutation provided as input to Pangolin or SpliceAI.

### affected_nucleotide_position

The position within the specified window (default: ±500 bp) that may be affected by the variant under investigation. This includes potential splice site alterations such as:

From **SpliceAI**:

- acceptor gain

- donor gain

- acceptor loss

- donor loss

From **Pangolin**:

- splice gain

- splice loss 

N.B The position is calculated relative to the genomic coordinate of the variant. For example, "+50" indicates that the affected site is 50 nucleotides downstream of the variant.

### Δ_tpye

The type of splicing alteration caused by the variant on functional splice sites.

- **acceptor gain:** the nucleotide specified in the "affected_nucleotide_position" gains splicing acceptor activity

- **donor gain:** the nucleotide specified in the "affected_nucleotide_position" as gains splicing donor activity

- **acceptor loss:** the nucleotide specified in the "affected_nucleotide_position" looses splicing acceptor activity

- **donor loss:** the nucleotide specified in the "affected_nucleotide_position" looses splicing donor activity

 For Pangolin the consequence typologies are the following:

- **splice gain:** the nucleotide specified in the "affected_nucleotide_position" aquires splicing activity

- **splice loss:** the nucleotide specified in the "affected_nucleotide_position" looses the splicing activity

### Δ_score,

The score represents the difference between the REF and ALT scores obtained from SpliceAI and Pangolin for a specific position. These scores reflect the likelihood of the position (specified in the "affected_nucleotide_position" value) becoming an acceptor or a donor site, or losing its acceptor or donor activity upon mutation (SpliceAI), or becoming a splicing site or loosing its splicing site properties(Pangolin).

### REF_score

SpliceAI's computed probability that the given position is a splice acceptor or donor based on the reference haplotype sequence. 
Pangolin's computed probability that the given position is a splice site based on the reference haplotype sequence.

### ALT_score

SpliceAI's computed probability that the given position is a splice acceptor or donor based on the alternate haplotype sequence containing the mutation (as specified in the "genomic_coordinates" value). 
Pangolin's computed probability that the given position is a splice site based on the alternate haplotype sequence containing the mutation (as specified in the "genomic_coordinates" value).

### gene_id

The ENSG (Ensembl) code representing the gene that contains the mutation provided in the input.

### transcript_id,

The ENST (Ensembl Transcript) code representing the transcript that contains the mutation provided in the input.

### ref_seq_id

The RefSeq identifier of the protein (if applicable) affected by the mutation, assuming the transcript encodes a protein.

The script aggregates this information into a CSV file, grouping the entries by gene_id, transcript_id, and ref_seq_id. The resulting CSV also includes the original mutations on the protein sequence provided in the input.


## Input

The script requires the following flags:

- **-i, --input_flag:** csv file containing the Mutations and the corresponding genomic coordinates
- **-d, --distance:** defines the genomic region analyzed for splicing alterations caused by the mutation.
- **-g, --genome_build_path:** defines the path containing the sequenced genome in .fa format for the 38 and 19 (37) genome build
- **-m, --mavisp** flag to be specified in case the input file follow the mavsip format
- **-f, --maf:** flag to be specified in case the input file follow the MAF format.
- **-c, --config_file:** config file specifying the name of the columns required by the scripts (only when -f,--maf flag is specified)

-m and -f columns are mutually exclusive.

An optional flag can be specified in order to modifiy the waiting time between two consecutive API requests:

- **-t,--time_delay:** Specifies the delay between two consecutive API queries (Default: 20 seconds).

### Mavisp format

The csv file in the MAVISp format must be provided along with the -m --mavisp flag and it must contain the following columns:

- **Mutation:** mutation in one-letter code  (i.e S2P)

- **HGVSg:** genomic coordinate of the Mutation expressed in the following format ("genome_version","chromosome number":g."genomic_coordinate""WT_nucleotide">"altered_nucleotide")

Here an example of input file in the MAVISp format:

|Mutation|HGVSg|
|--------|-----|
|S2P|"hg19,12:g.133263898A>G, hg38,12:g.132687312A>G"|
|S2Y|"hg19,12:g.133263897G>T, hg38,12:g.132687311G>T"|
|L3M|"hg19,12:g.133263895G>T"|
|L3P|"hg19,12:g.133263894A>G"|
|P346L|"hg19,5:g.79970811_79970812delinsTT"|

N.B the genomic coordinates must be provided as strings. If multiple coordinates (coming from two different genome annotation versions) are annotated for a single mutation, include them in the same string, separated by a comma and a space (", "). See the input file.

### MAF format

The csv file in MAF format be provided along with the -f --maf flag nad it must contain the following columns:

- **NCBI_Build:** information about the genome reference. it must be GRCh38 or GRCh37.

- **Chromosome:** a str or int object indicating a number between 1 to 22 or X or Y for not autosomal chromosomes

- **Start_Position:**  a str or int object representing a number indicating the starting position of the mutation on the genome

- **End_Position:** a str or int object representing a number indicating the ending position of the mutation on the genome

- **Reference_Allele:** str indicating the WT state of the nucleotide/s affected by the mutation

- **Tumor_Seq_Allele2:** str indicating the mutated nucleotide/s affected 

- **Transcript_ID:** str indicating the trancript id to which the variant under investigation belongs.

- **HGSVp_Short:** str indicating the mutation expressed in short HGVSp format (i.e p.T462M)

Here an example of input file in MAF format:
|Transcript_ID|NCBI_Build|Chromosome|Start_Position|End_Position|Reference_Allele|Tumor_Seq_Allele2|HGVSp_Short|
|-------------|----------|----------|--------------|------------|----------------|-----------------|-----------|
|ENST00000265316|GRCh38|2|219214390|219214390|G|A|p.T462M|
|ENST00000331849|GRCh38|16|20437278|20437278|G|A|p.G483R|
|ENST00000361673|GRCh37|18|58517065|58517065|C|T|p.G1928E|
|ENST00000379370|GRCh38|1|1042601|1042601|A|AGAGAG|-|p.?ins|
|ENST00000379370|GRCh38|1|1042466|1042469|GGGC|G|p.?del|

N.B  In MAF files, the status of both alleles for each mutation is reported. Remember to specify in the configuration file the correct column name corresponding to the allele carrying the mutation. In the example, the mutation was present only in Allele2 (column Tumor_Seq_Allele2).

The name of the columns can be different and needs to be specified through a config file

### Config file

Here an example of config file

```
maf_file_columns:
  NCBI_Build: "NCBI_Build"
  Chromosome: "Chromosome"
  Start_Position: "Start_Position"
  End_Position: "End_Position"
  wt: "Reference_Allele"
  mutant: "Tumor_Seq_Allele2"
  Transcript_ID: "Transcript_ID"
  HGSVp_Short: "HGVSp_Short"
```
Modify only the strings after the ":"


## Output

The script generates three types of output files:

- **spliceai_df_output.csv:** a csv file containing the output file form SpliceAI tools

- **pangolin_df_output.csv:** a csv file containing the output file form SpliceAI tools

- **not_found_variants.csv**

### spliceai_df_output.csv and angolin_df_output.csv

The first two files contain the following columns 

- **Mutation:** Input mutation

- **variant_coordinate:** genomic coordinate of the variant under investigation

- **affected_nucleotide_position:** position under investigation affected by the mutation

- **Δ_tpye** consequence on the position under investigation in terms of splicing acceptor or donor with gain or loss for SpliceAI, or spilicing gain or loss for Pangolin

- **Δ_score** The proability associated to that specific consequence (it's the difference between ALT_score and REF_score)

- **REF_score** The probability that a specific position in the wild-type sequence is a splicing acceptor or donor, as predicted by SpliceAI, or a splicing site, as predicted by Pangolin.

- **ALT_score** The probability that a specific position in the mutated sequence is a splicing acceptor or donor, as predicted by SpliceAI, or a splicing site, as predicted by Pangolin.

- **gene_id**
The Ensembl id of the gene in which the mutation has been found

- **transcript_id**

The Ensembl id of the transcrit in which the mutation has been found

- **ref_seq_id**

the ref_seq of the protein in which the mutation has been found

Here an example of output file from **SpliceAI**:

|Mutation|variant_coordinate|affected_nucleotide_position|Δ_tpye|Δ_score|REF_score|ALT_score|gene_id|transcript_id|ref_seq_id|
|--------|------------------|----------------------------|------|-------|---------|---------|-------|-------------|----------|
|S2Y|"hg38,12:g.132687311G>T"|-137bp|Acceptor_Gain|0.00|0.00|0.00|ENSG00000177084.19|"ENST00000320574.10 ENST00000535270.5"|"None ['NM_006231.4']"|
|S2Y|"hg38,12:g.132687311G>T"|-190bp|Acceptor_Loss|0.00|0.00|0.00|ENSG00000177084.19|"ENST00000320574.10, ENST00000535270.5"|"None, ['NM_006231.4']"|
|S2Y|"hg38,12:g.132687311G>T"|-190bp|Donor_Loss|0.01|0.52|0.52|ENSG00000177084.19|ENST00000320574.10|['NM_006231.4']
|S2Y|"hg38,12:g.132687311G>T"|-190bp|Donor_Loss|0.01|0.59|0.58|ENSG00000177084.19|ENST00000535270.5|None|


And from **Pangolin**:

|Mutation|variant_coordinate|affected_nucleotide_position|Δ_tpye|Δ_score|REF_score|ALT_score|gene_id|transcript_id|ref_seq_id|
|--------|------------------|----------------------------|------|-------|---------|---------|-------|---------------|----------|
|S2Y|"hg38,12:g.132687311G>T"|-190bp|Splice_Loss|-0.00|0.13|0.12|ENSG00000177084.19|ENST00000320574.10|['NM_006231.4']|
|S2Y|"hg38,12:g.132687311G>T"|-57bp|Splice_Gain|0.02,0.81,0.82|ENSG00000177084.19|ENST00000320574.10|['NM_006231.4']|
|S2Y|"hg19,12:g.133263897G>T"|-190bp|Splice_Loss|-0.00,0.13,0.12|ENSG00000177084.19_18|ENST00000320574.10_8|['NM_006231.4']|

N.B the output file are aggregated based on the gene_id, transcript_id and ref_seq_id. In case of a input file in the MAF format the output file will be filtered keeping only those variants belonging to the transcript ID contained in the input file.

The API and webserver may not be always synchronized. Therefore, predictions for some transcript IDs might occasionally be missing when querying the API. However, predictions for the MANE transcript are always guaranteed to be available.

### not_found_variants.csv

If some variants have been wrongly annotated or SpliceAI or Pangolin returned some errors during the analysis, they are annotated in the entry_not_found.csv containing the following columns:

- **Mutation:** the input mutation as reported in the input file.

- **HGVSg:** the genomic coordinate as reported in the input file.

- **note:** it reports the ERROR for which SpliceAI or Pangolin failed.

Here an example of entry_not_found.csv file:

|Mutation|HGVSg|note|
|--------|-----|----|
|S2Y|"hg19,12:g.133263890G>T"|ERROR: The SpliceAI model did not return any scores for chr12-133263890-G-T. This may be due to the variant falling outside of all Gencode exons and introns.|
|S2Y|"hg19,12:g.133263890G>T"|ERROR: Pangolin was unable to compute scores for this variant|


## Usage 

### API activation

Activate the local API.

In order to open the local API open a new terminal and run the following commands:
```
source /data/user/marnaudi/spliceai_lookup/spliceai-env/bin/activate
conda activate /home/marnaudi/.conda/envs/bioenv/
```
run the following command to check if there are some running process:

```
lsof -i :8080 
``` 
you need to kill all the processes with your name in the port 8080

To kill them, check the PID code and run the following command for each process to kill
```
kill -9 PID
```
Check that no process is active with
```
lsof -i :8080
```
Activate the local API server with
```
cd /data/user/marnaudi/spliceai_lookup/SpliceAI-lookup/
bash start_local_server.sh 

```
### splice_lookup.py 
In a new terminal in the folder with the input files, config file and the script:

```
source /data/user/marnaudi/spliceai_lookup/spliceai-env/bin/activate
conda activate /home/marnaudi/.conda/envs/bioenv/
python splice_lookup.py -i mavisp_input.csv -g /data/databases/genome_annotation/ -d 500 -t 9 -m 
python splice_lookup.py -i maf_input.csv -g /data/databases/genome_annotation/ -d 500 -t 9 -f -c config_maf.yaml
```

N.B Be sure that the required columns are contained in the input file

You can try an example located in the example folder, which contains sample runs for both MAF and MAVISP formats.
Inside maf_file or mavisp file folder:
```
bash run.sh
```


