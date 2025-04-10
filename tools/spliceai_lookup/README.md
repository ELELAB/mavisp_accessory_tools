# SpliceAI_lookup local API

splice_lookup.py is a python script designed to query the SpliceAI_lookup API (https://spliceailookup.broadinstitute.org) on a local server retriving the predictions about splicing alterations upon mutations from SpliceAI and Pangolin tools. 

## Requirements

The script requires the following Python version:

Python 3.10
The script also requires the following Python libraries:

pandas
request
time
argparse

In order to run the Splice_lookup locally install the requirements specified in the README.md
for the API_local_server installation

## Description

splice_lookup.py processes a CSV file containing:

- A "Mutation" column with protein sequence mutations in the format S2P.

- An "HGVSg" column with the corresponding genomic coordinates in the format hg19,12:g.133263898A>G.

The script performs the following steps:

- Check that the input file contains the required columns

- Reads the CSV file and identifies mutations with valid genomic coordinates.

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

- **acceptor gain** the nucleotide specified in the "affected_nucleotide_position" gains splicing acceptor activity

- **donor gain** the nucleotide specified in the "affected_nucleotide_position" as gains splicing donor activity

- **acceptor loss** the nucleotide specified in the "affected_nucleotide_position" looses splicing acceptor activity

- **donor loss** the nucleotide specified in the "affected_nucleotide_position" looses splicing donor activity

 For Pangolin the consequence typologies are the following:

- **splice gain** the nucleotide specified in the "affected_nucleotide_position" aquires splicing activity

- **splice loss** the nucleotide specified in the "affected_nucleotide_position" looses the splicing activity

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

- **-i, --input_flag** csv file containing the Mutations and the corresponding genomic coordinates
- **-d, --distance** defines the genomic region analyzed for splicing alterations caused by the mutation.

An optional column can be customized in order to modifiy the waiting time between two consecutive API requests:

- **-t,--time_delay** Specifies the delay between two consecutive API queries (Default: 20 seconds).

The csv file must contain the following columns:

- **Mutation** mutation in one-letter code  (i.e S2P)

- **HGVSg** genomic coordinate of the Mutation expressed in the following format ("genome_version","chromosome number":g."genomic_coordinate""WT_nucleotide">"altered_nucleotide")

Here an example of input file:

|Mutation|HGVSg|
|--------|-----|
|S2P|"hg19,12:g.133263898A>G, hg38,12:g.132687312A>G"|
|S2Y|"hg19,12:g.133263897G>T, hg38,12:g.132687311G>T"|
|L3M|"hg19,12:g.133263895G>T"|
|L3P|"hg19,12:g.133263894A>G"|

N.B the genomic coordinates must be provided as strings. If multiple coordinates (coming from two different genome annotation versions) are annotated for a single mutation, include them in the same string, separated by a comma and a space (", "). See the input file.


## Output

The script generates three types of output files:

- **spliceai_df_output.csv** a csv file containing the output file form SpliceAI tools

- **pangolin_df_output.csv** a csv file containing the output file form SpliceAI tools

- **not_found_variants.csv**

The first two files contain the follwoing columns 

- **Mutation** Input mutation

- **variant_coordinate** genomic coordinate of the variant under investigation

- **affected_nucleotide_position** position under investigation affected by the mutation

- **Δ_tpye** consequence on the position under investigation in terms of splicing acceptor or donor with gain or loss for SpliceAI, or spilicing gain or loss for Pangolin

- **Δ_score** The proability associated to that specific consequence (it's the difference between ALT_score and REF_score)

- **REF_score** The probability that a specific position in the wild-type sequence is a splicing acceptor or donor, as predicted by SpliceAI, or a splicing site, as predicted by Pangolin.

- **ALT_score**

he probability that a specific position in the mutated sequence is a splicing acceptor or donor, as predicted by SpliceAI, or a splicing site, as predicted by Pangolin.

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

N.B the output file are aggregated based on the gene_id, transcript_id and ref_seq_id

If some variants have been wrongly annotated or SpliceAI or Pangolin returned some errors during the analysis, they are annotated in the entry_not_found.csv contianing the following columns:

- **Mutation** the input mutation as reported in the input file.

- **HGVSg** the genomic coordinate as reported in the input file.

- **note** it reports the ERROR for which SpliceAI or Pangolin failed.

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
```
module load python/3.10/modulefile
python splice_lookup.py -i test.csv -d 500 -t 9
```

N.B Be sure that "Mutation" and "HGVSg" column are contained in the input file

You can run an example within the example folder
```
bash run.sh
```


