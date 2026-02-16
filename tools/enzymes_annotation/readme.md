# Enzyme annotation

The get_enzyme_catalytic_site.py is a Python script designed to query the M-CSA database to retrieve EC information and the catalytic residues associated with a given UniProt accession (uniprot_ac) provided as input.
Note: As of April 9, 2024, the M-CSA database contains 1,005 manually curated entries. However, there are homologous predictions of catalytic residues for most enzymes in the UniProt database according to these guidelines (https://www.ebi.ac.uk/thornton-srv/m-csa/documentation/#2.1.b%20-%20Navigating%20M-CSA%20using%20Statistics). Depending on the specified flags, the script reports predictions for the target protein derived from alignments and, if available, the manually curated annotations in the database.
Predictions are stored in a separate file and queried differently from the manually curated entries in the database. Depending on the specific case, the script produces output files with structures tailored to the type of information retrieved. Some entries may still have annotations in the database even if they are not part of the 1,005 manually curated entries. For these cases, the output file will follow the structure of the manually annotated entries but will include a note indicating that the information is a prediction. This is reported in the "is_manually_curated" column of the output file.
If there is neither an annotation nor a prediction for a target protein, the script will annotate only the enzyme class to which the protein belongs by querying the UniProt database directly.

## Input
uniprot_ac of protein of interest

## Usage

```
moudule load python/3.10/modulefile
python get_enzyme_catalytic_site.py -ac uniprot_ac -o output_file.csv -m
```

The flag -o is optional (Default uniprot_ac_catalytic_residues.csv)
The flag -m is optional and is used to retrieve predicted information from M-CSA database

## Output
The output files name and strucutres depend on if the target entry is annotated in the database and how it's annotated:
so if it's a prediction or a manual annotation:

In case the inuput is manually annotated in the database (or it is annotated as prediction):

**Uniprot_ac**: input uniprot_ac 
**target_gene_name**: gene name of the input protein 
**ec_number**: EC code for the enzime classification
**catalytic residues**: catalytic residues in the active site
**is_manually_curated**: Yes if it's a manually curated entry else is No
**reference_uniprot**: uniprot_ac of the manually curated entry representative of that specific EC class
**reference_enzyme_name**: enzyme name of the manually curated entry representative of that specific EC class
**reference_gene_name**: gene name of the manually curated entry representative of that specific EC class

In case the input is reported as prediction in the homologus file associated to the database:

uniprot_ac_catalytic_sites_homologus_prediction.csv

**Uniprot_ac**: input uniprot_ac 
**target_gene_name**: gene name of the input protein 
**target_residue**: potential residue in the target protein aligned with the residue of the reference enzyme
**mcsa_id_reference_entry**: id of the enzyme manually curated entry representative of that specific EC class
**reference_residue**: residue in the reference enzyme aligned with the residue in the target protein
**reference_uniprot_ac**: uniprot_ac of the manually curated entry representative of that specific EC class 
**reference_gene_name**: gene name of the manually curated entry representative of that specific EC class
**ec_number**: EC code for the enzime classification

In case the input uniprot_ac is not neither annotated in the database nor in the homologus file associated with it

uniprot_ac_ec_number_annotation.csv

The columns are the same of the previous output file but empty except for Uniprot_ac,target_gene_name,ec_number.




