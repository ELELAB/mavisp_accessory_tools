# Accessory Tools for MAVISp-Related Research  

This repository provides accessory tools associated with papers and activities of our research group that are not part of the core **MAVISp** framework.  

## Current Tools  

- **pfam_scrape**: A tool for defining non-redundant datasets, developed for a study currently under review at *Journal of Molecular Biology (JMB)*:  
  *“Deciphering Long-Range Effects of Mutations: An Integrated Approach Using Elastic Network Models and Protein Structure Networks.”*  
- **splice_lookup.py**: is a python script designed to query the SpliceAI_lookup API (https://spliceailookup.broadinstitute.org) on a local server retriving the predictions about splicing alterations upon mutations from SpliceAI and Pangolin tools. developed for a study currently under review at "XXX":
  *"From Structure to Function: Interpreting the Effects of DNA Polymerase Variants in Cancer"*
- **get_second_sphere_residues.py** Python script designed to identify the second coordination sphere residues of catalytic residues provided as input.
  *"From Structure to Function: Interpreting the Effects of DNA Polymerase Variants in Cancer"*
- **ArpeggioInteractionPipeline.py**  Python script for analyzing atomic interactions in protein and molecular structures using Arpeggio. It supports conversion of PDB files to CIF format, validation of residues from input files, computation of intra- and inter-molecular contacts, filtering of low-quality or clashing contacts, and generation of single-point or saturation mutation lists for residues involved in interactions. 
- **alphafill_cofactor_survey** Python script querying AlphaFill for a Uniprot AC or a list of Uniprot AC's to find and list potential cofactors.
- **stability_benchmark_sim_length** uses a Python script that visualizes the performance of different simulation times against a gold standard using confusion matrices, calculating key metrics such as sensitivity, specificity, accuracy, precision, F1 score, and structure length for one or multiple genes.
- **consensus_benchmark_sim_length** is a Python script that generates confusion matrix plots and calculates performance metrics for MAVISP stability classification across ensemble datasets from different simulation times, comparing them to a gold-standard time and producing both gene-specific and aggregated results with supporting CSV files. 
- **globcheck**: A tool for collecting 3 metrics of globularity (ashpericity, normalised radius of gyration, packing density) of a protein structure. Developed for a study currently under review at *Journal of Molecular Biology (JMB)*:  
  *“Deciphering Long-Range Effects of Mutations: An Integrated Approach Using Elastic Network Models and Protein Structure Networks.”*  
- **coco_sum**: A tool for collecting and parsing coiled-coil predictions using external predictor CoCoNat (https://github.com/BolognaBiocomp/coconat). Developed for a study currently under review at *Journal of Molecular Biology (JMB)*:  
  *“Deciphering Long-Range Effects of Mutations: An Integrated Approach Using Elastic Network Models and Protein Structure Networks.”*  
- **af\_tools**: A directory containing utilities to generate or process AlphaFold model structures for downstream local interaction analysis. More specifically:
	1. mavisp\_af3: a toolkit for processing the .zip package obtained from the AlphaFold3 server. It converts the predicted structures from .cif to .pdb, computes pdockq2 scores for multi-chain models, and outputs a summary table collecting all available confidence metrics for each generated structure.
	2. mavisp\_afmultimer: a toolkit that generates protein complex structures with AlphaFold-Multimer from provided FASTA sequences, computes pdockq2 scores, and outputs a summary table of confidence metrics in the same format as mavisp\_af3.
- **uniprot2refseq**: a script to retrieve the correct refseq ID given a uniprot ID - to use before the mavisp_automatization step to ensure the proper match between identifiers
  
## Citations  

If you use this repository, please cite the following works:  

- Arnaudi, Matteo, et al. *"MAVISp: A Modular Structure-Based Framework for Protein Variant Effects."* bioRxiv (2024).  
- Krzesińska, Karolina, et al. *"Deciphering Long-Range Effects of Mutations: An Integrated Approach Using Elastic Network Models and Protein Structure Networks."* Manuscript under review.
- Arnaudi, Matteo et al. *"From Structure to Function: Interpreting the Effects of DNA Polymerase Variants in Cancer"*.Manuscript under review  
- Kishore, Jaganathan et al. *"Predicting Splicing from Primary Sequence with Deep Learning"*
- Zeng, Tony et al. *"Predicting RNA splicing from DNA sequence using Pangolin"*
