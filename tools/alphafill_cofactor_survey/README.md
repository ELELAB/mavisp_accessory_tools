# **Alphafill cofactor survey**

## **Description**  
The `find_cofactor.py` script queries ALphafill for a Uniprot AC or a list of Uniprot ACs. It then extracts the list of potential cofactors for this/these protein(s). 
Only Alphafill results with identity >=30% are included

To run the script, the user must provide a single Uniprot AC or a MAVISp database index file (csv) including the column "UniProt AC".
The user can adjust the name of the output folder with an optional last argument

### **Output**
| Protein          | Heteroatom                |
|------------------|---------------------------|
| `Uniprot AC`     | `heteroatom1 heteroatom2` |
| `P00325` 	   | `ACN APR CHD CNH CO`      |

The output contains only heteroatoms found for PDB entries with >=30% identity. 


If the flag -f is used and you provide a txt file for filtering, the output will only contain heteroatoms from the list. All other heteroatoms will be ignored. 
The filtered output will have the same format:

| Protein          | Cofactor                |
|------------------|-------------------------|
| `Uniprot AC`     | `cofactor1 cofactor2  ` |
| `P00325`         | `ACN APR CHD CNH CO`    |



Finally, a txt file is provided with a list of unique heteroatoms found for all proteins in the input file. 


## **Requirements**
module load python/3.10/modulefile 

## **Usage**
The script uses 4 flags. Mandatory to use either -u or -i
-u: input uniprot AC
-i: input file (csv) containing a column "UniProt AC"
-o: output file prefix. If not provided, default is summary_output
-f: Filtering file (txt). To use if only specific cofactors/heteroatoms should be included in the output. 


# Run for one Uniprot AC
python find_cofactor.py -u P00533

# Run for file with more than one Uniprot AC
python find_cofactor.py -i index.csv

# Run with cofactor filter list (to only keep heteroatoms that are in that list)
python find_cofactor.py -i index.csv -f cofactors.txt -o summary_output



## **Example**
cd example/find_cofactor
python ../../find_cofactor.py -i index.csv -f ../annotate_heteroatoms/all_cofactors.txt -o summary_output

# **Annotate cofactors**
The annotate_cofactors.py script uses a json file to annotate cofactors from the PDB. The json file is a dictionary of cofactors in the PDB, and is compared to a user defined list of heteroatoms. 
The json dictionary can be downloaded from here: https://www.ebi.ac.uk/pdbe/api/pdb/compound/cofactors
This resource was described in: Mukhopadhyay, A., Borkakoti, N., Pravda, L., Tyzack, J. D., Thornton, J. M., & Velankar, S. (2019). Finding enzyme cofactors in Protein Data Bank. Bioinformatics, 35(18), 3510–3511.https: //doi.org/10.1093/bioinformatics/btz115
The script then provides a list of all the heteroatoms found in the json file (txt file), as well as a csv file of all heteroatoms together with their annotation  and metadata. 
Queries https://data.rcsb.org for metadata
Yana Rose, Jose M. Duarte, Robert Lowe, Joan Segura, Chunxiao Bi, Charmi Bhikadiya, Li Chen, Alexander S. Rose, Sebastian Bittrich, Stephen K. Burley, John D. Westbrook. RCSB Protein Data Bank: Architectural Advances Towards Integrated Searching and Efficient Access to Macromolecular Structure Data from the PDB Archive, Journal of Molecular Biology, 2020. DOI: 10.1016/j.jmb.2020.11.003


## **Usage**
The script can be run with the flags 
-d: dictionary file (json)
-c: input file with a list of cofactors/heteroatoms to match (txt). If this is not provided, the output will contains group names from the json list only.

python annotate_heteroatoms.py -d cofactors_dict.json -c output_unique_compounds.txt

## **Example** 
cd example/annotate_heteroatoms
python ../../annotate_heteroatoms.py -d ../../cofactors_dict.json -c ../find_cofactor/summary_output_unique_heteroatoms.txt


## **Output**
1) A csv with all your heteroatoms annotated with the json dictionary, "all_heteroatoms_annotated.csv ". If the heteroatom is not in the list it is annotated as "not a cofactor". The remaining columns contain annotations downloaded from PDB.
2) A list of only the heteroatoms that were cofactors (according to the json dictionary), "cofactor_only.txt".  

