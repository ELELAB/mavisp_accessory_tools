# **Alphafill cofactor survey**

## **Description**  
The `find_cofactor.py` script queries ALphafill for a Uniprot AC or a list of Uniprot ACs. It then extracts the list of potential cofactors for this/these protein(s). 
Only Alphafill results with identity >=30% are included

To run the script, the user must provide a list of uniprot AC's or a single Uniprot AC.
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
-i: input file (txt) list of Uniprot AC's to query
-o: output file prefix. If not provided, default is summary_output
-f: Filtering file (txt). To use if only specific cofactors/heteroatoms should be included in the output. 


# Run for one Uniprot AC
python find_cofactor.py -u P00533

# Run for file with more than one Uniprot AC
python find_cofactor.py -i uniprotACs.txt

# Run with cofactor filter list (to only keep heteroatoms that are in that list)
python find_cofactor.py -i uniprotACs.txt -f cofactors.txt -o summary_output



## **Example**
cd example/find_cofactor
python ../../find_cofactor.py -i uniprotACs.txt -f ../annotate_heteroatoms/cofactor_only.txt -o summary_output

# **Annotate cofactors**
The annotate_cofactors.py script uses a json file to annotate cofactors from the PDB. The json file is a dictionary of cofactors in the PDB, and is compared to a user defined list of heteroatoms. 
The json dictionary can be downloaded from here: https://www.ebi.ac.uk/pdbe/api/pdb/compound/cofactors
This resource was described in: Mukhopadhyay, A., Borkakoti, N., Pravda, L., Tyzack, J. D., Thornton, J. M., & Velankar, S. (2019). Finding enzyme cofactors in Protein Data Bank. Bioinformatics, 35(18), 3510–3511.https: //doi.org/10.1093/bioinformatics/btz115
The script then provides a list of the heteroatoms that were cofactors, as well as a csv file of all heteroatoms together with their annotation. 



## **Usage**
The script can be run with the flags 
-d: dictionary file (json)
-c: input file with a list of cofactors/heteroatoms to match (txt)

python annotate_heteroatoms.py -d cofactors_dict.json -c output_unique_compounds.txt

## **Example** 
cd example/annotate_heteroatoms
python ../../annotate_heteroatoms.py -d ../../cofactors_dict.json -c ../find_cofactor/summary_output_unique_heteroatoms.txt


## **Output**
1) A csv with all your heteroatoms annotated with the json dictionary, "all_heteroatoms_annotated.csv ". If the heteroatom is not in the list it is annotated as "not a cofactor"  
2) A list of only the heteroatoms that were cofactors (according to the json dictionary), "cofactor_only.txt".  

