# **Alphafill cofactor survey**

## **Description**  
The `find_cofactor.py` script queries ALphafill for a Uniprot AC. It then extracts the list of potential cofactors for this protein. 
Only Alphafill results with identity >=30% are included

The bash script `run_cofactor.sh` runs `find_cofactor.py` for a list of Uniprot AC's and provides a summary output file telling which co-factors were found in each entry.
To run the script, the user must provide a list of uniprot AC's.
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
The script uses 3 flags, one of which is mandatory
-i: input file (txt) list of Uniprot AC's to query
-o: output file prefix. If not provided, default is summary_output
-f: Filtering file (txt). To use if only specific cofactors/heteroatoms should be included in the output. 

./run_cofactor.sh -i inputfile.txt -o outputfile -f filterfile.txt

## **Example**
cd example
../run_cofactor.sh -i uniprotACs.txt -o summary_output -f cofactor_only.txt



# **Annotate cofactors**
The annotate_cofactors.py script uses a json file to annotate cofactors from the PDB. The json file is a dictionary of cofactors in the PDB, and is compared to a user defined list of heteroatoms. The script then provides a list of the heteroatoms that were cofactors, as well as a csv file of all heteroatoms together with their annotation. 

## **Usage**
The script can be run with the flags 
-d: dictionary file (json)
-c: input file with a list of cofactors/heteroatoms to match (txt)

python annotate_heteroatoms.py -d cofactors_dict.json -c output_unique_compounds.txt

## **Example** 
cd example
python ../annotate_heteroatoms.py -d ../cofactors_dict.json -c summary_output_unique_heteroatoms.txt

## **Output**
1) A csv with all your heteroatoms annotated with the json dictionary, "all_heteroatoms_annotated.csv ". If the heteroatom is not in the list it is annotated as "not a cofactor"  
2) A list of only the heteroatoms that were cofactors (according to the json dictionary), "cofactor_only.txt".  
