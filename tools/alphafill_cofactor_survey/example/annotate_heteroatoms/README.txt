# **Annotate cofactors**
The annotate_cofactors.py script uses a json file to annotate cofactors from the PDB. The json file is a dictionary of cofactors in the PDB, and is compared to a user defined list of heteroatoms. 
The json dictionary can be downloaded from here: https://www.ebi.ac.uk/pdbe/api/pdb/compound/cofactors
This resource was described in: Mukhopadhyay, A., Borkakoti, N., Pravda, L., Tyzack, J. D., Thornton, J. M., & Velankar, S. (2019). Finding enzyme cofactors in Protein Data Bank. Bioinformatics, 35(18), 3510–3511.https: //doi.org/10.1093/bioinformatics/btz115
The script then provides a list of all the heteroatoms found in the json file (txt file), as well as a csv file of all heteroatoms together with their annotation.



## **Usage**
The script can be run with the flags 
-d: dictionary file (json)
-c: input file with a list of cofactors/heteroatoms to match (txt). If this is not provided, the output will contains group names from the json list only.

python annotate_heteroatoms.py -d cofactors_dict.json -c output_unique_compounds.txt

## **Example** 
#Here we run the script to annotate cofactors that were listed with the find_cofactor.py script. 
python ../../annotate_heteroatoms.py -d ../../cofactors_dict.json -c ../find_cofactor/summary_output_unique_heteroatoms.txt

