In this example we used the script find_cofactor.py for a list of UniprotACs (txt) together with a filterfile to only recognize heteroatoms in that list. 


The script uses 4 flags. Mandatory to use either -u or -i
-u: input uniprot AC
-i: input file (txt) list of Uniprot AC's to query
-o: output file prefix. If not provided, default is summary_output
-f: Filtering file (txt). To use if only specific cofactors/heteroatoms should be included in the output. 

# Run for one Uniprot AC
#python ../../find_cofactor.py -u P00533

# Run for file with more than one Uniprot AC
#python ../../find_cofactor.py -i uniprotACs.txt

# Run with cofactor filter list (to only keep heteroatoms that are in that list)
python ../../find_cofactor.py -i uniprotACs.txt -f ../annotate_heteroatoms/cofactor_only.txt -o summary_output




