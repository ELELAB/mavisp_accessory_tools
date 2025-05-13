In this example we used the script find_cofactor.py for an index file with mavisp metadata (index.csv) together with a filterfile to only recognize heteroatoms in that list. 


The script uses 4 flags. Mandatory to use either -u or -i
-u: input uniprot AC
-i: input file (csv) mavisp metadata with a column "UniProt AC"
-o: output file prefix. If not provided, default is summary_output
-f: Filtering file (txt). To use if only specific cofactors/heteroatoms should be included in the output. 

# Run for one Uniprot AC
#python ../../find_cofactor.py -u P00533

# Run for file with more than one Uniprot AC
#python ../../find_cofactor.py -i index.csv

# Run with cofactor filter list (to only keep heteroatoms that are in that list)
python ../../find_cofactor.py -i index.csv -f ../annotate_heteroatoms/cofactor_only.txt -o summary_output




