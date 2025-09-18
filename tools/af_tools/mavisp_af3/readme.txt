# run with AF3: https://alphafoldserver.com/fold/4e0f61a712b82b87
# date_of_run: 05072024  ###PLEASE COMPLETE THIS

# this template is designed to work on the AlphaFold 3 output downloaded from
# https://alphafoldserver.com . The data is usually available in the form of a
# zip file containing a fold_* directory. Download the file from the server and
# copy it in the current folder. 

# notice that the template might not work correctly if other
# files are present, especially files from a previous
# run of the same template.

# run:

bash get_scores.sh fold_*.zip

# if you have a model complex, you might want to calculate the PdockQ2 score
# by running

bash get_pdockq2.sh

# finally, clean up

bash clean.sh

# notes
# for conversion to PDB (if needed for future purpouses)
# BeEM my_cif.cif -o my_pdb_file_name

