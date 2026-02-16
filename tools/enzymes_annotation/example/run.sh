module load python/3.10/modulefile
uniprot_ac=$1

python ../get_enzyme_catalytic_site.py -ac $uniprot_ac -m
