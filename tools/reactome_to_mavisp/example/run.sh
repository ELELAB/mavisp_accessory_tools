source /data/user/marnaudi/reactome/reactome2py/bin/activate
python ../reactome_to_mavisp.py -uf uniprot_list.txt -o reactome_outputs -s -r 
python ../reactome_post_process.py -i reactome_outputs -o reactome_outputs/summary