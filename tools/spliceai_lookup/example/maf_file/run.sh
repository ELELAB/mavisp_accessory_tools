#delay 20s
source /data/user/marnaudi/spliceai_lookup/spliceai-env/bin/activate
conda activate /home/marnaudi/.conda/envs/bioenv/
python ../../splice_lookup.py -i maf_input.csv -c maf_config.yaml  -f -g /data/databases/genome_annotation/ -d 500 -t 9
