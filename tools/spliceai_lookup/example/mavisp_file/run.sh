#delay 20s
source /data/user/marnaudi/spliceai_lookup/spliceai-env/bin/activate
conda activate /home/marnaudi/.conda/envs/bioenv/
python ../../splice_lookup.py -i mavisp_input.csv -m -g /data/databases/genome_annotation/ -d 500 -t 9
