# Notice that a few Docker containers need to be active for this to work
# please check instructions in the README file, in the main directory of the repository

. /usr/local/envs/spliceai/bin/activate

../../splice_lookup -i maf_input.csv -c config_maf.yaml -f -g /data/databases/genome_annotation/ -d 500 -t 9
