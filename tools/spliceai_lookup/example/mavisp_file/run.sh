# Notice that a few Docker containers need to be active for this to work
# please check instructions in the README file, in the main directory of the repository

. /usr/local/envs/spliceai/bin/activate

../../splice_lookup -i mavisp_input.csv -m -g /data/databases/genome_annotation/ -d 500 -t 9
