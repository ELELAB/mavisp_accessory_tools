# set ENST id of interest in this variable - to be changed
enst=$1

# link the AF aminoacids substitution scores file in this directory
ln -s /data/databases/alphamissense/AlphaMissense_isoforms_aa_substitutions.tsv.gz .

# write header to output file
zcat AlphaMissense_isoforms_aa_substitutions.tsv.gz | head -n 4 | tail -n 1 > am.tsv

# find lines containing UniProt AC of interest and save them to file
zgrep -P "${enst}\t" AlphaMissense_isoforms_aa_substitutions.tsv.gz >> am.tsv

# compress output file
gzip am.tsv
