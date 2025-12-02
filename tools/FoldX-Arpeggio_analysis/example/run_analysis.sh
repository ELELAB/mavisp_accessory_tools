module load python/3.10/modulefile
cp -r ../Snakefile .
cp -r ../templates .
snakemake -s Snakefile -c 1
