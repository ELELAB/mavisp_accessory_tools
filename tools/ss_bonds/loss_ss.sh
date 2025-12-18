module load python
python  ss_bonds_4_parallel.py wt $1  \
  --out-csv native_disulfides.csv \
  --native-mode LOOSE --ss-loose-min 1.8 --ss-loose-max 3.5 --no-angle \
  --variants-file $2 \
    --rotamer-scan \
  --rot-threshold 10.0 --dump-pairs all_pairs.csv --debug \
  --out-mutlist mutlist_disulfide_disruptions.tsv
