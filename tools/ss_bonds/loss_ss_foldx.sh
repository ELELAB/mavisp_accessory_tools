module load python
python  ss_bonds_4_parallel-foldx.py wt $1  \
  --out-csv native_disulfides_foldx.csv \
  --ss-loose-min 1.8 --ss-loose-max 3.5 \
  --variants-file $2 \
  --debug \
  --foldx-wt wt_self \
  --foldx-mode LOOSE \
  --native-source foldx \
  --out-mutlist mutlist_disulfide_disruptions_foldx.tsv
