# do not change - AF installation dir
AF=/usr/local/alphafold

# do not change - AF database dir
DATA_DIR=/scratch/databases/alphafold 

# input fasta file
# fasta files should have all the sequences that we want to
# fold/bind in the same run
FASTA=input.fasta

# output directory. Default is "here" but modify accordingly
OUTPUT_DIR=$(pwd) 

# maximum date to consider for templates. Format is YYYY-MM-DD
# this should be an older date if we want to exclude some templates
# i.e. if we want to model a structure without alphafold using its
# corresponding experimental structure(s). Otherwise it should be
# today or a date in the future
MAX_TEMPLATE_DATE=2023-03-17

# do not change - activate virtualenv
. $AF/venv/bin/activate

# run
python3 $AF/docker/run_docker.py \
  --fasta_paths=./$FASTA \
  --max_template_date=$MAX_TEMPLATE_DATE \
  --model_preset=multimer \
  --db_preset=full_dbs \
  --data_dir=$DATA_DIR \
  --output_dir=$OUTPUT_DIR

# notice that you can select which GPU(s) to run on by adding the
# --gpu-devices=0,1 option. In this case, using "0,1" selects both GPUs,
# or you can select only one using the appropriate option
