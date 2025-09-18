mkdir unrelaxed_models
mv unrelaxed*pdb unrelaxed_models
tar -cvzf unrelaxed_models.tar.gz unrelaxed_models
mkdir relaxed_models
mv relaxed*pdb relaxed_models
tar -cvzf relaxed_models.tar.gz relaxed_models
cp /data/raw_data/computational_data/alphafold_data/templates/* .
rm -r relaxed_models unrelaxed_models
. /usr/local/envs/py37/bin/activate
./extract_pae . -n 25 -v -a AF2
mkdir pae_ranked
mv pae_ranked* pae_ranked
rm *.pkl
tar -cvzf msas
