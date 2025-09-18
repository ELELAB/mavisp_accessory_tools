set -euo pipefail

mkdir relaxed_models
mv *pdb relaxed_models
tar cvzf relaxed_models.tar.gz relaxed_models
rm -r relaxed_models

mkdir pae_ranked
mv pae_ranked*.* pae_ranked

rm *.json *.cif terms_of_use.md

