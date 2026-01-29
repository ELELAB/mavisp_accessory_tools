This script implements the Gaussian Mixture Model (GMM) procedure used in the ProteoCast framework (DOI: 10.1101/2025.02.09.637326) to classify GEMME evolutionary scores into three functional categories:

neutral,mild, impactful

It is written to use in input a MAVISp csv file including GEMME scores
module load python
python gemme_GMM.py -i POLE-ensemble_mode.csv -o POLE-ensemble_mode_GMM.csv  --random-state 42
