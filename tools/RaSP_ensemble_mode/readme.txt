#requirements: 
# 1) Scripts prepare_clusters.sh, prepare_25frames.sh, run_analysis.sh
# 2) The directory with PDBs from clustering OR 25_frames_A.pdb


#Cleaning of PDB included in the script:
# for each cluster - 
# 1) Add chain if missing
# 2) Remove solvent and hydrogens
# 3) Fix aminoacid specific names: HIE/HIS, HID/HIS, HIP/HIS, LYN/LYS, ASH/ASP, GLH/GLU, CYX/CYS).
# 4) Ensure that each line contains the correct number of characters. 
# 5) Remove END/TER/ENDMOL lines. 

#To run:
#Run either prepare_clusters.sh OR prepare_25frames.sh depending on the desired input
./prepare_clusters.sh <input_directory> <num_clusters>
# Example: ./prepare_clusters.sh /data/user/shared_projects/mavisp/TPMT/simulations_analysis/free/AF2_18-245/replicate1/CHARMM36/md/6.clustering/ 4
./prepare_25frames.sh <25_frames_A.pdb>
# Example: ./prepare_25frames.sh /data/user/shared_projects/mavisp/TPMT/simulations_analysis/free/AF2_18-245/replicate1/CHARMM36/md/5.ensemble_for_ddg_mutatex/25_frames_A.pdb

#Exit the server
exit
#Re-enter the server
#Activate environment manually
source /usr/local/miniconda3/etc/profile.d/conda.sh
conda activate /usr/local/envs/RaSP_workflow

#Run analysis script 
tsp -N 4 bash run_analysis.sh chain cores frames
#example
#tsp -N 4 bash run_analysis.sh  A 4 26
