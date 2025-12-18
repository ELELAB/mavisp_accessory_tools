module load python

#before running verify with localization.py that your protein is eligible for this module depending on the cellular localization

python localization.py --id uniprotID --fasta uniprotID.fasta --deeploc always 

#have the trimmed wt pdb in the folder

#cp or link the mutation list from cancermuts most recent run
e.g, cp ../cancermuts_saturation/mutlist_04062024.txt .

#link your mutations folder from the mutatex scan for cysteine variants
cp -r /data/raw_data/computational_data/mutatex_data/marinara/pold1/ZN/AF_1-1107/cysteine_pdbs/mutations/AF_POLD1a_1-1107_model0_checked_Repair .

mkdir wt_self 

#mv the wild-type folders with cysteines to a folder called wt_self 
mv AF_POLD1a_1-1107_model0_checked_Repair/CA* wt_self/*

#loss of ss-bond - run 
bash loss-ss_foldx.sh WT_PDB_FILE mutlist.txt

e.g., bash loss_ss_foldx.sh Q07864_1-2286.pdb mutlist_04062024.txt

#denovo SS

bash denovo_ss.sh

