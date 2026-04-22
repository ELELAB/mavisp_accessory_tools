#take mutlist for ddg2summary for only PTMs
cp ../../cancermuts_saturation/mutlist_mutatex_P*.txt .
#symbolic links to mutatex scans: pdb file checked, final_averages dir
run ddg2summary
bash run_ddgs.sh Q7Z695_88-626_model0_checked.pdb mutlist_mutatex_P.txt
