rmsd_matrix.py

Computes all-atom and Cα RMSD between the initial PDB structures used in the
FoldX5 and FoldX5.1 datasets. Structures are matched by UniProt ID and protein
name. RMSD is computed only on overlapping residues when the two versions cover
different regions of the same protein. Alignment is performed using the QCP
algorithm via MDAnalysis.


REQUIREMENTS
Python 3.10
MDAnalysis >= 2.7.0
numpy, pandas, matplotlib, seaborn

USAGE
python rmsd_matrix.py \
    -f /path/to/foldx5_initial_structures \
    -i /path/to/data_collection_foldx5.1 \
    -o ./rmsd_results


OUTPUT
------
rmsd_results.csv                table with all-atom and Cα RMSD per structure pair
rmsd_distributions.png          histogram of RMSD values with mean and median lines
rmsd_scatter_allatom_vs_ca.png  scatter of all-atom vs Cα RMSD per structure
rmsd_barplot_sorted.png         bar plot of all structures sorted by Cα RMSD
rmsd_violin.png                 violin plot of RMSD distribution across all protein
