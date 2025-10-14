module load python/3.10/modulefile
csv=$1
# change the csv path and the number of mutations you want to plot on the x-axis
# suggested not to use more than 50, otherwise you need to chage the figuresize;
# for more information about available options: python dot_plot.py -h
python piechart.py -i $csv -d  clinvar_interpretation_internal_dictionary.txt -o pie.pdf
python piechart.py -i $csv -d clinvar_interpretation_internal_dictionary.txt  -o pie.png
