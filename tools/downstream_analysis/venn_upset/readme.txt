#have first the mavisp csv as a local copy in your folder ../mavisp_csv/date/
#do a symbolic link for the file so that we keep track of the last one used in this folder
ln -s ../mavisp_csv/15062023/ADCK1-simple_mode.csv 
Rscript venn_upset.R ADCK1-simple_mode.csv
