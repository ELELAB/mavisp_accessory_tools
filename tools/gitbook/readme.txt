#use the template called  config_gitbook.yaml and rename it to your protein of interest:
#e.g. cp config_gitbook.yaml config_ADCK1.yaml

#fill in the fields that applied and remove the ones that are optional and you didn't use - see the example of config_ADCK1.yaml in theis
folder

#run the script to generate the templated Markdown file for Gitbook
module load python


python generate_gitbook_report_v5.py --yaml config_GENENAME.yaml --output GENENAME.md 

#copy in this folder the figures you need to be uploaded from downstream_analysis - remember you might need to customize the dotplot
or lolliplot with your own run of the code 

#open the md file and add the extra text that you want to have and remove what does not apply and shouldn't be shown in the official
gitbook page

#for this inspection (for curators or reviewers, you can for example, Upload your .md file to:

https://stackedit.io
https://dillinger.io
#This shows you what GitBook will approximately render.


#open an issue on github for mavisp_data_collection with the path to this folder so that a first importing of the report at level 0 can be done

#after this when needed reviewers will be appointed and they will work on suggesting changes to be implemented here and a new issue open to move the entry to level 1 or 2 depending on the situation


