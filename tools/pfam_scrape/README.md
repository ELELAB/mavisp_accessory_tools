# **Pfam Scrape**

## **Description**  
The `pfam_scrape` script automates the extraction and compilation of protein domain
information from MAVISp-generated CSV files. Its primary goal is to generate a dataset
of proteins for benchmarking, by ensuring that only **unique Pfam domain annotations**
are retained, reducing redundancy.

The script allows the user to ensure all processed proteins csv files contain a given
column via command-line arguments.

This script is designed to be executed within:  
`/data/raw_data/computational_data/mavisp_database_saturation/`  
Meaning the script expects the file and folder structure of this directory, however it can be executed anywhere, as specified by option `-p`. 

The script will skip and alert if any protein file does not contain the column 'Pfam domain classification'. 

Users specify the **database version** (i.e., `{base_dir}`) and the **MAVISp mode** (`{mode_dir}`)
via command-line arguments. The script then processes all **final MAVISp CSV protein files** within
the dataset and retrieves the UniProt IDs from **index.csv file**.

If any required files or directories are missing, the script will log an error and exit.  

## **Output**  
The script generates a file called:  **`mavisp_unique_pfam.csv`**  

Each row in this output represents a **trimmed MAVISp protein domain** along with its associated details:  

| Column Name          | Description |
|----------------------|-------------|
| `protein`           | Name of the protein extracted from the filename |
| `mavisp_structure`  | The extracted domain range (e.g., `29-110`) |
| `total_aa`         | The total number of residues in the domain |
| `uniprot_id`       | The corresponding Uniprot Accession ID |
| `pfam_fold`        | The unique Pfam domain annotations for the domain |

### **Filtering & Deduplication**  
- The script **removes redundant Pfam domain annotations** (`pfam_fold`).  
- It retains only **one entry per unique Pfam annotation**, with duplicates being filtered in either:  
  - **Alphabetical order** (if `-a` is specified).  
  - **Random order** (default behavior).  

## **Options**  
```shell
  -h, --help            show this help message and exit
  -p PATH, --path PATH  Path to the base directory.
  -m {simple_mode,ensemble_mode}, --mode {simple_mode,ensemble_mode}
                        MAVISp mode directory to use for parsing deposited CSV files.
  -a, --alphabetical    Flag to drop duplicates by keeping the first domain alphabetically, sorted by protein name. default=Random.
  -c COLUMN, --column COLUMN
                        Option to report only proteins containing a specified column. Write in single quotes ''. Default is all protein CSV are considered in the specified directory
```

## **Requirements**
module load python/3.10/modulefile 

## **Usage**
pfam_scrape -p <base_directory> -m <mode> [-a] [-c 'column_name']

## **Example**
cd example
../pfam_scrape -p ./mavisp_db_date/ -m simple_mode -a
