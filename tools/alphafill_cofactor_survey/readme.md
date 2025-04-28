# **Alphafill cofactor survey**

## **Description**  
The `find_cofactor.py` script queries ALphafill for a Uniprot AC. It then extracts the list of potential cofactors for this protein. 
Only Alphafill results with identity >=30% are included

The bash script `run_cofactor.sh` runs `find_cofactor.py` for a list of Uniprot AC's and provides a summary output file telling which co-factors were found in each entry.
To run the script, the user must provide a list of uniprot AC's.
The user can adjust the name of the output folder with an optional last argument

### **Output**
| Protein          | Cofactor              |
|------------------|-----------------------|
| `Uniprot AC`     | `cofactor1 cofactor2` |
| `P00325` 	   | `ACN APR CHD CNH CO`  |

The output contains only cofactors found for PDB entries with >=30% identity. 




## **Requirements**
module load python/3.10/modulefile 

## **Usage**
./run_cofactor.sh inputfile.txt outputfile

## **Example**
cd example
../run_cofactor.sh uniprotACs.txt summary_output


