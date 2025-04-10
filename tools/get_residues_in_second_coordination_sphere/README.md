# get_second_sphere_residues.py 
The get_second_sphere_residues.py is a Python script designed to identify the second coordination sphere residues of catalytic residues provided as input.

## Requirements

### language

The script requires Python 3.10

### libraries

The script requires the following libraries
- Open Babel >= 3.0.0
- Biopython >= 1.80
- gemmi >= 0.5.8
- os
- re
- pandas 
- json
- argparse
- yaml
- logging

## Description

This script processes structural data to identify interaction networks for specified residues within a protein structure, leveraging the Arpeggio software for contact analysis. The script requires a PDB file containing the target structure, a CSV file listing residues of interest, and a configuration file defining parameters for the Arpeggio run. Initially, it validates the residues specified in the CSV, ensuring they match the structure in the PDB file, and checks the configuration file for correct formatting.


Once validated, the script processes the PDB file by converting it into a CIF file compatible with the Arpeggio run. It then performs an Arpeggio analysis for each target residue, identifying all interacting residues within a specified cutoff distance as outlined in the configuration. The contacts are categorized by interaction type, and the results are saved in an organized directory as JSON files.


At the conclusion of the analysis, the pipeline generates a CSV file that catalogs the target residues, their interacting residues, the specific atoms involved in each contact, and the nature of each interaction. Prior to finalization, the CSV is preprocessed to exclude clashing contacts, which are attributed to poorly modeled residues; these are instead listed in a separate file, clashes_contact.csv. Additionally, unless the -k flag is specified, contacts classified as "proximal" are removed, as these are deemed non-critical according to the Arpeggio framework. The residues of interest, along with their interactors, are summarized in a text file named catalytic_and_second_sphere_residues.txt. Furthermore, the analysis produces an additional file, catalytic_and_second_sphere_residues_saturation.txt, which enumerates all possible mutations associated with the target residues and their interactors identified during the analysis.


## Input
The script requires the following inputs, specified as flags:

**-p,--pdb_file**, Path to the PDB file containing the structure for which interacting residues will be identified.


**-cr,--catalytic_residues** CSV file specifying residues of interest, with chain ID and residue position.


**-o,--output** Output file name. Default is second_sphere_catalytic_residues.csv.


**-k,--keep_proximal** keep the 'proximal' contacts in the output file. Default (they are removed from the output).


**-y,--config_file** Specifies the configuration file in YAML format, which defines Arpeggio parameters and includes the chem_comp field to be added to the CIF file.


The catalytic_residues file specifies the target residues for which the script should calculate contacts. This file must follow the structure below:

|catalytic_residue|chain_id|
|-----------------|--------|
|Thr179|A|
|Asn146|A|

where:

**catalytic_residue**: Indicates the residue of interest, including the residue name and position (e.g., "Thr179". the nomenclature is three letters).


**chain_id**: Specifies the chain identifier for the residue, as given in the input PDB file.


An example of input file is the following:

```
catalytic_residue,chain_id
Thr179,A
Lys143,A
Asn146,A
Asp159,A
Asp141,A
Asp169,A
```
The configuration file in yaml format, must specify the following parameters to customize the Arpeggio run and the chem_comp parameter essential for the run to add in the cif file:

**write_hydrogenated** Boolean; if True, writes an MMCIF file with added hydrogens


**minimise_hydrogens** Boolean; if True, performs energy minimization of added hydrogens using OpenBabel.


**minimisation_steps** Integer; specifies the number of minimization steps for hydrogen optimization.


**minimisation_forcefield** String; forcefield for hydrogen minimization, options are "MMFF94", "UFF", and "Ghemical" (note: the latter is not recommended).


**minimisation_method** String; specifies the minimization method, with options "DistanceGeometry", "SteepestDescent", and "ConjugateGradients" (recommended).


**vdw_comp** Float; compensation factor for Van der Waals radius-based interactions.


**interacting_cutoff** Float; distance cutoff (in Å) for identifying interacting residues.


**ph** Float; pH level for hydrogen addition.


**include_sequence_adjacent** Boolean; if True, includes non-bonding interactions between sequential residues (default: False).


**use_ambiguities** Boolean; enables ambiguity handling for interactions with ambiguous contacts.


**output** String; specifies the directory for storing Arpeggio output files.


**chem_comp** String multiline; specify the chem_comp parameters to be added to the CIF file.



Note: The script uses Biopython to convert the PDB file to a CIF format. However, the library does not add the chem_comp field, which is essential for Arpeggio to run properly. For this reason, the script performs a post-processing step to add this field. The chem_comp field contains information such as the name, formula, molecular weight, etc., of the residues present in the PDB file.



The configuration file provided with the script includes all necessary information for the 20 natural amino acids, as well as the following cofactors: ZN, MG, CA, ADP, ATP, DA, DT, DC, DG, and H2O.



If a PDB file contains any residues or atoms not specified in the chem_comp field, please update the configuration file accordingly.



here an example of config file

```
write_hydrogenated: False
minimise_hydrogens: False
minimisation_steps: 50
minimisation_forcefield: "MMFF94"
minimisation_method: "ConjugateGradients"
vdw_comp: 0.1
interacting_cutoff: 5.0
ph: 7.4           
include_sequence_adjacent: False
use_ambiguities: False
output: "arpeggio_output"
chem_comp: |
  #
  loop_
  _chem_comp.formula
  _chem_comp.formula_weight
  _chem_comp.id
  _chem_comp.mon_nstd_flag
  _chem_comp.name
  _chem_comp.pdbx_synonyms
  _chem_comp.type
  "C3 H7 N O2"         89.093  ALA y ALANINE                              ? "L-PEPTIDE LINKING"
  "C6 H15 N4 O2"       175.209 ARG y ARGININE                             ? "L-PEPTIDE LINKING"
  "C4 H8 N2 O3"        132.118 ASN y ASPARAGINE                           ? "L-PEPTIDE LINKING"
  "C4 H7 N O4"         133.103 ASP y "ASPARTIC ACID"                      ? "L-PEPTIDE LINKING"
  "C3 H7 N O2 S"       121.158 CYS y CYSTEINE                             ? "L-PEPTIDE LINKING"
  "C5 H10 N2 O3"       146.144 GLN y GLUTAMINE                            ? "L-PEPTIDE LINKING"
  "C5 H9 N O4"         147.129 GLU y "GLUTAMIC ACID"                      ? "L-PEPTIDE LINKING"
  "C2 H5 N O2"         75.067  GLY y GLYCINE                              ? "PEPTIDE LINKING"
  "C6 H10 N3 O2"       156.162 HIS y HISTIDINE                            ? "L-PEPTIDE LINKING"
  "C6 H13 N O2"        131.173 ILE y ISOLEUCINE                           ? "L-PEPTIDE LINKING"
  "C6 H13 N O2"        131.173 LEU y LEUCINE                              ? "L-PEPTIDE LINKING"
  "C6 H15 N2 O2"       147.195 LYS y LYSINE                               ? "L-PEPTIDE LINKING"
  "C5 H11 N O2 S"      149.211 MET y METHIONINE                           ? "L-PEPTIDE LINKING"
  "C9 H11 N O2"        165.189 PHE y PHENYLALANINE                        ? "L-PEPTIDE LINKING"
  "C5 H9 N O2"         115.130 PRO y PROLINE                              ? "L-PEPTIDE LINKING"
  "C3 H7 N O3"         105.093 SER y SERINE                               ? "L-PEPTIDE LINKING"
  "C4 H9 N O3"         119.119 THR y THREONINE                            ? "L-PEPTIDE LINKING"
  "C11 H12 N2 O2"      204.225 TRP y TRYPTOPHAN                           ? "L-PEPTIDE LINKING"
  "C9 H11 N O3"        181.189 TYR y TYROSINE                             ? "L-PEPTIDE LINKING"
  "C5 H11 N O2"        117.146 VAL y VALINE                               ? "L-PEPTIDE LINKING"
  "Zn 2"               65.409  ZN  . "ZINC ION"                           ? "non-polymer"
  "MG 2"               24.305  MG  . "MAGNESIUM ION"                      ? "non-polymer"
  "Ca 2"               40.078  CA  . "CALCIUM ION"                        ? "non-polymer"
  "H2 O"               18.015  HOH . "WATER"                              ? "non-polymer" 
  "C10 H16 N5 O13 P3"  507.181 ATP . "ADENOSINE-5'-TRIPHOSPHATE"          ? "non-polymer"        
  "C10 H15 N5 O10 P2"  427.201 ADP n "ADENOSINE-5'-DIPHOSPHATE"           ? "non-polymer"
  "C10 H14 N5 O6 P"    331.222 DA  y "2'-DEOXYADENOSINE-5'-MONOPHOSPHATE" ? "DNA linking"
  "C9 H14 N3 O7 P"     307.197 DC  y "2'-DEOXYCYTIDINE-5'-MONOPHOSPHATE"  ? "DNA linking"
  "C10 H14 N5 O7 P"    347.221 DG  y "2'-DEOXYGUANOSINE-5'-MONOPHOSPHATE" ? "DNA linking"
  "C10 H15 N2 O8 P"    322.208 DT  y "THYMIDINE-5'-MONOPHOSPHATE"         ? "DNA linking"  
  #
```

## Output

The script generate the cif file from the pdb in the same location where the script has been run. Additionally it stores the Arpeggio output files in the directory specified in the configuration file. The output file is a JSON file containing a list of all calculated contacts. Below is an example of the output for a specific contact:

```
{
        "bgn": {
            "auth_asym_id": "A",
            "auth_atom_id": "NZ",
            "auth_seq_id": 143,
            "label_comp_id": "LYS",
            "label_comp_type": "P",
            "pdbx_PDB_ins_code": " "
        },
        "contact": [
            "vdw_clash",
            "ionic",
            "polar"
        ],
        "distance": 2.66,
        "end": {
            "auth_asym_id": "A",
            "auth_atom_id": "OD2",
            "auth_seq_id": 141,
            "label_comp_id": "ASP",
            "label_comp_type": "P",
            "pdbx_PDB_ins_code": " "
        },
        "interacting_entities": "INTRA_NON_SELECTION",
        "type": "atom-atom"
    },
```
In this example, the "bgn" and "end" fields specify the two residues under investigation.


The script automatically parses the Arpeggio output, generating a CSV file that collects the contacts for each residue specified in the input file. The output file has the following structure:

|catalytic_residue|catalytic_residue_atom|second_sphere_residue|second_sphere_residue_atom|contact|distance|interacting_entities|type|
|-----------------|----------------------|---------------------|--------------------------|-------|--------|--------------------|----|
|ASN146|CB|ASP141|O|proximal|3.94|INTER|atom-atom|
|ASN146|CG|ASP141|C|proximal|4.93|INTER|atom-atom|
|ASN146|CG|ASP141|CG|proximal|4.54|INTER|atom-atom|
|ASN146|CG|ASP141|O|proximal|3.73|INTER|atom-atom| 

where:

**catalytic_residue** The target residue specified in the input file.


**catalytic_residue_atom** The atom of the target residue establishing a contact.


**second_sphere_residue** The residue identified in contact with the target residue.


**second_sphere_residue_atom** The atom of the residue in contact with the target one.


**contact** Indicates the type of contact.


**distance** The distance between the two atoms.


**interacting_entities** The type of entities that are interacting.


**type** The type of interactions.



The script generates a txt file, "catalytic_and_second_sphere_residues.txt", containing the residues specified as input and their interactors

```
F172
D169
E48
D141
K143
A145
```
The script generates an additional txt file, "catalytic_and_second_sphere_residues_saturation_mutlist.txt", listing all the possible mutations associated with the residues specified as input and their interactors.

```
L225A
L225C
L225D
L225E
L225F
L225G
L225H
L225I
L225M
L225N
L225P
L225Q
L225R
L225S
L225T
L225V
L225Y
L225W
L225K
L225s
L225y
L225p
I166A
```

### Type of contacts ("contact" and "type" columns)

**atom-atom contact**

|Key|Interaction|Description|
|---|-----------|-----------|
|clash|Clash|Denotes if the atom is involved in a steric clash.|
|covalent|Covalent|Denotes if the atom appears to be covalently bonded.|
|vdw_clash|VdW|Clash|Denotes if the van der Waals radius of the atom is clashing with one or more other atoms.|
|vdw|VdW|Denotes if the van der Waals radius of the the atom is interacting with one or more other atoms.|
|proximal|Proximal|Denotes if the atom is > the VdW interaction distance, but within 5 Angstroms of other atom(s).|
|hbond|Hydrogen Bond|Denotes if the atom forms a hydrogen bond.|
|weak_hbond	Weak Hydrogen Bond	Denotes if the atom forms a weak hydrogen bond.|
|xbond|Halogen Bond|Denotes if the atom forms a halogen bond.|
|ionic|Ionic|Denotes if the atom may interact via charges.|
|metal|Metal Complex|Denotes if the atom is part of a metal complex.|
|aromatic|Aromatic|Denotes an aromatic ring atom interacting with another aromatic ring atom.|
|hydrophobic|Hydrophobic|Denotes hydrophobic interaction.|
|carbonyl|Carbonyl|Denotes a carbonyl-carbon:carbonyl-carbon interaction.|
|polar|Polar|Less strict hydrogen bonding (without angle terms).|
|weak_polar|Weak Polar|Less strict weak hydrogen bonding (without angle terms).|


**atom-plane interactions**


|Key|Interaction|Description|
|---|-----------|-----------|
|CARBONPI|Carbon-PI|Weakly electropositive carbon atom - Π interactions|
|CATIONPI|Cation-PI|Cation - Π interactions|
|DONORPI|Donor-PI|Hydrogen Bond donor - Π interactions|
|HALOGENPI|Halogen-PI|Halogen bond donors - Π|
|METSULPHURPI|Sulphur-PI|Methionine sulphur - Π ring interactions [ref]|


**plane-plane interactions**


Follows nomenclature established by https://www.sciencedirect.com/science/article/pii/S0079610707000442?via%3Dihub


**group-group/plane interactions**



|Key|Interaction|Description|
|---|-----------|-----------|
|AMIDEAMIDE|amide - amide|https://onlinelibrary.wiley.com/doi/10.1002/jcc.21212|
|AMIDERING|amide - ring|https://doi.org/10.1016/0014-5793(86)80730-X|


### Interacting entities

|Key|Meaning|
|---|-------|
|INTER|Between an atom from the user's selection and a non-selected atom|
|INTRA_SELECTION|Between two atoms both in the user's selection|
|INTRA_NON_SELECTION|Between two atoms that are both not in the user's selection|
|SELECTION_WATER|Between an atom in the user's selection and a water molecule|
|NON_SELECTION_WATER|Between an atom that is not in the user's selection and a water molecule|
|WATER_WATER|Between two water molecules|

for more information regarding the contacts please consult https://github.com/PDBeurope/arpeggio

The script generates an additional csv file called "clashes_contacts.csv" collecting all the clashing contacts of the residues specified as input and their interactors:

```
|catalytic_residue|catalytic_residue_atom|second_sphere_residue|second_sphere_residue_atom|contact|distance|interacting_entities|type|
|-----------------|----------------------|---------------------|--------------------------|-------|--------|--------------------|----|
|ASN146|ND2|ASP141|O|vdw_clash_polar|2.81|INTER|atom-atom|
|ASN146|ND2|ASP141|OD2|vdw_clash_polar|2.91|INTER|atom-atom|
|ASN146|O|GLY158|N|vdw_clash_polar|2.89|INTER|atom-atom|
|ASN146|N|LYS143|O|vdw_clash_polar|2.94|INTER|atom-atom|
```

Here an exmaple of output tree:
```
├── P51955.cif
├── P51955.pdb
├── arpeggio_output
│   ├── P51955_ASN_146.json
│   ├── P51955_ASP_141.json
│   ├── P51955_ASP_159.json
│   ├── P51955_ASP_169.json
│   ├── P51955_LYS_143.json
│   └── P51955_THR_179.json
├── catalytic_residues.csv
├── catalytic_and_second_sphere_residues.txt
├── catalytic_and_second_sphere_residues_saturation_mutlist.txt 
├── clashes_contacts.csv 
├── run.sh
├── second_sphere_catalytic_residues.csv
```

## Usage

```
module load python/3.10/modulefile
python ../get_second_sphere_residues.py -p P51955.pdb -cr catalytic_residues.csv -y ../config.yaml 
```

In the example folder in order to reproduce the analysis run:

``
bash run.sh
``




