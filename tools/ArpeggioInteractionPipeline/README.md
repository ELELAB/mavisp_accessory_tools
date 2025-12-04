# ArpeggioInteractionPipeline.py
Arpeggio Interaction Pipeline is a Python script for analyzing atomic interactions in protein and molecular structures using Arpeggio. It supports conversion of PDB files to CIF format, validation of residues from input files, computation of intra- and inter-molecular contacts, filtering of low-quality or clashing contacts, and generation of single-point or saturation mutation lists for residues involved in interactions. The script can run on entire structures, specific residues, or defined residue pairs, and produces structured CSV and JSON outputs along with mutation lists.

## Requirements

### language

The script requires Python >= 3.8

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
- sys
- arpeggio


## Description

This script processes structural data to identify interaction networks among all residues and additional molecular entities, such as cofactors or ligands, within a protein structure. Alternatively, it can focus on specific residues or molecules provided through dedicated CSV files. The analysis leverages the Arpeggio software to compute atomic contacts and classify them by interaction type.
The script requires a PDB file containing the target structure, a YAML configuration file defining parameters for the Arpeggio run, and—if a focused analysis is desired—one or more CSV files listing the residues or residue pairs of interest. When CSV files are provided, the script first validates the listed residues or molecules, ensuring that they match those present in the input PDB file, and checks the configuration file for proper formatting.
After validation, the PDB file is converted into a CIF file compatible with Arpeggio, and the analysis proceeds in one of three modes:

1) **Full analysis**, Computes all intra- and inter-chain contacts across the entire structure, classifying them by interaction type and saving the results as organized JSON files. Upon completion, the pipeline generates a comprehensive CSV file listing all identified interactions, including information such as residue name, chain ID, position, atom names, and interaction type. The column names can be customized through a dedicated section in the configuration file for improved readability.

2) **Single-molecule analysis**, Uses a CSV file specifying individual residues or molecules. For each, it identifies all interacting partners within the defined cutoff distance, categorizes the contacts by interaction type, and stores them as JSON files. Upon completion, the pipeline generates a comprehensive CSV file listing all identified interactions, including information such as residue name, chain ID, position, atom names, and interaction type. The column names can be customized through a dedicated section in the configuration file for improved readability.


3) **Pairwise analysis** Uses a CSV file specifying residue or molecule pairs. For each pair, the interactions between the two entities are computed, categorized, and saved as JSON files.
Upon completion, the pipeline generates a comprehensive CSV file listing all identified interactions, including information such as residue name, chain ID, position, atom names, and interaction type. The column names can be customized through a dedicated section in the configuration file for improved readability.

Before finalization, clashing contacts—typically arising from poorly modeled residues are excluded and listed separately in clashes_contact.csv. The -i flag can be used to restrict the output to inter-chain (interface) contacts only, while the -k flag allows retention of “proximal” contacts, which are otherwise filtered out as non-significant according to the Arpeggio framework.
Finally, the pipeline produces two additional summary files:
interacting_entities_list.txt, listing all residues and their interacting partners.
interacting_entities_saturation_mutlist.txt, enumerating all possible mutations for the residues identified during the analysis.
These files will not include molecules different from protein residues.


## Input
The script requires the following inputs, specified as flags:

**-p,--pdb_file**, Path to the PDB file containing the structure for which interacting residues will be identified.


**-i, --interface**, Keep only inter-chain interactions (interface between different chain IDs).


**-sm, --single_molecule_file**, CSV for 'single_molecule' mode (first_molecule,first_molecule_position,first_molecule_chain_id).


**-pm, --pairs_of_molecules_file**, CSV for 'pairs_of_molecules' mode (first_...,second_...).


**-o,--output** Output file name. Default is arpeggio_output.csv.


**-k,--keep_proximal** keep the 'proximal' contacts in the output file. Default (they are removed from the output).


**-y,--config_file** Specifies the configuration file in YAML format, which defines columns names in the output csv,the Arpeggio parameters and includes the chem_comp field to be added to the CIF file.

N.B: The -sm and -pm flags are mutually exclusive. If neither of them is specified, the analysis is automatically performed on all molecular entities present in the structure.

### single molecule analysis

A CSV file specifying the target molecular entities for which the script should calculate contacts can be provided to perform an analysis focused on specific molecules within the cutoff distance defined in the configuration file.
This file must follow the structure below:

|first_molecule|first_molecule_position|first_molecule_chain_id|
|--------------|-----------------------|-----------------------|
|ATP|401|A|
|GLN|100|A|

where:

**first_molecule**: The molecule of interest. Both one-letter and three-letter formats are supported (e.g., Thr, THR, and T are all valid).


**first_molecule_position**: The residue or molecule position, as given in the input PDB file.


**first_molecule_chain_id**: The chain identifier of the molecule, as given in the input PDB file.



### pairs of molecules analysis

A CSV file specifying pairs of molecular entities can be provided to perform a focused analysis on specific interactions between two molecules.
This file must follow the structure below:

|first_molecule|first_molecule_position|first_molecule_chain_id|second_molecule|second_molecule_position|second_molecule_chain_id|
|--------------|-----------------------|-----------------------|---------------|------------------------|------------------------|
|ASP|699|B|THR|703|B|
|Gln|700|B|THR|491|A|

where:

**first_molecule**: The first molecule of interest (one-letter and three-letter formats supported).


**first_molecule_position**: Position of the first molecule in the PDB file.


**chain_id**: Chain identifier of the first molecule.


**second_molecule**: The second molecule of interest (one-letter and three-letter formats supported).


**second_molecule_position**: Position of the second molecule in the PDB file.


**second_molecule_chain_id**: Chain identifier of the second molecule


### Configuration file

The configuration file, provided in YAML format, defines both the custom output column names for the final CSV and the Arpeggio runtime parameters, including the essential chem_comp definitions required for proper execution.

#### Output columns customization

**col_input_residue**, Column name for the target molecule listed in the CSV file (useful when performing single-molecule or pairwise analyses; e.g., “catalytic_residue”).

**col_input_residue_pos**, Column name for the target molecule’s position.

**col_input_residue_chain_id**, Column name for the target molecule’s chain ID.

**col_input_residue_atom**, Column name for the target atom involved in the contact.

**col_residue_in_contact**, Column name for the molecule interacting with the target molecule.

**col_residue_in_contact_atom**, Column name for the atom of the interacting molecule.

**col_residue_in_contact_pos**, Column name for the position of the interacting molecule

**col_residue_in_contact_chain_id**, Column name for the chain ID of the interacting molecule.

#### Arpeggio run parameters

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
output_columns_customization:
    col_input_residue: "catalytic_residue"
    col_input_residue_pos: "catalytic_residue_pos"
    col_input_residue_chain_id: "catalytic_residue_chain_id"
    col_input_residue_atom: "catalytic_residue_atom"
    col_residue_in_contact: "second_sphere_residue"
    col_residue_in_contact_atom: "second_sphere_residue_atom"
    col_residue_in_contact_pos: "second_sphere_pos"
    col_residue_in_contact_chain_id: "second_sphere_residue_chain_id"
arpeggio_run_parameters:
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

The script automatically parses the Arpeggio output and generates a CSV file summarizing all contacts identified for each residue specified in the input file.
The column names are defined according to the entries provided in the configuration file.
The output file has the following structure:

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


The script automatically parses the Arpeggio output, generating a CSV file that collects the contacts for each residue specified in the input file. The output file has the following structure with the columns name called after as speciifed in the config file:

|catalytic_residue|catalytic_residue_pos|catalytic_residue_chain_id|catalytic_residue_atom|second_sphere_residue|second_sphere_pos|second_sphere_residue_atom|second_sphere_residue_chain_id|contact|distance|interacting_entities|type|
|-----------------|---------------------|--------------------------|----------------------|---------------------|-----------------|--------------------------|------------------------------|-------|--------|--------------------|----|
|THR|703|B|N|ASP|699|O|B|"vdw, hbond, polar"|3.11|INTRA_SELECTION|atom-atom|
|THR|703|B|CB|ASP|699|O|B|"proximal, weak_polar"|3.47|INTRA_SELECTION|atom-atom|
|THR|703|B|CG2|ASP|699|O|B|"vdw, weak_polar"|3.26|INTRA_SELECTION|atom-atom|

where:

**catalytic_residue** The target molecule specified in the input file, or automatically identified during the full analysis performed on all entities in the structure.


**catalytic_residue_pos** The position of the target molecule, as provided in the input file or determined during the full analysis.


**catalytic_residue_chain_id**  The chain identifier of the target molecule specified in the input file or identified in the full structure analysis.


**catalytic_residue_atom** The atom of the target molecule that establishes the contact.


**second_sphere_residue** The molecule identified as being in contact with the target molecule.


**second_sphere_residue_pos** The position of the molecule identified in contact with the target molecule.


**second_sphere_residue_chain_id** The chain identifier of the molecule found in contact with the target molecule.


**second_sphere_residue_atom** The atom of the interacting molecule involved in the contact.


**contact** Indicates the type of contact.


**distance** The distance between the two atoms.


**interacting_entities** The type of molecular entities that are interacting (e.g., INTRA_SELECTION, INTER_CHAIN).


**type** The overall interaction category (e.g., atom–atom).



The script also generates a text file, interacting_entities_list.txt, listing each target molecule (or residue) and all its identified interactors.

```
F172
D169
E48
D141
K143
A145
```
The script generates an additional txt file, "interacting_entities_saturation_mutlist.txt", listing all the possible mutations associated with the residues specified as input and their interactors.

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
python ../ArpeggioInteractionPipeline.py -p P51955.pdb -sm catalytic_residues.csv -y ../config.yaml  # analysis on specified molecules
python ../ArpeggioInteractionPipeline.py -p P51955.pdb -sm catalytic_residues.csv -y ../config.yaml  # analysis on specified on pairs of molecules
python ../ArpeggioInteractionPipeline.py -p P51955.pdb -y ../config.yaml -i # analysis on all the entities in the strucutre keeping only the inter-chain contacts 
```

In the example folder there is a subfolder for each analysis supported by the script with the appropiate exmaple  
in order to reproduce the analysis in every folder run:

```
bash run.sh
```