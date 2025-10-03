import yaml
import argparse

# Load YAML
def load_yaml(file_path):
    with open(file_path, 'r') as file:
        return yaml.safe_load(file)

# Helper with defaults and optional warning
def get(data, key, default, warn=False):
    if warn and key not in data:
        print(f"[WARNING] Missing optional field: '{key}'")
    return data.get(key, default)

# Generate Markdown with original static text and YAML variables
def generate_markdown(data):

    md = f"""\
## Metadata

- **Gene Name**: {get(data, 'gene_name', 'N.A.')}
- **UniProt ID**: {get(data, 'uniprot_ac', 'N.A.')}
- **RefSeq ID**: {get(data, 'refseq_id', 'N.A.')}
- **Date of Aggregated Dataset**: {get(data, 'date_aggregated_dataset', 'N.A.')}
- **Date of EVE Run**: {get(data, 'date_eve_run', 'N.A.')}
- **MAVISp Mode**: {', '.join(get(data, 'mavisp_mode', [], warn=True))}
- **Structure Source (STABILITY & LONG_RANGE)**: {get(data, 'structure_source', 'N.A.')}
- **Trimmed Models (Default)**: {', '.join(get(data, 'trimmed_models_default', [], warn=True))}
- **Trimmed Models (AlloSigma2)**: {', '.join(get(data, 'trimmed_models_allosigma', [], warn=True))}
- **Structure Source (LOCAL_INTERACTION)**: {get(data, 'structure_local_interaction', 'N.A.')}
- **Date of AF Model Download**: {get(data, 'date_af_model_download', 'N.A.')}
- **Date of Cancermuts Run**: {get(data, 'date_cancermuts_run', 'N.A.')}
- **Date of Mentha2PDB Database**: {get(data, 'date_mentha2pdb_database', 'N.A.')}

## Curators:
"""
    for curator, info in get(data, 'curators', {}, warn=True).items():
        md += f"- {curator}, {'; '.join(info.get('affiliation', []))}\n"

    md += "\n## Reviewers:\n"
    for reviewer, info in get(data, 'reviewers', {}, warn=True).items():
        md += f"- {reviewer}, {'; '.join(info.get('affiliation', []))}\n"

    refs = get(data, 'reference_databases', {}, warn=True)
    md += f"""

## Reference databases
- [GeneCards]({refs.get('GeneCards', 'N.A.')})
- [Uniprot]({refs.get('UniProt', 'N.A.')})
- [NCBI]({refs.get('NCBI', 'N.A.')})
- [Human Protein Atlas]({refs.get('HumanProteinAtlas', 'N.A.')})
- [Alliance of Genome Resources]({refs.get('AllianceGenome', 'N.A.')})
- [MAVE Database]({refs.get('MAVEDatabase', 'N.A.')})
- [G2P]({refs.get('G2P', 'N.A.')})


## Overview on variants included in MAVISp
We included data on all the possibile variants in the trimmed region and annotations for the ones reported in COSMIC, CbioPortal or Clinvar (if any). The results can be seen in the figure below.
[Upload upset plot]


## Methods

### Source of variants
We retrieve the variants with [Cancermuts](https://github.com/ELELAB/cancermuts) to aggregate the data from ClinVar with the variants from COSMIC and cBioPortal.

"""

    # PROCHECK section
    procheck = get(data, 'procheck', {}, warn=True)
    rama = procheck.get('ramachandran', {})
    bad_contacts = procheck.get('bad_contacts', 'N.A.')
    if rama:
        md += f"""

### Procheck - model evaluation
The initial structure has been evaluated with ProCheck before using it for further analyses. Here the reported results. We aim to use structures that do not have residues in the disallowed regions of the Ramachandran plot and high 'bad contact' score.
Ramachandran plot: {rama.get('core', 'N.A.')} core {rama.get('allow', 'N.A.')} allow {rama.get('gener', 'N.A.')} gener {rama.get('disall', 'N.A.')} disall
Bad contacts: {bad_contacts}


### INFORMATION ON COFACTORS - ALPHAFILL
 [OPTIONAL]

### SIMPLE MODE - STABILITY module

We trimmed the model to remove regions with low pLDDT score or high PAE values. The folding free energy calculations have been carried out with FoldX5 using [MutateX](https://github.com/ELELAB/mutatex), cartddg2020 protocol and ref2015 energy function with [RosettaDDGPrediction](https://github.com/ELELAB/RosettaDDGPrediction) and [RaSP_workflow](https://github.com/ELELAB/RaSP_workflow). Templates and config files for MutateX, RosettaDDGPrediction and RaSP_workflow are provided here: https://osf.io/y3p2x/.

### EFOLDMINE

The EFOLDMINE module, predicts residues with propensity to belong to early folding regions using EfoldMine [10.1038/s41598-017-0836]. Trained on residue-level hydrogen/deuterium exchange nuclear magnetic resonance data extracted from the Start2Fold database [https://doi.org/10.1093/nar/gkv1185]. MAVISp is available in the simple mode of MAVISp and allow to identify mutation sites in the early folding regions, using a threshold of 0.169 and only regions with a minimum length of three early folding residues to avoid isolated peaks. It can be used to select variants to compare with the wild-type in the ensemble mode using biomolecular simulations to study the folding process.

"""

    # ENSEMBLE MODE and LONG_RANGE section
    md += f"""
### ENSEMBLE MODE - Molecular Dynamics and Stability

- Source of structure: {', '.join(get(data, 'ensemble_sources', ['N.A.'], warn=True))}
- Trimmed models for ensemble mode: {', '.join(get(data, 'ensemble_trimmed_regions', [], warn=False))}
- Force field: {get(data, 'simulation_force_field', ['N.A.'], warn=True)[0]}
- Simulation length: {get(data, 'simulation_length', ['N.A.'], warn=True)[0]}
- OSF trajectory data: [your_link_here]
- Number of frames for FoldX5 run: {', '.join(map(str, get(data, 'ensemble_size_foldx', ['N.A.'], warn=True)))}
- Number of clusters for Rosetta or RaSP runs: {', '.join(map(str, get(data, 'ensemble_size_rosetta', ['N.A.'], warn=True)))}

On the selected frames from the simulation for calculations of changes in free energy of folding with MutateX (foldX5 energy function) and the representative structures of the main clusters after clustering of the main-chain RMSD matrix with the GROMOS algorithm for calculations with Rosetta (cartddg2020 protocol and ref2015 energy function) or RaSP.

### ENSEMBLE MODE - CABS-flex

We used the same structure applied in the simple mode and carried out CABS-flex calculations with all-atom reconstruction to collect 20 representative models using the SnakeMake workflow for [CABS-flex data collection](https://github.com/ELELAB/MAVISp_CABSflex_pipeline) and the mavisp_template config files that can be downloaded from the MAVISp-extended data OSF repository: https://osf.io/ag7th.

### INTERACTOME

We use [PPI2PDB](https://github.com/ELELAB/PPI2PDB) to retrieve interactors reported with a Mentha score higher than 0.2. Mentha2PDB also queries experimentally validated dimeric complexes generated with AlphaFold2 from the HuRI and HuMAP databases by Burke et al. (https://doi.org/10.1038/s41594-022-00910-8). Additionally, we use STRING2PDB to extract interactors from the physical subnetwork of STRING, requiring a STRING score > 0.15 and evidence of interaction supported by either curated database annotation (database score > 0) or experimental data (experimental score > 0). Additionally, we retrieved the known complexes with missing residues from the Protein Data Bank using [PDMiner](https://github.com/ELELAB/PDBminer). STRINGS2PDB always use the latest version of the STRING database (i.e., version 12 which was released July 2023).

### ALPHAFOLD MULTIMER and ALPHAFOLD3

In the case the model of interaction or an experimental structure of the complex is not available, we produced models with a standalone version of Alphafold-Multimer or AlphaFold3 for the interactors reported in the Table 1 below in the Result Section. We verified the PAE scores and pDockQ scores, as well as the predicted interfaces with visual inspection with PyMOL and retained only the candidates with good quality prediction for further analyses. We also verified if any of the SLiMs reported by ELM and retrieved with Cancermuts  could link the protein with the possible interactors reported by Mentha2PDB.  In the case of AF3, we also used the TM scores in the assessment.

### SLiM-PROTEIN MODELING

We retrieved SLiMs in the proximity of the mutation sites identified by ELM using Cancermuts. If an interactor is known to recognize that class of motif this information can be used to design a model of their interaction. At the moment we support only modeling for LIR-LC3 proteins, BH3-BCL proteins, UIM-Ubiquitin(Ub) and BRCT domains.

#### LIR-LC3, BH3-BCL2-like and UIM-Ub MODELING

We used [phosphoiLIR](https://github.com/ELELAB/phospho-iLIR) and then apply [SLiMfast](https://github.com/ELELAB/SLiMfast) to filter out entries that are not satisfying structural criteria for LIRs. We retained only phospho-variants that have at least a NetPhos score of 0.4 for the corresponding phospho-site. The candidate LIRs can be used (in absence of a structure of the complex) to model the interaction with LC3 proteins following these attempts: i) AFmultimer, ii)AF-linker, iii) SLiM de novo pipeline. In the result we report if the LIRs are only predicted or updates on their experimental validation. We used a search based on regular expression for the BH3 motif according to the definition in: 10.1371/journal.pcbi.1007485. Then we use the SLiMfast protocol to evaluate if the motif can be retained for further analyses based on a structural filter. A similar protocol is applied to UIM motifs, following the definition in:  10.3389/fmolb.2021.676235 and 10.1016/s0968-0004(01)01835-7.

#### BRCT domains - TARGET MODELING

At first we used [DeepLoc2](https://services.healthtech.dtu.dk/services/DeepLoc-2.0/) to evaluate if the protein has a predicted nuclear localization (unless we know from experiments or literature). If it there is no evidence for it we do not continue in the analyses. If there are evidences we verify If a match is identified between a phospho-site on the target protein and an interactor with a BCRT domain, as well as if the phospho-site has an accessibility higher than 20% for its sidechain (as estimated by NACCESS). In absence of available structures, we use AFmulti to predict the complex of interaction and in presence of a reasonable PAE score, we compare the model to known complexes BRCT-protein structures before using it for mutational scans.

### PROTEIN-PROTEIN INTERACTIONS

We carried out binding free energy calculations upon mutation using FoldX5 with MutateX and RosettaDDGPrediction (flexddg protocol and talaris2014 energy function). Templates and config files for MutateX, RosettaDDGPrediction are provided here:  https://osf.io/y3p2x/.

### LOCAL INTERACTIONS WITH DNA
Source of structure:
Trimmed regions:
We used MutateX to collect FoldX-based predictions on the changes in binding free energy upon mutation.

[STATIC TEXT FOR ENSEMBLE MODE CONTINUES HERE...]

### LONG_RANGE MODULE

- long-range Allosigma2 cutoff applied:
- rSAS cutoff = {get(data, 'allosigma_parameters', {}).get('rSAS_cutoff', 'N.A.')}%
- allosteric free energy cutoff = +/- {get(data, 'allosigma_parameters', {}).get('free_energy_cutoff', 'N.A.')} kcal/mol
- distance cutoff = {get(data, 'allosigma_parameters', {}).get('distance_cutoff', 'N.A.')} Å
- distance method = {get(data, 'allosigma_parameters', {}).get('distance_method', 'N.A.')}

We use AlloSigma2 (10.1093/nar/gkaa338) in the simple mode of MAVISp to calculate the allosteric signaling map for the protein of interest. The data are processed using the [AlloSigma2 workflow](https://github.com/ELELAB/MAVISp_allosigma2_workflow). The analysis can be done on a list of functional response sites (i.e., catalytic site, cofactor binding sites...) or to find long-range effects in putative binding interfaces as calculated by Fpocket (10.1093/nar/gkq383). We retain only the response sites with a changes in an absolute change in allosteric free energy higher than a certain threshold in absolute value and in sites that are solvent exposed (only in the case of the pocket analysis). Additionally, we apply a distance cutoff distance to remove false positive response sites that are or could become (upon conformational changes) in direct contact with the target site.  This cutoff  has been benchmarked from the analysis of proteins with different folds that are included in the database. There is also the possibility to fine-tune the cutoff on the specific protein of interest to improve the signal.
In the ensemble mode, we apply a method based on path analysis from a atomic contact Protein Structure Network obtained with PyInteraph2 (10.1021/acs.jcim.3c00574). It allows to further validate the results from the coarse grain model AlloSigma2, retaining only the mutation sites for which a path of communication to their response site has been found in the PSN map.

[STATIC TEXT FOR LONG_RANGE CONTINUES HERE...]

### FUNCTIONAL SITES

- active site: {', '.join(get(data, 'functional_sites', {}).get('active_site', []))}
- substrate binding site: {', '.join(get(data, 'functional_sites', {}).get('substrate_binding', []))}
- cofactor binding site (cofactor: {get(data, 'functional_sites', {}).get('cofactor_binding', {}).get('cofactor', 'N.A.')}): {', '.join(get(data, 'functional_sites', {}).get('cofactor_binding', {}).get('residues', []))}

The module includes an assessment of the effect of the mutations occurring at sites important for the catalytic activity, substrate recognition or cofactor binding of an enzyme or a protein based on structural proximity. Damaging variants are the ones that occurrs at these sites or in the second sphere of interactions for the site, as evaluated by Arpeggio ( 10.1016/j.jmb.2016.12.004).

[STATIC TEXT FOR FUNCTIONAL SITES CONTINUES HERE...]

## RESULTS

### INTERACTOME

We identified XXX interactors with PPI2PDB. Below we summarize the status of the currently processed interactions. According to the MAVISp defaults, we applied a priority list based on the interactors with Mentha scores > 0.4 or reported by STRING and PDBMiner.  For each interactor,  in case of not available PDB structure or precalculated model, we first gather the available publications reported by Mentha or through  literature searches with Google Scholar to collect evidence of a physical interaction and on binding interfaces to guide the modeling. The table below is currently not complete since there are many potential interactors for this protein with high scores.


- Table 1. Summary of interactors from PPI2PDB processed within MAVISp LOCAL_INTERACTION module

Last update 13/09/2025


"""

    # SLiMs section
    md += """
### SLiMs ANALYSIS
"""
    for slim_type, slim_info in get(data, 'slims', {}, warn=True).items():
        md += f"\n{slim_type} motifs:\n"
        if slim_info.get('present'):
            for motif in slim_info.get('motifs', []):
                md += f"- positions: {motif.get('positions', 'N.A.')}\n"
                if 'accessibility' in motif:
                    md += f"  accessibility: {motif['accessibility']}\n"
                if 'nuclear_localization_probability' in motif:
                    md += f"  nuclear localization probability: {motif['nuclear_localization_probability']}\n"
                if 'phosphosites' in motif:
                    for site in motif['phosphosites']:
                        md += f"  phosphosite position: {site.get('position', 'N.A.')}, kinases: {', '.join(site.get('kinases', []))}\n"
                if 'interactors' in motif:
                    md += f"  interactors: {', '.join(motif['interactors'])}\n"
                if 'notes' in motif:
                    md += f"  notes: {motif['notes']}\n"
        else:
            md += "- None identified\n"

    md += """


### MAVISp CONSEQUENCES

The predicted consequences  for each mutation can be consulted from the MAVISp webserver with a dotplot representation:  https://services.healthtech.dtu.dk/services/MAVISp-1.0/
Additional information are provided in the following

#### VARIANTS in ClinVar and VUS

This section provides information on variants retrieved from ClinVar, prediction on pathogenicity and their predicted consequences according to MAVISp.  For the prediction of damaging effects we here report primarily the results from AlphaMissense.  The full dataset can be navigated by the MAVISp database.
[file from piechart here if relevant]

[insert dotplot and lolliplot after the statements below if relevant (see ACVR1B example)]

The results for the pathogenic and benign variants of ClinVar are reported below:

The results for the VUS or Variants with Conflicting Evidence of ClinVar are reported below:    

#### OTHER PREDICTED PATHOGENIC VARIANTS WITH ASSOCIATED MECHANISTIC INDICATOR

We then evaluated the variants in COSMIC for which AlphaMissense provided a pathogenic prediction, along with there is a predicted loss-of-function or gain-of-function effects. For these we identified as possible mechanism with MAVISp (combining simple and ensemble mode, in case both information are available).
[a similar analysis can be done for other interesting subsets of mutations, see examples in ACVR1B]
[IMPORT LOLLIPLOT HERE]

We carried out a similar analysis for Cbioportal and the results of the predicted pathogenic variants with mechanistic indicators can be seen in the lolliplot below
[IMPORT LOLLIPLOT HERE]

## ADDITIONAL ANALYSES

[to be used if there are interesting variants or data you want to report, see example of ACVR1B]

variant_name:
Investigation:

## ADDITIONAL NOTES FROM CURATORS

include additional notes here that are not covered by the sections above

## PUBLICATIONS INCLUDING THESE DATA

#to customize below with the proper refs - if the proteins are in Table S1 for the first paper you can use for now the biorxiv ref:

- Arnaudi M., Beltrame L., Degn K., Utichi M. et al. biorxiv, https://doi.org/10.1101/2022.10.22.513328]\n"

"""
    return md

# Entry point
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate MAVISp report from a single YAML file')
    parser.add_argument('--yaml', required=True, help='Path to combined YAML file')
    parser.add_argument('--output', default='mavisp_report.md', help='Output Markdown file name')

    args = parser.parse_args()

    data = load_yaml(args.yaml)
    markdown_report = generate_markdown(data)

    with open(args.output, 'w') as md_file:
        md_file.write(markdown_report)
