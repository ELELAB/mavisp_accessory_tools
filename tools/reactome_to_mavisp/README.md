# Reactome UniProt Reaction Workflow

This pipeline runs an automated Reactome analysis for one or more UniProt accessions. For each target protein, the workflow uses a local Reactome release to retrieve NORMAL and DISEASE events, reconstruct pathway context and ordering from BioPAX, extract reaction/protein/complex/disease-variant annotations, and write cleaned CSV outputs in MAVISp-supported format for downstream analysis.

The code is organized as a small Python package. The command-line entry point is reactome_to_mavisp.py, while the core logic is split across the reactome_pipeline/ modules.

---

## Requirements

### Python


Python >= 3.8


### Python packages

Required packages:


pybiopax

pandas

numpy

networkx

requests

matplotlib


The workflow also uses standard-library modules such as argparse, os, shutil, re, json, zipfile, urllib.request, datetime, collections, pathlib, and typing.

Example installation:


pip install pybiopax pandas numpy networkx requests matplotlib


---

## Description

### Main files

| File | Role |

|---|---|

| reactome_to_mavisp.py | Main command-line entry point. Parses arguments, validates or refreshes the local Reactome release files, initializes one shared local Reactome database for all targets, writes metadata.json, runs the workflow, and records accessions that do not produce valid output in entries_not_in_reactome.csv. |

| reactome_pipeline/workflow.py | Contains the ReactomeScript class and orchestrates the full local Reactome workflow for one UniProt accession, including NORMAL/DISEASE contexts, pathway de-duplication, local pathway ordering, and final result.csv generation. |

| reactome_pipeline/local_reactome.py | Loads UniProt2Reactome_PE_Reactions.txt and Homo_sapiens.owl, indexes BioPAX events and pathways, reconstructs pathway ancestry, and extracts local PathwayStep ordering relationships. |

| reactome_pipeline/event_context.py | Builds a generic local EventContext for one Reactome event, including event roles, physical entities, complexes, sequence features, controls, pathway chains, and event-specific BioPAX annotations. |

| reactome_pipeline/context_factory.py | Builds independent NORMAL and DISEASE EventContexts for a UniProt target. NORMAL reactions come from UniProt2Reactome_PE_Reactions.txt; DISEASE reactions come from disease_variant_ewas_mapping.tsv. |

| reactome_pipeline/disease_variants.py | Loads and indexes disease_variant_ewas_mapping.tsv by UniProt accession and disease reaction, preserving variant, disease, disease identifier, functional status, literature, and normal/disease reaction/pathway metadata. |

| reactome_pipeline/legacy_adapter.py | Converts the generic EventContext representation into the nested structure expected by DataProcessingFunctions, without querying Reactome. |

| reactome_pipeline/data_processing.py | Flattens nested reaction/pathway/protein annotations into a cleaned pandas.DataFrame, matches disease metadata to the exact mutant physical entity, expands variant × disease associations, derives the compact mutation field, reorders columns, removes duplicates, and propagates complex information. |

| reactome_pipeline/graph_utils.py | Builds and processes directed pathway graphs from local BioPAX PathwayStep relationships and identifies start/end nodes and cycle membership. |

| reactome_pipeline/uniprot_utils.py | Contains UniProt helper functions for gene-to-accession conversion, accession-to-protein-name retrieval, and parsing protein names. |

| reactome_post_process.py | Optional post-processing script that merges individual result.csv files into summary CSV tables and generates the highest_pathways.pdf and disease_targets.pdf plots. |

| reactome_pipeline/__init__.py | Marks reactome_pipeline/ as a Python package. It can remain empty. |

The scripts are organized as follow:


project/

├── example/
│   ├── reactome_outputs/
│   │   └── Q8N726/
│   ├── readme.txt
│   ├── run.sh
│   └── uniprot_list.txt
├── reactome_data/
│   ├── disease_variant_ewas_mapping.tsv
│   ├── Homo_sapiens.owl
│   └── UniProt2Reactome_PE_Reactions.txt
├── reactome_outputs/
│   ├── metadata.json
│   ├── P04637/
│   │   ├── pathways_order/
│   │   └── result.csv
│   ├── Q8N726/
│   │   ├── pathways_order/
│   │   └── result.csv
│   └── summary/
│       ├── disease_single_sequence_site.csv
│       ├── disease_targets.pdf
│       ├── highest_pathways.pdf
│       ├── merged_highest_pathways.csv
│       └── merged_reaction.csv
├── reactome_pipeline/
│   ├── context_factory.py
│   ├── data_processing.py
│   ├── disease_variants.py
│   ├── event_context.py
│   ├── graph_utils.py
│   ├── __init__.py
│   ├── legacy_adapter.py
│   ├── local_reactome.py
│   ├── uniprot_utils.py
│   └── workflow.py
├── reactome_post_process.py
├── reactome_to_mavisp.py
├── README.md
└── uniprot_list.txt


Python cache files such as __pycache__/ or *.pyc are generated automatically and should not be edited or tracked manually.

For each UniProt accession, the workflow performs the following steps.

### 1. Initial Reactome mapping check

Before running the per-target analysis, reactome_to_mavisp.py validates the three required local Reactome files:


reactome_data/UniProt2Reactome_PE_Reactions.txt
reactome_data/Homo_sapiens.owl
reactome_data/disease_variant_ewas_mapping.tsv


UniProt2Reactome_PE_Reactions.txt is filtered to human entries and indexed as UniProt accession → Reactome reaction IDs. disease_variant_ewas_mapping.tsv is independently indexed as UniProt accession → disease reaction IDs and reaction ID → disease-variant records.

The local BioPAX model is loaded once and shared across all input targets.

A target can therefore produce output through NORMAL events, DISEASE events, or both.

If neither branch produces any valid local Reactome event/output, the accession is written to:


entries_not_in_reactome.csv


with the status:


no_valid_reactome_output


This avoids repeated web queries and prevents the large Homo_sapiens.owl model from being reloaded for every UniProt accession.

---

### 2. Retrieve Reactome pathways

NORMAL reaction discovery starts from UniProt2Reactome_PE_Reactions.txt, while DISEASE reaction discovery starts independently from disease_variant_ewas_mapping.tsv.

For every event, pathway context is reconstructed from the local Homo_sapiens.owl model. LocalReactomeDatabase indexes Reactome pathways and parent relationships and EventContextBuilder derives the pathway chains associated with each event.

This allows the final output to report not only the specific low-level pathway containing a reaction, but also the broader pathway hierarchy in which that reaction occurs.

NORMAL and DISEASE contexts remain independent even when they refer to related normal/disease biological processes.

---

### 3. Optional pathway ordering

Unless --skip_pathway_order is used, the workflow attempts to infer the order of reactions within each lowest-level pathway.

Ordering is reconstructed entirely from the local BioPAX model. For each pathway, LocalReactomeDatabase.get_pathway_reaction_order() uses explicit BioPAX PathwayStep relationships, including:


reaction.step_process_of
PathwayStep.next_step
PathwayStep.next_step_of
PathwayStep.step_process


These relationships are used to build a directed graph where:

nodes represent Reactome event stable IDs;

edges represent explicit next/previous PathwayStep relationships;

isolated reactions are retained as graph nodes;

start/end nodes are defined from graph in-degree/out-degree;

cycle membership is recorded;

non-redundant start-to-end paths are written to ordered_paths.csv when at least one linear path exists.

When --skip_pathway_order is used, the workflow skips graph construction and still produces the main result.csv output with ordered=False.

---

### 4. Resolve target information

Before writing the final result, the workflow collects target-level information for the input UniProt accession. This is done in:


ReactomeScript.resolve_target_information()


Reactome event discovery remains fully local. The only target-level external lookup used by the current workflow is the existing UniProt helper used to obtain a readable protein name.

The final target_name is selected using the following priority:

protein name retrieved from UniProt;

the original UniProt accession as fallback.

The target information dictionary contains:

| Field | Meaning |

|---|---|

| target_name | Final target name written in result.csv. Preferentially retrieved from UniProt; otherwise set to the input UniProt accession. |

---

### 5. Collect candidate reactions

The current workflow discovers target events directly from the two local Reactome mappings.

For NORMAL contexts:


UniProt accession
    -> UniProt2Reactome_PE_Reactions.txt
    -> NORMAL reaction IDs


For DISEASE contexts:


UniProt accession
    -> disease_variant_ewas_mapping.tsv
    -> DISEASE reaction IDs


Disease reactions are removed from the NORMAL reaction set so that the two event branches remain independent.

The disease index also stores the disease-variant records associated with each disease reaction, including exact mutant physical-entity stable IDs.

---

### 6. Identify target reactions

The script builds EventContexts independently for the NORMAL and DISEASE reaction IDs associated with the requested UniProt accession.

Each EventContext is created from the same local BioPAX model and contains:

the Reactome event;

event roles and participants;

a registry of physical entities and complexes;

sequence features and cellular locations;

controls and catalysis;

pathway chains;

event-specific BioPAX metadata.

DISEASE contexts additionally receive only the disease-variant records belonging to the current target UniProt accession.

During tabular processing, mutation metadata are matched to a protein row by exact equality between the row physical-entity stable ID and variant_entity_stable_id. This prevents disease/mutation metadata from being copied onto unrelated partners in the same disease reaction.

A real non-mutant partner of a DISEASE reaction is therefore retained as a DISEASE-context protein row, but its variant-specific fields remain empty.

---

### 7. Extract BioPAX annotations

For each local Reactome EventContext, the workflow extracts detailed reaction-level and physical-entity annotations from Homo_sapiens.owl.

The extracted information includes:

protein display name;

UniProt accession;

cellular location;

sequence intervals;

sequence sites;

modification type;

direct complex membership and complex stable IDs;

explicit BioPAX component stoichiometry;

EntitySet/protein-family membership;

pathway hierarchy;

reaction name and Reactome stable ID;

biochemical left/right participants;

conversion direction;

reaction EC numbers;

explicit catalytic EC numbers, when present;

regulatory controllers;

NORMAL/DISEASE event context;

disease and disease identifiers;

disease/normal pathway and reaction mappings;

exact disease-variant stable IDs;

compact mutation labels;

variant modification metadata;

functional status;

PubMed references.

Complex semantics are kept explicit: complex_of is derived from direct BioPAX component relationships, while complex_entity_set represents outer EntitySet-like membership through member_physical_entity. Missing stoichiometry is not inferred.

This step produces nested annotation dictionaries that are converted by LegacyEventAdapter into the structure used by DataProcessingFunctions.

---

### 8. Build and write output tables

The nested annotation dictionaries are flattened into a pandas.DataFrame.

The final table is then:

cleaned;

deduplicated;

column-ordered;

expanded to one row per variant × disease association when one variant has multiple disease associations;

matched so variant-specific metadata are attached only to the exact mutant physical entity;

filtered to remove protein-family rows;

filtered to remove redundant pathway representations while keeping NORMAL and DISEASE contexts independent;

optionally reordered using locally generated pathway-ordering files;

written to result.csv.

Protein-family rows are removed so that the final output focuses on individual protein entries rather than broad family-level Reactome entities.

The final result.csv therefore contains both NORMAL and DISEASE event contexts. Disease reactions can contain mutant rows with disease/variant metadata and genuine non-mutant partner rows with those variant-specific columns empty.

### Optional post-processing

After running the main workflow for multiple UniProt accessions, the optional script reactome_post_process.py merges individual result.csv files, harmonizes target columns when needed, removes duplicate summary rows, writes three post-processed CSV files, and generates two PDF plots.

merged_reaction.csv contains a compact reaction-level summary across all analysed UniProt accessions. It keeps the target accession, target name, highest pathway, disease annotation, dynamic intermediate pathway hierarchy, lowest-level pathway, reaction name, reaction ID, left/right participants, and reaction direction.

merged_highest_pathways.csv contains a simplified pathway-level summary. It reports the highest-level Reactome pathways associated with each target protein, together with the corresponding Reactome pathway ID, UniProt accession, and target name.

disease_single_sequence_site.csv contains a filtered subset of the concatenated results. It keeps only rows where highest_pathway == Disease and SequenceSite contains one numeric residue position.

highest_pathways.pdf shows a target × highest-pathway matrix. Duplicate target/pathway combinations are removed before plotting.

disease_targets.pdf shows target × disease associations. functional_status values are summarized as loss of function, gain of function, mixed, or unknown. Rows without a disease_name, including non-mutant partners of DISEASE reactions, are not plotted.

---

## Input

Input for reactome_to_mavisp.py script:

| Argument | Description |

|---|---|

| -u, --uniprot_ac | UniProt accession to analyze. Default: Q8N726. Ignored when --uniprot_file is supplied. |

| -uf, --uniprot_file | Text file containing one UniProt accession per line. Blank lines and lines starting with # are ignored. Duplicate accessions are removed while preserving input order. |

| -o, --output_dir | Main output directory. Default: reactome_outputs. A subfolder is created for each accession. |

| -s, --skip_pathway_order | Skip local pathway-order inference. result.csv is still produced and reactions are marked ordered=False. |

| -r, --refresh_reactome_data | Download the current Reactome release files before analysis, replacing the local reaction mapping, disease-variant mapping, and Homo_sapiens.owl. The current Reactome release number is queried and stored in metadata. |

| --reaction_map_file | Path to UniProt2Reactome_PE_Reactions.txt. Default: reactome_data/UniProt2Reactome_PE_Reactions.txt. |

| --biopax_file | Path to Homo_sapiens.owl. Default: reactome_data/Homo_sapiens.owl. |

| --disease_variant_file | Path to disease_variant_ewas_mapping.tsv. Default: reactome_data/disease_variant_ewas_mapping.tsv. |

| --reactome_release | Reactome release label stored in metadata.json when existing local files are used. |

| --reactome_download_date | Download date stored in metadata.json when supplied. If unknown, the BioPAX file modification date is used as fallback metadata. |

Here an example of file with a list of uniprot ac


P04637

Q8N726

Q9Y2X3


Input for reactome_post_process.py script:

Arguments:

| Argument | Description |

|---|---|

| -i, --input_dir | Main Reactome output directory containing UniProt-specific folders. |

| -o, --output_dir | Output directory for merged tables and plots. Default: same as --input_dir. |

| --result_filename | Name of the result file inside each UniProt folder. Default: result.csv. |

---

## Output

By default, outputs are written under:


reactome_outputs/


For accessions such as P04637 and Q8N726, together with optional post-processing, the output structure is:


reactome_outputs/

├── metadata.json
├── entries_not_in_reactome.csv        # Only created when at least one accession produces no valid output
├── P04637/
│   ├── result.csv
│   └── pathways_order/                # Only populated when pathway ordering is enabled
│       └── <Reactome_pathway_ID>/
│           ├── graph_edges.csv
│           ├── graph_nodes.csv
│           └── ordered_paths.csv      # Only written when at least one linear path exists
├── Q8N726/
│   ├── result.csv
│   └── pathways_order/
└── summary/
    ├── merged_reaction.csv
    ├── merged_highest_pathways.csv
    ├── disease_single_sequence_site.csv
    ├── highest_pathways.pdf
    └── disease_targets.pdf


metadata.json records the Reactome release metadata and the local filenames used for BioPAX, UniProt/reaction mapping, and disease-variant mapping.

### result.csv

result.csv is the main output of the workflow. Each row corresponds to an individual protein entry in one Reactome event context involving the target. NORMAL and DISEASE contexts are explicitly separated by event_context_type.

Common columns include:

| Column | Description |

|---|---|

| target_uniprot_ac | Input UniProt accession used as the analysis target. |

| target_name | Final target protein name. Preferentially retrieved from UniProt; otherwise the target UniProt accession. |

| highest_pathway | Highest-level Reactome pathway in the reconstructed hierarchy. |

| highest_pathway_id | Reactome stable ID of the highest-level pathway. |

| pathway_1, pathway_2, ... | Dynamic intermediate pathway hierarchy levels. The number of levels depends on the pathway. |

| pathway_1_id, pathway_2_id, ... | Reactome stable IDs corresponding to the intermediate pathway levels. |

| lowest_pathway | Lowest-level pathway associated with the event. |

| lowest_pathway_id | Reactome stable ID of the lowest-level pathway. |

| reaction_name | Reactome event/reaction display name. |

| reaction_id | Reactome stable event/reaction ID. |

| reaction_Left | Left-side participants for conversion-like BioPAX reactions, when applicable. |

| reaction_Right | Right-side participants for conversion-like BioPAX reactions, when applicable. |

| reaction_Conversion_Direction | BioPAX conversion direction, when available. |

| reaction_EC_Number | EC number explicitly associated with the BioPAX biochemical reaction. |

| catalytic_EC_Number | EC number explicitly associated with a Catalysis record, when present. It is not copied from reaction_EC_Number. |

| Controller_of_reaction_ACTIVATION | Controllers explicitly annotated as activation of the reaction. |

| Controller_of_reaction_INHIBITION | Controllers explicitly annotated as inhibition of the reaction. |

| protein | Display name of the protein/PhysicalEntity represented by the row. |

| uniprot_ac | UniProt accession(s) of the protein represented by the row. This can differ from target_uniprot_ac for reaction partners. |

| cellular_location | BioPAX cellular location of the protein entry. |

| SequenceInterval | Sequence interval annotation, when available. |

| SequenceSite | Sequence-site positions derived from BioPAX features. This is not the compact mutation label. |

| Modification_type | BioPAX modification annotation, when available. |

| complex_of | Direct BioPAX Complex(es) containing the protein through component relationships. |

| complex_of_stid | Reactome stable ID(s) corresponding to complex_of. |

| complex_entity_set | Outer EntitySet-like parent(s) containing a direct complex through member_physical_entity. |

| stoichiometry | Explicit BioPAX stoichiometric coefficient for the corresponding direct complex membership. Missing coefficients are not inferred. |

| is_a_protein_family | Boolean flag indicating whether the BioPAX protein entry represents an EntitySet/protein family. Such rows are removed from the final result.csv. |

| member_physical_entity_of | Parent protein EntitySet/family of the row protein, when available. |

| event_context_type | Event branch: NORMAL or DISEASE. |

| disease_name | Disease associated with the exact disease variant. Empty for NORMAL rows and for non-mutant partners in a DISEASE reaction. |

| disease_cross_reference | Disease cross-reference(s) from disease_variant_ewas_mapping.tsv, for example MONDO terms when present. |

| disease_identifier | Disease identifier paired positionally with disease_name, for example a DOID. |

| disease_pathway_id | Disease pathway stable ID(s) associated with the disease-variant record. |

| disease_pathway_name | Disease pathway name(s) associated with the disease-variant record. |

| disease_reaction_id | Disease reaction stable ID(s) from the disease-variant mapping. |

| disease_reaction_name | Disease reaction name(s) from the disease-variant mapping. |

| normal_pathway_id | Corresponding normal pathway stable ID(s) from the disease-variant mapping. |

| normal_pathway_name | Corresponding normal pathway name(s). |

| normal_reaction_id | Corresponding normal reaction stable ID(s). |

| normal_reaction_name | Corresponding normal reaction name(s). |

| mutation | Compact mutation label extracted from the Reactome disease variant display name, for example R81G, R98L;R99S, or V22Pfs*46. Multiple mutation labels are retained when present. |

| variant | Full Reactome disease-variant display name, including contextual text such as the protein name and cellular compartment. |

| variant_entity_stable_id | Reactome stable ID of the exact mutant physical entity. This field is used to match disease metadata to the correct protein row. |

| variant_uniprot_ac | UniProt accession associated with the disease-variant record. |

| variant_modification_class | Reactome modification class of the variant, such as ReplacedResidue or FragmentReplacedModification. |

| variant_modification_description | Detailed Reactome description of the sequence change. |

| functional_status | Functional-status annotation from the disease-variant mapping, for example loss_of_function, gain_of_function, or combined status strings. |

| variant_literature_pubmed | PubMed identifier(s) supporting the disease-variant annotation. |

| ordered | Boolean flag indicating whether the reaction belongs to the selected locally reconstructed ordered path. False when no applicable ordering is available or pathway ordering is skipped. |

The exact number of pathway_N / pathway_N_id columns is dynamic and depends on the deepest pathway hierarchy present in the analysed results.

For disease variants associated with more than one disease, the workflow writes one row per variant × disease pair, preserving positional pairing between disease_name and disease_identifier.

### skipped_reactions.csv

The current fully local workflow does not generate skipped_reactions.csv.

The previous web-service-based implementation used this file to record reactions that failed online metadata/BioPAX retrieval. In the current implementation, reaction/event discovery and BioPAX annotation are obtained from the validated local Reactome release, so those retry-based skip categories are no longer part of the workflow.

If a target produces no valid NORMAL or DISEASE output after local processing, the target is instead recorded at the global level in entries_not_in_reactome.csv.

### entries_not_in_reactome.csv

This file is written at the global output-directory level when at least one input accession does not produce a valid output.

Possible statuses:

| Status | Meaning |

|---|---|

| no_valid_reactome_output | No NORMAL or DISEASE Reactome events produced a valid final output from the local release for the UniProt accession. |

Common columns:

| Column | Description |

|---|---|

| uniprot_ac | UniProt accession. |

| status | Failure/status category. |

| reason | Explanation of why no final output was produced. |

### pathways_order/

This directory is created when pathway ordering is enabled.

For each lowest-level pathway, the workflow can write:

| File | Description |

|---|---|

| graph_edges.csv | Directed edges between Reactome event stable IDs reconstructed from explicit BioPAX PathwayStep ordering. |

| graph_nodes.csv | Reaction/event nodes with is_start, is_end, in_degree, out_degree, and in_cycle information. |

| ordered_paths.csv | Non-redundant start-to-end reaction paths inferred from the directed graph. This file is written only when at least one linear path exists. |

When --skip_pathway_order is used, pathway ordering is not generated and reactions in result.csv are marked as unordered.

### post process analysis

The post-processing script writes:

| File | Description |

|---|---|

| merged_reaction.csv | Deduplicated reaction-level summary across all analyzed UniProt accessions, including target, pathway hierarchy, disease name, reaction identifiers, reaction participants, and direction. |

| merged_highest_pathways.csv | Deduplicated table of highest-level pathways per target. |

| disease_single_sequence_site.csv | Subset of Disease highest-pathway rows where SequenceSite is a single numeric residue position. |

| highest_pathways.pdf | Target × highest-pathway overview plot. |

| disease_targets.pdf | Target × disease plot. Functional status is summarized as LoF, GoF, mixed, or unknown. |

---

## Run

### Run one UniProt accession


python reactome_to_mavisp.py -u P04637


### Run one UniProt accession and skip pathway ordering


python reactome_to_mavisp.py -u P04637 -s


Skipping pathway ordering is faster because it avoids building local pathway-level reaction graphs.

### Run with input list


python reactome_to_mavisp.py -uf uniprot_list.txt -o reactome_outputs


### Run post process analysis


python reactome_post_process.py -i reactome_outputs -o reactome_outputs/summary


---

##  example


# Refresh the three local Reactome files to the current release
python reactome_to_mavisp.py -r -uf uniprot_list.txt -o reactome_outputs

# Single protein, faster run without pathway ordering
python reactome_to_mavisp.py -u P04637 -s

# Multiple proteins
python reactome_to_mavisp.py -uf uniprot_list.txt -o reactome_outputs -s

# Merge results and generate the two plots
python reactome_post_process.py -i reactome_outputs -o reactome_outputs/summary