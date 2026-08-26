# Reactome UniProt Reaction Workflow

This pipeline runs an automated Reactome analysis for one or more UniProt accessions. For each target protein, the workflow retrieves human Reactome pathways, identifies reactions that contain the target, parses BioPAX-level reaction/protein/complex annotations, and writes cleaned CSV outputs in MAVISp supported format for downstream analysis.

The code is organized as a small Python package. The command-line entry point is `reactome_to_mavisp.py`, while the core logic is split across the `reactome_pipeline/` modules.

---

## Requirements

### Python

```text
Python >= 3.8
```

### Python packages

Required packages:

```text
reactome2py
pybiopax
pandas
numpy
networkx
requests
```

The workflow also uses standard-library modules such as `argparse`, `os`, `shutil`, `re`, `time`, `copy`, `contextlib`, `io`, `collections`, `pathlib`, and `typing`.

Example installation:

```bash
pip install reactome2py pybiopax pandas numpy networkx requests
```

---

## Description

### Main files

| File | Role |
|---|---|
| `reactome_to_mavisp.py` | Main command-line entry point. Parses arguments, checks whether a UniProt accession is mapped in Reactome, runs the workflow, and writes `entries_not_in_reactome.csv` when needed. |
| `reactome_pipeline/workflow.py` | Contains the `ReactomeScript` class and orchestrates the full Reactome workflow for one UniProt accession. |
| `reactome_pipeline/reactome_analysis.py` | Contains Reactome/BioPAX helper functions, safe retry logic for Reactome calls, BioPAX parsing, disease-link extraction, and target-reaction checks. |
| `reactome_pipeline/data_processing.py` | Flattens nested reaction/pathway/protein annotations into a cleaned `pandas.DataFrame`, reorders columns, removes duplicates, and propagates complex information. |
| `reactome_pipeline/graph_utils.py` | Builds and processes directed reaction graphs used for optional pathway-ordering analysis. |
| `reactome_pipeline/uniprot_utils.py` | Contains UniProt helper functions for gene-to-accession conversion, accession-to-protein-name retrieval, and parsing protein names. |
| `reactome_post_process.py` | Optional post-processing script that merges individual `result.csv` files into summary CSV tables. |
| `reactome_pipeline/__init__.py` | Marks `reactome_pipeline/` as a Python package. It can remain empty. |

The scripts are organized as follow:

```text
project/
├── reactome_to_mavisp.py
├── reactome_post_process.py
├── reactome_pipeline/
│   ├── __init__.py
│   ├── workflow.py
│   ├── reactome_analysis.py
│   ├── data_processing.py
│   ├── graph_utils.py
│   └── uniprot_utils.py
└── README.md
```

Python cache files such as `__pycache__/` or `*.pyc` are generated automatically and should not be edited or tracked manually.

For each UniProt accession, the workflow performs the following steps.

### 1. Initial Reactome mapping check

Before running the full workflow, `reac_classified.py` first checks whether the input UniProt accession is mapped to at least one human Reactome pathway.

This check is performed using:

```python
rc.content.mapping(
    id=uniprot_ac,
    resource="UniProt",
    species="9606",
    by="pathways"
)
```

This initial step is used only to decide whether the accession should be processed further.

If Reactome returns no mapped human pathways, the accession is not analysed and is written to:

```text
entries_not_in_reactome.csv
```

with the status:

```text
not_found_in_reactome
```

This prevents the workflow from spending time on accessions for which Reactome has no pathway-level annotation.

---

### 2. Retrieve Reactome pathways

For accessions that pass the initial mapping check, the workflow retrieves all human Reactome pathways associated with the UniProt accession using:

```python
rc.content.mapping(
    id=uniprot_ac,
    resource="UniProt",
    species="9606",
    by="pathways"
)
```

Each retrieved pathway is later expanded into its full Reactome hierarchy using ancestor information.

This allows the final output to report not only the specific low-level pathway containing a reaction, but also the broader pathway context in which that reaction occurs.

---

### 3. Optional pathway ordering

Unless `--skip_pathway_order` is used, the workflow attempts to infer the order of reactions within each retrieved pathway.

For each pathway, the workflow downloads the corresponding BioPAX model and extracts next/previous reaction relationships. These relationships are used to build a directed graph where:

- nodes represent reactions;
- edges represent next/previous relationships between reactions;
- inferred reaction paths are written to `ordered_paths.csv`.

The pathway-ordering step is optional because it can be slow for proteins associated with many Reactome pathways. This is due to the need for multiple BioPAX downloads and graph operations.

When `--skip_pathway_order` is used, the workflow skips this step and still produces the main `result.csv` output.

---

### 4. Resolve target information

Before filtering reactions, the workflow collects target-level information for the input UniProt accession. This is done in:

```python
ReactomeScript.resolve_target_information()
```

The purpose of this step is to prepare the identifiers that will later be used to decide whether a candidate Reactome reaction actually contains the input protein.

The workflow collects:

- the readable name of the target protein, which will be written in `result.csv`;
- optional Reactome stable identifiers for the target protein;
- alternative Reactome forms of the same protein;
- Reactome complexes associated with the UniProt accession;
- the identifiers needed for target-reaction detection.

The final `target_name` is selected using the following priority:

1. protein name retrieved from UniProt;
2. Reactome display/name retrieved from `search_fireworks`, if that optional call works;
3. the original UniProt accession as fallback.

The `search_fireworks` call is used only as an optional extra source of Reactome-specific protein identifiers. These identifiers can help detect reactions containing specific Reactome protein forms of the target.

However, `search_fireworks` is not required for the workflow to continue. If Reactome returns a server error or no valid entries, the workflow continues using:

- Reactome complexes associated with the UniProt accession;
- direct UniProt accession matching inside BioPAX reaction models.

The target information dictionary contains:

| Field | Meaning |
|---|---|
| `target_name` | Final target name written in `result.csv`. Preferentially retrieved from UniProt, otherwise from Reactome, otherwise set to the input UniProt accession. |
| `target_protein_stId` | Reactome stable IDs corresponding to the target protein, when available from `search_fireworks`. Used as an optional way to recognise reactions containing the target. |
| `target_protein_stId_other` | Alternative Reactome forms of the target protein, retrieved from Reactome when target stable IDs are available. |
| `all_target_protein_stId` | Combined list of direct and alternative Reactome target IDs. Used later during reaction filtering. |
| `target_protein_complexes_names` | Reactome complexes associated with the UniProt accession. Used later to identify reactions where the target appears as part of a complex. |

---

### 5. Collect candidate reactions

After retrieving the Reactome pathways, the workflow collects the reactions that may potentially involve the input protein.

For each lowest-level pathway, the workflow retrieves the contained events and keeps only true Reactome reactions.

Pathway containers are excluded to avoid analysing higher-level pathway objects as if they were individual biochemical reactions.

At this stage, the reactions are still considered **candidate reactions**. This means that they belong to Reactome pathways associated with the input UniProt accession, but they have not yet been confirmed to directly contain the target protein.

---

### 6. Identify target reactions

The script takes:

- the candidate reactions collected from Reactome pathways;
- the target information prepared by `resolve_target_information()`.

For each candidate reaction, the workflow checks whether the reaction actually contains the input protein.

A reaction is considered a **target reaction** if at least one of the following checks is true:

- the reaction contains a Reactome complex associated with the target UniProt accession;
- the reaction contains a Reactome protein stable ID corresponding to the target protein or one of its alternative forms;
- the BioPAX model of the reaction directly contains the input UniProt accession.

The third check is the most robust one because it searches BioPAX protein and entity-reference cross-references directly for the input UniProt accession.

Only reactions that pass this filtering step are kept for detailed annotation extraction.

---

### 7. Extract BioPAX annotations

For each reaction confirmed to contain the target protein, the workflow downloads or reuses the corresponding BioPAX model and extracts detailed reaction-level annotations.

The extracted information includes:

- protein display name;
- UniProt accession;
- cellular location;
- sequence intervals;
- sequence sites;
- modification type;
- complex membership;
- stoichiometry;
- parent protein family or physical entity;
- pathway hierarchy;
- reaction name and Reactome stable ID;
- biochemical left/right participants;
- conversion direction;
- regulatory controllers;
- disease links.

This step produces nested annotation dictionaries describing the target-containing reactions.

---

### 8. Build and write output tables

The nested annotation dictionaries are flattened into a `pandas.DataFrame`.

The final table is then:

- cleaned;
- deduplicated;
- column-ordered;
- filtered to remove protein-family rows;
- optionally reordered using pathway-ordering files;
- written to `result.csv`.

Protein-family rows are removed so that the final output focuses on individual protein entries rather than broad family-level Reactome entities.

The final `result.csv` therefore contains only Reactome reactions that passed the target-reaction filtering step and for which BioPAX-level annotations could be extracted.


### Optional post-processing

After running the main workflow for multiple UniProt accessions, the optional script `reactome_post_process.py` can merge individual `result.csv` files. It  concatenates the available result tables, harmonizes the target columns when needed, removes duplicate rows, and writes three post-processed output files:


**merged_reaction.csv** contains a compact reaction-level summary across all analysed UniProt accessions. It keeps the target accession, target name, pathway hierarchy, disease annotation, lowest-level pathway, reaction name, reaction ID, reaction participants, and reaction direction. This file is useful for quickly comparing which reactions are associated with each protein across the full Reactome output.

**merged_highest_pathways.csv** contains a simplified pathway-level summary. It reports the highest-level Reactome pathways associated with each target protein, together with the corresponding Reactome pathway ID, UniProt accession, and target name. This file is useful for obtaining a non-redundant overview of the major biological areas covered by the analysed proteins.

**disease_single_sequence_site.csv** contains a filtered subset of the concatenated results. It keeps only rows belonging to the highest-level Reactome pathway Disease and where the sequence-site annotation corresponds to a single numeric residue position. This output is useful for downstream inspection of disease-associated reactions involving specific modified or annotated residue sites

---

## Input

Input for reactome_to_mavisp.py script:

| Argument | Description |
|---|---|
| `-u`, `--uniprot_ac` | UniProt accession to analyze. Default: `Q8N726`. Ignored when `--uniprot_file` is supplied. |
| `-uf`, `--uniprot_file` | Text file containing one UniProt accession per line. Blank lines and lines starting with `#` are ignored. |
| `-o`, `--output_dir` | Main output directory. Default: `reactome_outputs`. A subfolder is created for each accession. |
| `-s`, `--skip_pathway_order` | Skip pathway-order inference. This speeds up the analysis and avoids writing `pathways_order/` files. |

Here an example of file with a list of uniprot ac 


```text
P04637
Q8N726
Q9Y2X3
```
Input for reactome_post_process.py script:

Arguments:

| Argument | Description |
|---|---|
| `-i`, `--input_dir` | Main Reactome output directory containing UniProt-specific folders. |
| `-o`, `--output_dir` | Output directory for merged tables. Default: same as `--input_dir`. |
| `--result_filename` | Name of the result file inside each UniProt folder. Default: `result.csv`. |

---

## Output

By default, outputs are written under:

```text
reactome_outputs/
```

For a single accession such as `P04637`, the output structure is:

```text
reactome_outputs/
├── entries_not_in_reactome.csv        # Only created when at least one accession fails
└── P04637/
    ├── result.csv
    ├── skipped_reactions.csv          # Only created when reactions are skipped
    └── pathways_order/                # Only created when pathway ordering is enabled
        └── <Reactome_pathway_ID>/
            ├── graph_edges.csv
            ├── graph_nodes.csv
            └── ordered_paths.csv
```


### `result.csv`

`result.csv` is the main output of the workflow. Each row corresponds to a protein entry in a Reactome reaction involving the target protein.

Common columns include:

| Column | Description |
|---|---|
| `target_uniprot_ac` | Input UniProt accession. |
| `target_name` | Final target name. Preferentially from UniProt, otherwise from Reactome, otherwise the UniProt accession. |
| `highest_pathway` | Highest-level Reactome pathway. |
| `highest_pathway_id` | Reactome stable ID of the highest-level pathway. |
| `pathway_1`, `pathway_2`, ... | Intermediate pathway hierarchy levels. |
| `pathway_1_id`, `pathway_2_id`, ... | Reactome stable IDs of intermediate pathway levels. |
| `lowest_pathway` | Lowest-level pathway containing the reaction. |
| `lowest_pathway_id` | Reactome stable ID of the lowest-level pathway. |
| `reaction_name` | Reactome reaction name. |
| `reaction_id` | Reactome stable reaction ID. |
| `protein` | Protein entry parsed from BioPAX. |
| `uniprot_ac` | UniProt accession associated with the parsed protein entry. |
| `cellular_location` | Cellular location of the protein entry. |
| `SequenceInterval` | Sequence interval annotation, when available. |
| `SequenceSite` | Sequence site annotation, when available. |
| `Modification_type` | Protein modification annotation, when available. |
| `is_a_protein_family` | Boolean flag indicating whether the BioPAX entry represents a protein family. Protein-family rows are removed from the final output. |
| `complex_of` | Complexes in which the protein participates. |
| `stoichiometry` | Stoichiometric coefficient of the protein within the corresponding complex. |
| `member_physical_entity_of` | Parent protein family or physical entity, when available. |
| `reaction_Left` | Left-side participants of the biochemical reaction. |
| `reaction_Right` | Right-side participants of the biochemical reaction. |
| `reaction_Conversion_Direction` | Biochemical reaction directionality. |
| `Controller_of_reaction_*` | Controllers, activators, or inhibitors associated with the reaction. |
| `disease_name` | Disease annotation/link when available. |
| `ordered` | Boolean flag indicating whether the reaction could be ordered using pathway-ordering files. |

Additional columns can appear depending on the BioPAX content returned by Reactome.


### `skipped_reactions.csv`

This file is written inside a UniProt-specific output folder when one or more reactions could not be processed.

Typical reasons include:

| Reason | Meaning |
|---|---|
| `query_id_returned_none_after_retries` | Reactome did not return usable metadata for a reaction after retries. |
| `biopax_unavailable_after_retries` | The BioPAX model for that reaction could not be downloaded after retries. |

Common columns:

| Column | Description |
|---|---|
| `uniprot_ac` | UniProt accession being analyzed. |
| `reaction_id` | Reactome stable reaction ID. |
| `pathway_id` | Reactome stable ID of the pathway associated with the reaction. |
| `pathway_name` | Pathway name associated with the reaction. |
| `reason` | Reason why the reaction was skipped. |


### `entries_not_in_reactome.csv`

This file is written at the global output-directory level when at least one input accession does not produce a valid output.

Possible statuses:

| Status | Meaning |
|---|---|
| `not_found_in_reactome` | The UniProt accession had no mapped human Reactome pathways. |
| `no_valid_reactome_output` | Reactome contained mapped pathways for the accession, but no valid reactions remained after filtering. |

Common columns:

| Column | Description |
|---|---|
| `uniprot_ac` | UniProt accession. |
| `status` | Failure/status category. |
| `reason` | Explanation of why no final output was produced. |


### `pathways_order/`

This directory is created only when pathway ordering is enabled.

For each pathway, the workflow writes:

| File | Description |
|---|---|
| `graph_edges.csv` | Directed edges between reaction nodes. |
| `graph_nodes.csv` | Reaction nodes with `is_start` and `is_end` flags. |
| `ordered_paths.csv` | Ordered reaction paths inferred from the directed graph. |

When `--skip_pathway_order` is used, this directory is not created and reactions in `result.csv` are marked as unordered.

### `post process analysis`

The post-processing script writes:

| File | Description |
|---|---|
| `merged_reaction.csv` | Deduplicated reaction-level summary across all analyzed UniProt accessions. |
| `merged_highest_pathways.csv` | Deduplicated table of highest-level pathways per target. |
| `disease_single_sequence_site.csv` | Subset of disease-pathway rows where the sequence-site annotation is a single numeric residue position. |

---

## Run

### Run one UniProt accession

```bash
python reactome_to_mavisp.py -u P04637
```

### Run one UniProt accession and skip pathway ordering

```bash
python reactome_to_mavisp.py -u P04637 -s
```
Skipping pathway ordering is faster because it avoids building pathway-level reaction graphs.

### Run with input list 

```bash
python reactome_to_mavisp.py -uf uniprot_list.txt -o reactome_outputs
```
### Run post process analysis

```bash
python reactome_post_process.py -i reactome_outputs -o reactome_outputs/summary
```

---

##  example

```bash
# Single protein, faster run without pathway ordering
python reactome_to_mavisp.py -u P04637 -s

# Multiple proteins
python reactome_to_mavisp.py -uf uniprot_list.txt -o reactome_outputs -s

# Merge results
python reactome_to_mavisp.py -i reactome_post_process.py -o reactome_outputs/summary
```

