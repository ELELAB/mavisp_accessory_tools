"""
Reactome analysis workflow for UniProt accessions.

Given one UniProt accession, or a file containing multiple UniProt accessions,
the script queries Reactome to retrieve associated pathways, identifies
candidate reactions involving the target protein, extracts BioPAX-level
protein/complex/reaction annotations, and writes cleaned CSV outputs.

Main outputs
------------
For each UniProt accession, the script creates a dedicated output folder
containing:

- result.csv
    Final cleaned table of Reactome reactions involving the target protein.

- skipped_reactions.csv
    Reactions that could not be processed because Reactome or BioPAX returned
    no usable data after retries.

- pathways_order/
    Optional pathway-ordering files generated when pathway ordering is not
    skipped.

At the global output level, the script can also write:

- entries_not_in_reactome.csv
    UniProt accessions that were not found in Reactome or did not produce a
    valid final output.

Usage
-----
Run the workflow for one UniProt accession:

    python reac_classified.py -u P04637

Run the workflow for multiple UniProt accessions:

    python reac_classified.py -uf uniprot_list.txt

Skip pathway ordering:

    python reac_classified.py -u P04637 -s
"""
from __future__ import annotations

import os
import argparse
import shutil
from typing import Dict, List

import pandas as pd

try:
    import reactome2py as rc
    from reactome2py import analysis, content, utils
except ImportError:
    rc = None  # type: ignore

from reactome_pipeline.reactome_analysis import ReactomeAnalysisFunctions
from reactome_pipeline.workflow import ReactomeScript
    
def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            'This script retrieves pathways associated with a UniProt accession from Reactome, '
            'orders the reactions, collects protein level data and writes several CSV outputs. '
            'It is a class‑based refactoring of the original reac.py script.'
        )
    )
    parser.add_argument(
        '-uf', '--uniprot_file',
        dest='uniprot_file',
        required=False,
        type=str,
        help='Text file containing one UniProt accession per line.'
        )

    parser.add_argument(
        '-o', '--output_dir',
        dest='output_dir',
        default='reactome_outputs',
        type=str,
        help='Main output directory.'
    )
    parser.add_argument(
        '-u', '--uniprot_ac', dest='uniprot_ac', default='Q8N726', type=str,
        help='UniProt accession of the protein of interest (default: Q8N726).'
    )
    parser.add_argument(
        '-s', '--skip_pathway_order', dest='skip_pathway_order', action='store_true',
        help='Skip the pathway ordering phase (analogous to supplying --pdb_file in the original script).'
    )
    return parser.parse_args()

def read_uniprot_accessions(uniprot_file: str) -> List[str]:
    accessions: List[str] = []

    with open(uniprot_file) as handle:
        for line in handle:
            accession = line.strip()

            if not accession:
                continue

            if accession.startswith("#"):
                continue

            accessions.append(accession)

    return accessions

def uniprot_has_reactome_entry(uniprot_ac: str) -> bool:
    """
    Return True if Reactome contains at least one mapped pathway for the UniProt accession.
    """
    if rc is None:
        raise ImportError("reactome2py is required to query Reactome.")

    reactome_pathways = ReactomeAnalysisFunctions.safe_reactome_call(
        rc.content.mapping,
        id=uniprot_ac,
        resource="UniProt",
        species="9606",
        by="pathways",
        default=[]
    )

    return bool(reactome_pathways)


def write_missing_reactome_entries(
    output_dir: str,
    missing_entries: List[Dict[str, str]]
) -> None:
    """
    Write UniProt accessions not producing a valid Reactome output.
    """
    if not missing_entries:
        return

    os.makedirs(output_dir, exist_ok=True)

    output_file = os.path.join(output_dir, "entries_not_in_reactome.csv")

    new_df = pd.DataFrame(missing_entries)

    if os.path.exists(output_file):
        old_df = pd.read_csv(output_file)
        final_df = pd.concat([old_df, new_df], ignore_index=True)
        final_df = final_df.drop_duplicates(subset=["uniprot_ac", "status"])
    else:
        final_df = new_df.drop_duplicates(subset=["uniprot_ac", "status"])

    final_df.to_csv(output_file, index=False)

def main() -> None:
    args = parse_arguments()

    if args.uniprot_file:
        uniprot_accessions = read_uniprot_accessions(args.uniprot_file)
    else:
        uniprot_accessions = [args.uniprot_ac]

    missing_reactome_entries: List[Dict[str, str]] = []

    for uniprot_ac in uniprot_accessions:
        print(f"[INFO] Running Reactome analysis for {uniprot_ac}")

        protein_output_dir = os.path.join(args.output_dir, uniprot_ac)

        # First check: is this UniProt accession known by Reactome?
        if not uniprot_has_reactome_entry(uniprot_ac):
            print(f"[INFO] {uniprot_ac} is not present in Reactome.")

            missing_reactome_entries.append({
                "uniprot_ac": uniprot_ac,
                "status": "not_found_in_reactome",
                "reason": "Reactome search_fireworks returned no entries"
            })

            continue

        workflow = ReactomeScript(
            uniprot_ac=uniprot_ac,
            skip_pathway_order=args.skip_pathway_order,
            output_dir=protein_output_dir
        )

        has_valid_output = workflow.run()

        if not has_valid_output:
            print(f"[INFO] Removing empty output folder for {uniprot_ac}")

            if os.path.isdir(protein_output_dir):
                shutil.rmtree(protein_output_dir)

            missing_reactome_entries.append({
                "uniprot_ac": uniprot_ac,
                "status": "no_valid_reactome_output",
                "reason": "Reactome entry found, but no valid reactions remained after filtering"
            })

    write_missing_reactome_entries(
        output_dir=args.output_dir,
        missing_entries=missing_reactome_entries
    )


if __name__ == '__main__':
    main()