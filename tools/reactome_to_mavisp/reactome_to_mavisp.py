#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (C) 2026, Matteo Arnaudi  <mata@cancer.dk>,<matarn@dtu.dk>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""
Fully local Reactome analysis workflow for UniProt accessions.

Given one UniProt accession, or a file containing multiple UniProt accessions,
the script uses local Reactome release files to retrieve NORMAL and DISEASE
events, extracts BioPAX-level annotations, and writes cleaned CSV outputs.

Required local Reactome files
-----------------------------
- reactome_data/UniProt2Reactome_PE_Reactions.txt
- reactome_data/Homo_sapiens.owl
- reactome_data/disease_variant_ewas_mapping.tsv

Main outputs
------------
For each UniProt accession, the script creates a dedicated output folder
containing:

- result.csv
    Final cleaned table of Reactome reactions involving the target protein.

At the global output level, the script can also write:

- entries_not_in_reactome.csv
    UniProt accessions that did not produce a valid local Reactome output.

Usage
-----
Run the workflow for one UniProt accession:

    python reactome_to_mavisp.py -u P04637 -s

Run the workflow for multiple UniProt accessions:

    python reactome_to_mavisp.py -uf uniprot_list.txt -s
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
from typing import Dict, List

import pandas as pd
import urllib.request
import zipfile
from reactome_pipeline.local_reactome import LocalReactomeDatabase
from reactome_pipeline.disease_variants import DiseaseVariantIndex
from reactome_pipeline.workflow import ReactomeScript
from datetime import datetime


DEFAULT_REACTION_MAP = (
    "reactome_data/UniProt2Reactome_PE_Reactions.txt"
)

DEFAULT_BIOPAX = (
    "reactome_data/Homo_sapiens.owl"
)

DEFAULT_DISEASE_VARIANTS = (
    "reactome_data/disease_variant_ewas_mapping.tsv"
)

REACTOME_DOWNLOAD_BASE = (
    "https://reactome.org/download/current"
)

REACTOME_REACTION_MAP_URL = (
    f"{REACTOME_DOWNLOAD_BASE}/"
    "UniProt2Reactome_PE_Reactions.txt"
)

REACTOME_DISEASE_VARIANTS_URL = (
    f"{REACTOME_DOWNLOAD_BASE}/"
    "disease_variant_ewas_mapping.tsv"
)

REACTOME_BIOPAX_URL = (
    f"{REACTOME_DOWNLOAD_BASE}/biopax.zip"
)

REACTOME_VERSION_URL = (
    "https://reactome.org/ContentService/data/database/version"
)

def get_file_date(path: str) -> str:
    timestamp = os.path.getmtime(path)
    return datetime.fromtimestamp(timestamp).strftime("%Y-%m-%d")

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Run the fully local Reactome-to-MAVISp workflow for one "
            "or more UniProt accessions."
        )
    )

    parser.add_argument(
        "-uf",
        "--uniprot_file",
        dest="uniprot_file",
        required=False,
        type=str,
        help="Text file containing one UniProt accession per line.",
    )

    parser.add_argument(
        "-o",
        "--output_dir",
        dest="output_dir",
        default="reactome_outputs",
        type=str,
        help="Main output directory.",
    )

    parser.add_argument(
        "-u",
        "--uniprot_ac",
        dest="uniprot_ac",
        default="Q8N726",
        type=str,
        help=(
            "UniProt accession of the protein of interest "
            "(default: Q8N726)."
        ),
    )

    parser.add_argument(
        "-s",
        "--skip_pathway_order",
        dest="skip_pathway_order",
        action="store_true",
        help=(
            "Skip pathway ordering. Recommended while using the "
            "fully local workflow unless local ordered_paths.csv "
            "files already exist."
        ),
    )

    parser.add_argument(
        "-r",
        "--refresh_reactome_data",
        dest="refresh_reactome_data",
        action="store_true",
        help=(
            "Download the current Reactome release files again, "
            "replacing the local reaction mapping, BioPAX, and "
            "disease-variant files before running the analysis."
        ),
    )

    parser.add_argument(
        "--reaction_map_file",
        default=DEFAULT_REACTION_MAP,
        type=str,
        help=(
            "Path to UniProt2Reactome_PE_Reactions.txt "
            f"(default: {DEFAULT_REACTION_MAP})."
        ),
    )

    parser.add_argument(
        "--biopax_file",
        default=DEFAULT_BIOPAX,
        type=str,
        help=(
            "Path to Homo_sapiens.owl "
            f"(default: {DEFAULT_BIOPAX})."
        ),
    )

    parser.add_argument(
        "--disease_variant_file",
        default=DEFAULT_DISEASE_VARIANTS,
        type=str,
        help=(
            "Path to disease_variant_ewas_mapping.tsv "
            f"(default: {DEFAULT_DISEASE_VARIANTS})."
        ),
    )

    parser.add_argument(
        "--reactome_release",
        default="unknown",
        type=str,
        help=(
            "Reactome release used to build the local database "
            "(for example: 94). Stored in metadata.json."
        ),
    )

    parser.add_argument(
        "--reactome_download_date",
        default="unknown",
        type=str,
        help=(
            "Date on which the local Reactome files were downloaded, "
            "preferably in YYYY-MM-DD format. Stored in metadata.json."
        ),
    )

    return parser.parse_args()


def read_uniprot_accessions(
    uniprot_file: str,
) -> List[str]:
    """Read one UniProt accession per non-empty, non-comment line."""
    accessions: List[str] = []

    with open(uniprot_file) as handle:
        for line in handle:
            accession = line.strip()

            if not accession:
                continue

            if accession.startswith("#"):
                continue

            accessions.append(accession)

    # Preserve input order while removing duplicates.
    return list(dict.fromkeys(accessions))

def get_current_reactome_release() -> str:
    """
    Retrieve the current Reactome database release number
    from the Reactome Content Service.
    """
    print(
        "[INFO] Retrieving current Reactome release..."
    )

    request = urllib.request.Request(
        REACTOME_VERSION_URL,
        headers={
            "Accept": "text/plain",
            "User-Agent": (
                "reactome-to-mavisp/1.0"
            ),
        },
    )

    try:
        with urllib.request.urlopen(
            request,
            timeout=30,
        ) as response:

            release = (
                response
                .read()
                .decode("utf-8")
                .strip()
            )

    except Exception as exc:
        print(
            "[WARNING] Could not retrieve current "
            f"Reactome release: {exc}"
        )
        return "unknown"

    if not release:
        print(
            "[WARNING] Reactome release endpoint "
            "returned an empty response."
        )
        return "unknown"

    print(
        f"[INFO] Current Reactome release: {release}"
    )

    return release


def validate_local_reactome_files(
    reaction_map_file: str,
    biopax_file: str,
    disease_variant_file: str,
) -> None:
    """Fail early when one of the required local Reactome files is missing."""
    required_files = {
        "Reactome UniProt/reaction mapping": reaction_map_file,
        "Reactome Homo sapiens BioPAX": biopax_file,
        "Reactome disease variant mapping": disease_variant_file,
    }

    missing = [
        f"{label}: {path}"
        for label, path in required_files.items()
        if not os.path.isfile(path)
    ]

    if missing:
        raise FileNotFoundError(
            "Missing required local Reactome file(s):\n- "
            + "\n- ".join(missing)
        )


def download_file(
    url: str,
    destination: str,
) -> None:
    """
    Download a file atomically.

    The download is first written to <destination>.tmp and replaces the
    existing file only after the transfer has completed successfully.
    """
    os.makedirs(
        os.path.dirname(destination) or ".",
        exist_ok=True,
    )

    temporary_file = destination + ".tmp"

    print(
        f"[INFO] Downloading:\n"
        f"       {url}\n"
        f"       -> {destination}"
    )

    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": (
                "reactome-to-mavisp/1.0 "
                "(local Reactome data downloader)"
            )
        },
    )

    try:
        with urllib.request.urlopen(request) as response:
            with open(temporary_file, "wb") as handle:
                shutil.copyfileobj(
                    response,
                    handle,
                    length=1024 * 1024,
                )

        os.replace(
            temporary_file,
            destination,
        )

    except Exception:
        if os.path.exists(temporary_file):
            os.remove(temporary_file)
        raise

def write_reactome_metadata(
    output_dir,
    reactome_release,
    reactome_download_date,
    reaction_map_file,
    biopax_file,
    disease_variant_file,
):
    os.makedirs(output_dir, exist_ok=True)

    if reactome_download_date in (None, "", "unknown"):
        reactome_download_date = get_file_date(biopax_file)

    metadata = {
        "reactome_release": reactome_release,
        "reactome_download_date": reactome_download_date,
        "biopax_file": os.path.basename(biopax_file),
        "reaction_mapping_file": os.path.basename(reaction_map_file),
        "disease_mapping_file": os.path.basename(disease_variant_file),
    }

    with open(
        os.path.join(output_dir, "metadata.json"),
        "w",
    ) as handle:
        json.dump(
            metadata,
            handle,
            indent=2,
            sort_keys=True,
        )

def refresh_reactome_data(
    reaction_map_file: str,
    biopax_file: str,
    disease_variant_file: str,
) -> None:
    """
    Download the current Reactome files required by the pipeline.

    UniProt2Reactome_PE_Reactions.txt and
    disease_variant_ewas_mapping.tsv are downloaded directly.

    Homo_sapiens.owl is extracted from Reactome's biopax.zip archive.
    """

    print("[INFO] Refreshing local Reactome data...")

    # ---------------------------------------------------------
    # UniProt -> Reactome reaction mapping
    # ---------------------------------------------------------

    download_file(
        REACTOME_REACTION_MAP_URL,
        reaction_map_file,
    )

    # ---------------------------------------------------------
    # Disease variant mapping
    # ---------------------------------------------------------

    download_file(
        REACTOME_DISEASE_VARIANTS_URL,
        disease_variant_file,
    )

    # ---------------------------------------------------------
    # BioPAX
    # ---------------------------------------------------------

    biopax_directory = (
        os.path.dirname(biopax_file)
        or "."
    )

    os.makedirs(
        biopax_directory,
        exist_ok=True,
    )

    archive_file = os.path.join(
        biopax_directory,
        "biopax.zip",
    )

    download_file(
        REACTOME_BIOPAX_URL,
        archive_file,
    )

    print(
        "[INFO] Extracting Homo_sapiens.owl "
        "from biopax.zip..."
    )

    temporary_biopax = biopax_file + ".tmp"

    try:
        with zipfile.ZipFile(
            archive_file,
            "r",
        ) as archive:

            candidates = [
                member
                for member in archive.namelist()
                if os.path.basename(member)
                == "Homo_sapiens.owl"
            ]

            if not candidates:
                raise FileNotFoundError(
                    "Homo_sapiens.owl was not found "
                    "inside the downloaded Reactome "
                    "biopax.zip archive."
                )

            if len(candidates) > 1:
                print(
                    "[WARNING] Multiple Homo_sapiens.owl "
                    "entries found in biopax.zip. "
                    f"Using: {candidates[0]}"
                )

            member = candidates[0]

            with archive.open(
                member,
                "r",
            ) as source:
                with open(
                    temporary_biopax,
                    "wb",
                ) as destination:
                    shutil.copyfileobj(
                        source,
                        destination,
                        length=1024 * 1024,
                    )

        os.replace(
            temporary_biopax,
            biopax_file,
        )

    except Exception:
        if os.path.exists(temporary_biopax):
            os.remove(temporary_biopax)
        raise

    finally:
        if os.path.exists(archive_file):
            os.remove(archive_file)

    print(
        "[INFO] Reactome data refresh completed."
    )

def write_missing_reactome_entries(
    output_dir: str,
    missing_entries: List[Dict[str, str]],
) -> None:
    """Write UniProt accessions not producing a valid Reactome output."""
    if not missing_entries:
        return

    os.makedirs(output_dir, exist_ok=True)

    output_file = os.path.join(
        output_dir,
        "entries_not_in_reactome.csv",
    )

    new_df = pd.DataFrame(missing_entries)

    if os.path.exists(output_file):
        old_df = pd.read_csv(output_file)
        final_df = pd.concat(
            [old_df, new_df],
            ignore_index=True,
        )
        final_df = final_df.drop_duplicates(
            subset=["uniprot_ac", "status"]
        )
    else:
        final_df = new_df.drop_duplicates(
            subset=["uniprot_ac", "status"]
        )

    final_df.to_csv(
        output_file,
        index=False,
    )


def main() -> None:
    args = parse_arguments()

    if args.refresh_reactome_data:

        reactome_release = (
            get_current_reactome_release()
        )

        refresh_reactome_data(
            reaction_map_file=args.reaction_map_file,
            biopax_file=args.biopax_file,
            disease_variant_file=args.disease_variant_file,
        )

        reactome_download_date = (
            datetime.now().strftime("%Y-%m-%d")
        )

    else:
        reactome_release = (
            args.reactome_release
        )

        reactome_download_date = (
            args.reactome_download_date
        )

    validate_local_reactome_files(
        reaction_map_file=args.reaction_map_file,
        biopax_file=args.biopax_file,
        disease_variant_file=args.disease_variant_file,
    )

    write_reactome_metadata(
        output_dir=args.output_dir,
        reactome_release=reactome_release,
        reactome_download_date=reactome_download_date,
        reaction_map_file=args.reaction_map_file,
        biopax_file=args.biopax_file,
        disease_variant_file=args.disease_variant_file,
    )

    if args.uniprot_file:
        uniprot_accessions = read_uniprot_accessions(
            args.uniprot_file
        )
    else:
        uniprot_accessions = [
            args.uniprot_ac
        ]

    print(
        "[INFO] Initializing shared local Reactome database..."
    )

    shared_db = LocalReactomeDatabase(
        reaction_map_file=args.reaction_map_file,
        biopax_file=args.biopax_file,
    )

    shared_disease_index = DiseaseVariantIndex(
        args.disease_variant_file
    )

    print(
        "[INFO] Loading and indexing shared BioPAX data..."
    )

    shared_db.load_biopax_model()
    shared_db.build_event_index()
    shared_db.build_pathway_index()
    shared_db.build_pathway_parent_index()

    missing_reactome_entries: List[
        Dict[str, str]
    ] = []

    for uniprot_ac in uniprot_accessions:
        print(
            f"[INFO] Running local Reactome analysis for "
            f"{uniprot_ac}"
        )

        protein_output_dir = os.path.join(
            args.output_dir,
            uniprot_ac,
        )

        workflow = ReactomeScript(
            uniprot_ac=uniprot_ac,
            skip_pathway_order=args.skip_pathway_order,
            output_dir=protein_output_dir,
            reaction_map_file=args.reaction_map_file,
            biopax_file=args.biopax_file,
            disease_variant_file=args.disease_variant_file,
            db=shared_db,
            disease_index=shared_disease_index,
        )

        has_valid_output = workflow.run()

        if not has_valid_output:
            print(
                f"[INFO] No valid local Reactome output "
                f"for {uniprot_ac}."
            )

            if os.path.isdir(
                protein_output_dir
            ):
                shutil.rmtree(
                    protein_output_dir
                )

            missing_reactome_entries.append(
                {
                    "uniprot_ac": uniprot_ac,
                    "status": (
                        "no_valid_reactome_output"
                    ),
                    "reason": (
                        "No NORMAL or DISEASE Reactome "
                        "events produced a valid final output "
                        "from the local release."
                    ),
                }
            )

    write_missing_reactome_entries(
        output_dir=args.output_dir,
        missing_entries=missing_reactome_entries,
    )


if __name__ == "__main__":
    main()
