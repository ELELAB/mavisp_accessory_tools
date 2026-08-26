#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd


MERGED_REACTION_BASE_COLUMNS = [
    "target_uniprot_ac",
    "target_name",
    "highest_pathway",
    "highest_pathway_id",
    "disease_name",
]

MERGED_REACTION_TAIL_COLUMNS = [
    "lowest_pathway",
    "lowest_pathway_id",
    "reaction_name",
    "reaction_id",
    "reaction_Left",
    "reaction_Right",
    "reaction_Conversion_Direction",
]

HIGHEST_PATHWAY_COLUMNS = [
    "highest_pathway",
    "highest_pathway_id",
    "target_uniprot_ac",
    "target_name",
]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Post-process Reactome result.csv files generated inside "
            "UniProt-specific folders. The script creates: "
            "merged_reaction.csv, merged_highest_pathways.csv, and "
            "disease_single_sequence_site.csv."
        )
    )

    parser.add_argument(
        "-i",
        "--input_dir",
        required=True,
        type=str,
        help="Main Reactome output directory containing UniProt-specific folders."
    )

    parser.add_argument(
        "-o",
        "--output_dir",
        required=False,
        default=None,
        type=str,
        help=(
            "Directory where merged output files will be written. "
            "Default: same as --input_dir."
        )
    )

    parser.add_argument(
        "--result_filename",
        required=False,
        default="result.csv",
        type=str,
        help="Name of the result file inside each UniProt folder. Default: result.csv"
    )

    return parser.parse_args()


def find_result_files(input_dir: Path, result_filename: str) -> List[Path]:
    """
    Find result files directly inside first-level UniProt subfolders.

    Example:
        input_dir/P04637/result.csv
        input_dir/Q8N726/result.csv
    """
    return sorted(input_dir.glob(f"*/{result_filename}"))


def infer_target_uniprot_ac(result_file: Path) -> str:
    """Infer the target UniProt accession from the parent folder name."""
    return result_file.parent.name


def harmonize_target_columns(df: pd.DataFrame, result_file: Path) -> pd.DataFrame:
    """
    Ensure target_uniprot_ac and target_name are present.

    Priority for target_uniprot_ac:
        1. existing target_uniprot_ac column
        2. existing uniprot_ac column
        3. parent folder name

    Priority for target_name:
        1. existing target_name column
        2. existing protein column
        3. empty string
    """
    df = df.copy()

    if "target_uniprot_ac" not in df.columns:
        if "uniprot_ac" in df.columns:
            df["target_uniprot_ac"] = df["uniprot_ac"]
        else:
            df["target_uniprot_ac"] = infer_target_uniprot_ac(result_file)

    if "target_name" not in df.columns:
        if "protein" in df.columns:
            df["target_name"] = df["protein"]
        else:
            df["target_name"] = ""

    return df


def load_result_files(result_files: List[Path]) -> pd.DataFrame:
    dfs: List[pd.DataFrame] = []

    for result_file in result_files:
        try:
            df = pd.read_csv(result_file)

            if df.empty:
                print(f"[INFO] Skipping empty file: {result_file}")
                continue

            df = harmonize_target_columns(df, result_file)
            df["source_result_file"] = str(result_file)
            df["source_folder"] = result_file.parent.name
            dfs.append(df)

        except Exception as error:
            print(f"[WARNING] Could not read {result_file}: {error}")

    if not dfs:
        return pd.DataFrame()

    return pd.concat(dfs, ignore_index=True)


def get_pathway_columns(df: pd.DataFrame) -> List[str]:
    """
    Return pathway_1/pathway_1_id ... pathway_n/pathway_n_id in numeric order.
    """
    pathway_numbers = set()

    for column in df.columns:
        match = re.fullmatch(r"pathway_(\d+)(_id)?", str(column))
        if match:
            pathway_numbers.add(int(match.group(1)))

    pathway_columns: List[str] = []

    for number in sorted(pathway_numbers):
        pathway_col = f"pathway_{number}"
        pathway_id_col = f"pathway_{number}_id"

        if pathway_col in df.columns:
            pathway_columns.append(pathway_col)

        if pathway_id_col in df.columns:
            pathway_columns.append(pathway_id_col)

    return pathway_columns


def select_existing_columns(df: pd.DataFrame, columns: List[str]) -> pd.DataFrame:
    """
    Select requested columns, creating missing columns as empty strings.

    This makes the script robust if one result.csv lacks a specific column.
    """
    df = df.copy()

    for column in columns:
        if column not in df.columns:
            df[column] = ""

    return df[columns]


def find_sequence_site_column(df: pd.DataFrame) -> Optional[str]:
    """
    Find sequence_site column independently of case/style.

    It accepts:
        sequence_site
        SequenceSite
    """
    candidates: Dict[str, str] = {
        str(col).lower(): str(col)
        for col in df.columns
    }

    for key in ["sequence_site", "sequencesite"]:
        if key in candidates:
            return candidates[key]

    return None


def is_single_numeric_sequence_site(value: object) -> bool:
    """
    True only for a single integer-like residue position.

    Accepted:
        330
        "330"

    Rejected:
        "330_331"
        "330-331"
        "330;331"
        empty / NaN
    """
    if pd.isna(value):
        return False

    value_str = str(value).strip()

    return bool(re.fullmatch(r"\d+", value_str))


def write_outputs(concatenated_df: pd.DataFrame, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    pathway_columns = get_pathway_columns(concatenated_df)

    merged_reaction_columns = (
        MERGED_REACTION_BASE_COLUMNS
        + pathway_columns
        + MERGED_REACTION_TAIL_COLUMNS
    )

    merged_reaction_df = select_existing_columns(
        concatenated_df,
        merged_reaction_columns
    ).drop_duplicates()

    merged_reaction_file = output_dir / "merged_reaction.csv"

    merged_reaction_df.to_csv(
        merged_reaction_file,
        index=False
    )

    merged_highest_pathways_df = select_existing_columns(
        concatenated_df,
        HIGHEST_PATHWAY_COLUMNS
    ).drop_duplicates()

    merged_highest_pathways_file = output_dir / "merged_highest_pathways.csv"

    merged_highest_pathways_df.to_csv(
        merged_highest_pathways_file,
        index=False
    )

    sequence_site_column = find_sequence_site_column(concatenated_df)

    if sequence_site_column is None:
        print(
            "[WARNING] No sequence_site/SequenceSite column found. "
            "Writing empty disease_single_sequence_site.csv."
        )
        disease_single_sequence_site_df = concatenated_df.iloc[0:0].copy()

    else:
        disease_single_sequence_site_df = concatenated_df[
            (
                concatenated_df["highest_pathway"]
                .astype(str)
                .str.strip()
                .eq("Disease")
            )
            & concatenated_df[sequence_site_column].apply(
                is_single_numeric_sequence_site
            )
        ].copy().drop_duplicates()

    disease_single_sequence_site_file = output_dir / "disease_single_sequence_site.csv"

    disease_single_sequence_site_df.to_csv(
        disease_single_sequence_site_file,
        index=False
    )

    print(f"[INFO] Output written: {merged_reaction_file}")
    print(f"[INFO] Rows: {len(merged_reaction_df)}")

    print(f"[INFO] Output written: {merged_highest_pathways_file}")
    print(f"[INFO] Rows: {len(merged_highest_pathways_df)}")

    print(f"[INFO] Output written: {disease_single_sequence_site_file}")
    print(f"[INFO] Rows: {len(disease_single_sequence_site_df)}")


def main() -> None:
    args = parse_arguments()

    input_dir = Path(args.input_dir).resolve()

    if args.output_dir:
        output_dir = Path(args.output_dir).resolve()
    else:
        output_dir = input_dir

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")

    if not input_dir.is_dir():
        raise NotADirectoryError(f"Input path is not a directory: {input_dir}")

    result_files = find_result_files(
        input_dir=input_dir,
        result_filename=args.result_filename
    )

    if not result_files:
        print(f"[WARNING] No {args.result_filename} files found inside {input_dir}")
        return

    concatenated_df = load_result_files(result_files)

    if concatenated_df.empty:
        print("[WARNING] No valid result files could be loaded.")
        return

    if "highest_pathway" not in concatenated_df.columns:
        raise ValueError("Missing required column: highest_pathway")

    write_outputs(
        concatenated_df=concatenated_df,
        output_dir=output_dir
    )

    print(f"[INFO] Found result files: {len(result_files)}")
    print(
        "[INFO] Concatenated rows before filtering/deduplication: "
        f"{len(concatenated_df)}"
    )


if __name__ == "__main__":
    main()