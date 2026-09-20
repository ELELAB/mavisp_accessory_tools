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

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


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
    return sorted(input_dir.glob(f"*/{result_filename}"))


def infer_target_uniprot_ac(result_file: Path) -> str:
    return result_file.parent.name


def harmonize_target_columns(
    df: pd.DataFrame,
    result_file: Path,
) -> pd.DataFrame:
    df = df.copy()

    if "target_uniprot_ac" not in df.columns:
        df["target_uniprot_ac"] = infer_target_uniprot_ac(result_file)

    if "target_name" not in df.columns:
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
    df = df.copy()

    for column in columns:
        if column not in df.columns:
            df[column] = ""

    return df[columns]


def find_sequence_site_column(df: pd.DataFrame) -> Optional[str]:
    candidates: Dict[str, str] = {
        str(col).lower(): str(col)
        for col in df.columns
    }

    for key in ["sequence_site", "sequencesite"]:
        if key in candidates:
            return candidates[key]

    return None


def is_single_numeric_sequence_site(value: object) -> bool:
    if pd.isna(value):
        return False

    value_str = str(value).strip()
    return bool(re.fullmatch(r"\d+", value_str))


def build_target_label(df: pd.DataFrame) -> pd.Series:
    return (
        df["target_uniprot_ac"].fillna("").astype(str).str.strip()
        + " | "
        + df["target_name"].fillna("").astype(str).str.strip()
    )


def save_highest_pathways_plot(
    concatenated_df: pd.DataFrame,
    output_dir: Path,
) -> None:
    required = {"target_uniprot_ac", "target_name", "highest_pathway"}
    if not required.issubset(concatenated_df.columns):
        print("[WARNING] Cannot create highest_pathways.pdf: missing required columns.")
        return

    plot_df = concatenated_df.copy()
    plot_df["target_label"] = build_target_label(plot_df)
    plot_df["highest_pathway"] = plot_df["highest_pathway"].fillna("").astype(str).str.strip()

    plot_df = plot_df[
        (plot_df["target_label"] != " | ")
        & (plot_df["highest_pathway"] != "")
    ][["target_label", "highest_pathway"]].drop_duplicates()

    if plot_df.empty:
        print("[WARNING] No data available for highest_pathways.pdf")
        return

    pathway_order = (
        plot_df["highest_pathway"]
        .value_counts()
        .index
        .tolist()
    )

    target_order = (
        plot_df["target_label"]
        .value_counts()
        .index
        .tolist()
    )

    x_map = {pathway: i for i, pathway in enumerate(pathway_order)}
    y_map = {target: i for i, target in enumerate(target_order)}

    x = plot_df["highest_pathway"].map(x_map)
    y = plot_df["target_label"].map(y_map)

    fig_width = max(8, len(pathway_order) * 0.6)
    fig_height = max(6, len(target_order) * 0.35)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.scatter(x, y, marker="s", s=80)

    ax.set_xticks(range(len(pathway_order)))
    ax.set_xticklabels(pathway_order, rotation=90)
    ax.set_yticks(range(len(target_order)))
    ax.set_yticklabels(target_order)

    ax.set_xlabel("Highest pathway")
    ax.set_ylabel("Target")
    ax.set_title("Target × highest pathway")
    ax.grid(True, linestyle=":", alpha=0.4)

    fig.tight_layout()
    out_file = output_dir / "highest_pathways.pdf"
    fig.savefig(out_file)
    plt.close(fig)

    print(f"[INFO] Output written: {out_file}")


def summarize_functional_status(values: pd.Series) -> str:
    """
    Aggregate functional-status annotations for one target/disease pair.
    """
    has_lof = False
    has_gof = False

    for value in values.dropna():
        value_str = str(value).strip().lower()

        if not value_str:
            continue

        normalized = re.sub(r"[_-]+", " ", value_str)
        normalized = re.sub(r"\s+", " ", normalized)

        if normalized == "lof" or "loss of function" in normalized:
            has_lof = True

        if normalized == "gof" or "gain of function" in normalized:
            has_gof = True

    if has_lof and has_gof:
        return "mixed"

    if has_lof:
        return "loss_of_function"

    if has_gof:
        return "gain_of_function"

    return "unknown"


def save_disease_pathways_plot(
    concatenated_df: pd.DataFrame,
    output_dir: Path,
) -> None:
    """
    Plot target x disease associations.

    X axis:
        disease_name

    Y axis:
        target UniProt accession | target protein name

    Symbol/color:
        LoF, GoF, mixed, or unknown.
    """
    required = {
        "target_uniprot_ac",
        "target_name",
        "disease_name",
        "functional_status",
    }

    if not required.issubset(concatenated_df.columns):
        missing = sorted(
            required.difference(concatenated_df.columns)
        )
        print(
            "[WARNING] Cannot create disease_targets.pdf. "
            f"Missing required columns: {', '.join(missing)}"
        )
        return

    plot_df = concatenated_df.copy()
    plot_df["target_label"] = build_target_label(plot_df)

    plot_df["disease_name"] = (
        plot_df["disease_name"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    plot_df = plot_df[
        (plot_df["target_label"] != " | ")
        & (plot_df["disease_name"] != "")
    ].copy()

    if plot_df.empty:
        print(
            "[WARNING] No disease annotations available "
            "for disease_targets.pdf"
        )
        return

    agg_df = (
        plot_df.groupby(
            [
                "target_label",
                "disease_name",
            ],
            dropna=False,
        )["functional_status"]
        .apply(summarize_functional_status)
        .reset_index(name="status")
    )

    if agg_df.empty:
        print(
            "[WARNING] No aggregated disease data available "
            "for disease_targets.pdf"
        )
        return

    disease_order = (
        agg_df["disease_name"]
        .value_counts()
        .index
        .tolist()
    )

    target_order = (
        agg_df["target_label"]
        .value_counts()
        .index
        .tolist()
    )

    x_map = {
        disease: index
        for index, disease in enumerate(disease_order)
    }

    y_map = {
        target: index
        for index, target in enumerate(target_order)
    }

    style_map = {
        "loss_of_function": {
            "marker": "o",
            "color": "tab:blue",
            "label": "LoF",
        },
        "gain_of_function": {
            "marker": "^",
            "color": "tab:red",
            "label": "GoF",
        },
        "mixed": {
            "marker": "s",
            "color": "tab:purple",
            "label": "Mixed",
        },
        "unknown": {
            "marker": "x",
            "color": "0.5",
            "label": "Unknown",
        },
    }

    fig_width = max(
        7.0,
        min(
            24.0,
            1.35 * len(disease_order) + 3.0,
        ),
    )

    fig_height = max(
        3.5,
        min(
            24.0,
            0.55 * len(target_order) + 2.5,
        ),
    )

    fig, ax = plt.subplots(
        figsize=(
            fig_width,
            fig_height,
        )
    )

    for status, style in style_map.items():
        subset = agg_df[
            agg_df["status"] == status
        ]

        if subset.empty:
            continue

        ax.scatter(
            subset["disease_name"].map(x_map),
            subset["target_label"].map(y_map),
            marker=style["marker"],
            s=95,
            c=style["color"],
            label=style["label"],
            zorder=3,
        )

    ax.set_xticks(
        range(len(disease_order))
    )

    ax.set_xticklabels(
        disease_order,
        rotation=45,
        ha="right",
    )

    ax.set_yticks(
        range(len(target_order))
    )

    ax.set_yticklabels(
        target_order
    )

    ax.set_xlabel("Disease")
    ax.set_ylabel("Target")
    ax.set_title("Target × disease")

    ax.set_xlim(
        -0.5,
        len(disease_order) - 0.5,
    )

    ax.set_ylim(
        -0.5,
        len(target_order) - 0.5,
    )

    ax.grid(
        True,
        linestyle=":",
        alpha=0.25,
        zorder=0,
    )

    ax.legend(
        title="Functional status",
        bbox_to_anchor=(1.02, 1.0),
        loc="upper left",
        borderaxespad=0.0,
    )

    fig.tight_layout()

    out_file = (
        output_dir
        / "disease_targets.pdf"
    )

    fig.savefig(
        out_file,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(
        f"[INFO] Output written: {out_file}"
    )


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
    merged_reaction_df.to_csv(merged_reaction_file, index=False)

    merged_highest_pathways_df = select_existing_columns(
        concatenated_df,
        HIGHEST_PATHWAY_COLUMNS
    ).drop_duplicates()

    merged_highest_pathways_file = output_dir / "merged_highest_pathways.csv"
    merged_highest_pathways_df.to_csv(merged_highest_pathways_file, index=False)

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

    save_highest_pathways_plot(concatenated_df, output_dir)
    save_disease_pathways_plot(concatenated_df, output_dir)


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