import os
from collections import defaultdict
from typing import Any, Dict, List, Set, Tuple, Optional
import networkx as nx

import pandas as pd

from reactome_pipeline.context_factory import ReactomeEventContextFactory
from reactome_pipeline.data_processing import DataProcessingFunctions
from reactome_pipeline.disease_variants import DiseaseVariantIndex
from reactome_pipeline.legacy_adapter import LegacyEventAdapter
from reactome_pipeline.local_reactome import LocalReactomeDatabase
from reactome_pipeline.uniprot_utils import UniprotFunctions
from reactome_pipeline.graph_utils import PathwayGraphFunctions


class ReactomeScript:
    """
    Fully local Reactome -> MAVISp workflow.

    Reactome data are read from the local release files:

    - UniProt2Reactome_PE_Reactions.txt
    - Homo_sapiens.owl
    - disease_variant_ewas_mapping.tsv

    NORMAL and DISEASE events are built independently but use the same
    EventContextBuilder through ReactomeEventContextFactory.
    """

    def __init__(
        self,
        uniprot_ac: str,
        skip_pathway_order: bool = False,
        output_dir: str = ".",
        reaction_map_file: str = (
            "reactome_data/UniProt2Reactome_PE_Reactions.txt"
        ),
        biopax_file: str = "reactome_data/Homo_sapiens.owl",
        disease_variant_file: str = (
            "reactome_data/disease_variant_ewas_mapping.tsv"
        ),
        db: Optional[LocalReactomeDatabase] = None,
        disease_index: Optional[DiseaseVariantIndex] = None,
    ) -> None:
        self.uniprot_ac = uniprot_ac
        self.skip_pathway_order = skip_pathway_order
        self.output_dir = output_dir

        self.reaction_map_file = reaction_map_file
        self.biopax_file = biopax_file
        self.disease_variant_file = disease_variant_file

        os.makedirs(
            self.output_dir,
            exist_ok=True,
        )

        if db is None:
            db = LocalReactomeDatabase(
                reaction_map_file=self.reaction_map_file,
                biopax_file=self.biopax_file,
            )

        self.db = db

        if disease_index is None:
            disease_index = DiseaseVariantIndex(
                self.disease_variant_file
            )

        self.disease_index = disease_index

        self.context_factory = ReactomeEventContextFactory(
            self.db,
            self.disease_index,
        )

        self.legacy_adapter = LegacyEventAdapter()

    def output_path(self, *parts: str) -> str:
        return os.path.join(self.output_dir, *parts)

    def resolve_target_information(self) -> Dict[str, Any]:
        """
        Resolve target-level metadata without using Reactome web services.

        The protein name is obtained through the existing UniProt helper when
        available. Reactome event discovery itself remains fully local.
        """
        try:
            target_name = (
                UniprotFunctions.uniprot_ac_to_protein_name(
                    self.uniprot_ac
                )
            )
        except Exception as exc:
            print(
                "[WARNING] Could not resolve UniProt protein name for "
                f"{self.uniprot_ac}: {exc}"
            )
            target_name = None

        return {
            "target_name": target_name or self.uniprot_ac,
        }

    def build_local_event_contexts(
        self,
    ) -> Dict[str, List[Dict[str, Any]]]:
        """
        Build independent NORMAL and DISEASE EventContexts for the target.
        """
        print(
            f"[INFO] Building local Reactome contexts for "
            f"{self.uniprot_ac}"
        )

        contexts = self.context_factory.build_all_contexts(
            self.uniprot_ac
        )

        print(
            f"[INFO] NORMAL EventContexts: "
            f"{len(contexts['normal'])}"
        )
        print(
            f"[INFO] DISEASE EventContexts: "
            f"{len(contexts['disease'])}"
        )

        return contexts

    def build_legacy_input(
        self,
        contexts: Dict[str, List[Dict[str, Any]]],
    ) -> List[Dict[Any, Any]]:
        """
        Convert all EventContexts into the nested legacy representation used
        by DataProcessingFunctions.

        NORMAL and DISEASE contexts remain independent records.
        """
        all_contexts = (
            contexts.get("normal", [])
            + contexts.get("disease", [])
        )

        return self.legacy_adapter.adapt_many(
            all_contexts
        )

    def remove_duplicate_pathway_reactions(
        self,
        result_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Remove duplicated reaction rows caused by overlapping Reactome pathway
        hierarchies.

        The historical pathway-specificity logic is preserved. NORMAL and
        DISEASE events are kept independent by grouping on context type as
        well as reaction name when that column is available.
        """
        grouping_columns = ["reaction_name"]

        if "event_context_type" in result_df.columns:
            grouping_columns.insert(
                0,
                "event_context_type",
            )

        grouped = result_df.groupby(
            grouping_columns,
            dropna=False,
        )

        paths_to_remove: List[
            Tuple[str, str, str]
        ] = []

        for group_key in grouped.groups.keys():
            reaction_df = grouped.get_group(group_key)

            if isinstance(group_key, tuple):
                if len(grouping_columns) == 2:
                    context_type, reaction_name = group_key
                else:
                    context_type = ""
                    reaction_name = group_key[-1]
            else:
                context_type = ""
                reaction_name = group_key

            lowest_pathway_ids = list(
                set(
                    reaction_df[
                        "lowest_pathway_id"
                    ].dropna().tolist()
                )
            )
            reaction_ids = list(
                set(
                    reaction_df[
                        "reaction_id"
                    ].dropna().tolist()
                )
            )

            reaction_dict: Dict[
                str,
                Set[str],
            ] = defaultdict(set)

            if len(lowest_pathway_ids) != len(
                reaction_ids
            ):
                for lowest_id, reaction_id in zip(
                    reaction_df[
                        "lowest_pathway_id"
                    ].tolist(),
                    reaction_df[
                        "reaction_id"
                    ].tolist(),
                ):
                    if pd.notna(lowest_id) and pd.notna(
                        reaction_id
                    ):
                        reaction_dict[
                            reaction_id
                        ].add(lowest_id)

                reaction_with_multiple_pathways = {
                    key: list(value)
                    for key, value in reaction_dict.items()
                    if len(value) > 1
                }

                for (
                    reaction_id_val,
                    pathways,
                ) in (
                    reaction_with_multiple_pathways.items()
                ):
                    paths_to_check: List[
                        Tuple[str, str]
                    ] = []

                    for pathway in pathways:
                        mask = (
                            (
                                result_df[
                                    "reaction_name"
                                ]
                                == reaction_name
                            )
                            & (
                                result_df[
                                    "reaction_id"
                                ]
                                == reaction_id_val
                            )
                            & (
                                result_df[
                                    "lowest_pathway_id"
                                ]
                                == pathway
                            )
                        )

                        if (
                            "event_context_type"
                            in result_df.columns
                        ):
                            mask = mask & (
                                result_df[
                                    "event_context_type"
                                ]
                                == context_type
                            )

                        df_for_comparison = (
                            result_df[
                                mask
                            ].copy()
                        )

                        pathway_cols = [
                            col
                            for col in result_df.columns
                            if "pathway" in col
                            and "id" in col
                        ]

                        pathway_cols.append(
                            "reaction_id"
                        )

                        path_strings = set()

                        for _, row in (
                            df_for_comparison.iterrows()
                        ):
                            values = [
                                str(value)
                                for value in row[
                                    pathway_cols
                                ].tolist()
                                if pd.notna(value)
                            ]

                            path_strings.add(
                                "_".join(values)
                            )

                        for path_string in path_strings:
                            paths_to_check.append(
                                (
                                    pathway,
                                    path_string,
                                )
                            )

                    if len(paths_to_check) < 2:
                        continue

                    processed = []

                    for lowest_id, path_string in (
                        paths_to_check
                    ):
                        path_parts = [
                            value
                            for value in (
                                path_string.split("_")
                            )
                            if value != "nan"
                        ]

                        processed.append(
                            (
                                lowest_id,
                                path_parts,
                            )
                        )

                    shortest_lowest_id, _ = min(
                        processed,
                        key=lambda item: len(item[1]),
                    )

                    paths_to_remove.append(
                        (
                            str(context_type),
                            str(shortest_lowest_id),
                            str(reaction_id_val),
                        )
                    )

        if not paths_to_remove:
            return result_df

        remove_set = set(paths_to_remove)

        def should_remove(row: pd.Series) -> bool:
            context_type = str(
                row.get(
                    "event_context_type",
                    "",
                )
            )

            key = (
                context_type,
                str(row["lowest_pathway_id"]),
                str(row["reaction_id"]),
            )

            return key in remove_set

        return result_df[
            ~result_df.apply(
                should_remove,
                axis=1,
            )
        ].copy()

    def order_reactions_by_pathway_files(
        self,
        result_df_filtered: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Reorder reactions using locally generated pathway-order files.

        When pathway ordering is enabled, ordered_paths.csv files are generated
        from the local Reactome BioPAX release before this method is called.
        If no linear ordering is available for a pathway, rows are retained and
        ``ordered`` is set to False.
        """
        if self.skip_pathway_order:
            result = result_df_filtered.copy()
            result["ordered"] = False
            return result

        ordered_result = pd.DataFrame()

        for stid in set(
            result_df_filtered[
                "lowest_pathway_id"
            ].dropna().tolist()
        ):
            ordered_paths_file = self.output_path(
                "pathways_order",
                stid,
                "ordered_paths.csv",
            )

            if os.path.exists(
                ordered_paths_file
            ):
                ordered_paths_df = pd.read_csv(
                    ordered_paths_file
                )

                df_sorted = (
                    DataProcessingFunctions
                    .reorder_dataframe_from_ordered_paths(
                        result_df_filtered,
                        stid,
                        ordered_paths_df,
                    )
                )
            else:
                df_sorted = result_df_filtered[
                    result_df_filtered[
                        "lowest_pathway_id"
                    ]
                    == stid
                ].copy()

                df_sorted["ordered"] = False

            ordered_result = pd.concat(
                [
                    ordered_result,
                    df_sorted,
                ],
                ignore_index=True,
            )

        # Keep rows that do not have a lowest pathway ID.
        no_pathway = result_df_filtered[
            result_df_filtered[
                "lowest_pathway_id"
            ].isna()
        ].copy()

        if not no_pathway.empty:
            no_pathway["ordered"] = False

            ordered_result = pd.concat(
                [
                    ordered_result,
                    no_pathway,
                ],
                ignore_index=True,
            )

        return ordered_result

    def build_and_write_result(
        self,
        legacy_data: List[Dict[Any, Any]],
        target_info: Dict[str, Any],
    ) -> bool:
        """
        Flatten local EventContexts, clean pathway overlaps and write result.csv.
        """
        result_df = (
            DataProcessingFunctions.process_data(
                legacy_data
            )
        )

        if result_df.empty:
            print(
                f"[INFO] No Reactome reactions found for "
                f"{self.uniprot_ac}."
            )
            return False

        sorted_columns = sorted(
            result_df.columns,
            key=(
                DataProcessingFunctions
                .column_sort_key
            ),
        )

        result_df = result_df[
            sorted_columns
        ]

        # Propagate complex and stoichiometry information only within the same
        # biological event context and reaction.
        if "complex_of" in result_df.columns:
            group_columns = [
                "reaction_id",
            ]

            if (
                "event_context_type"
                in result_df.columns
            ):
                group_columns.insert(
                    0,
                    "event_context_type",
                )

            result_df = (
                result_df.groupby(
                    group_columns,
                    group_keys=False,
                    dropna=False,
                )
                .apply(
                    DataProcessingFunctions
                    .fill_complex_info
                )
            )

        # Keep historical behavior: protein-family/entity-set rows are not
        # included in the final result.csv.
        if (
            "is_a_protein_family"
            in result_df.columns
        ):
            result_df = result_df.loc[
                result_df[
                    "is_a_protein_family"
                ]
                == False
            ]

        if result_df.empty:
            print(
                f"[INFO] No Reactome reactions left for "
                f"{self.uniprot_ac} after removing "
                "protein families."
            )
            return False

        result_df_filtered = (
            self.remove_duplicate_pathway_reactions(
                result_df
            )
        )

        if not self.skip_pathway_order:

            pathway_ids = (
                result_df_filtered[
                    "lowest_pathway_id"
                ]
                .dropna()
                .astype(str)
                .unique()
                .tolist()
            )

            self.write_local_pathway_ordering(
                pathway_ids
            )

        ordered_result_df_filtered = (
            self.order_reactions_by_pathway_files(
                result_df_filtered
            )
        )

        if ordered_result_df_filtered.empty:
            print(
                f"[INFO] Final Reactome output is empty "
                f"for {self.uniprot_ac}."
            )
            return False

        ordered_result_df_filtered.insert(
            0,
            "target_uniprot_ac",
            self.uniprot_ac,
        )

        ordered_result_df_filtered.insert(
            1,
            "target_name",
            target_info["target_name"],
        )

        ordered_result_df_filtered.to_csv(
            self.output_path("result.csv"),
            sep=",",
            index=False,
        )

        print(
            f"[INFO] Wrote "
            f"{len(ordered_result_df_filtered)} rows to "
            f"{self.output_path('result.csv')}"
        )

        return True

    def write_local_pathway_ordering(
        self,
        pathway_ids: List[str],
    ) -> None:
        """
        Build pathway graphs and ordered paths entirely from local BioPAX data.
        """

        for pathway_id in sorted(set(pathway_ids)):

            if not pathway_id:
                continue

            print(
                f"[INFO] Building local pathway ordering for "
                f"{pathway_id}"
            )
            pathway_output_dir = self.output_path(
                "pathways_order",
                pathway_id,
            )

            os.makedirs(
                pathway_output_dir,
                exist_ok=True,
            )

            ordered_paths_file = os.path.join(
                pathway_output_dir,
                "ordered_paths.csv",
            )

            # Remove ordering from a previous run.
            # A new file will be written below only if the current
            # Reactome graph contains at least one linear path.
            if os.path.exists(ordered_paths_file):
                os.remove(ordered_paths_file)

            reaction_order = (
                self.db.get_pathway_reaction_order(
                    pathway_id
                )
            )

            if not reaction_order:
                print(
                    f"[INFO] No pathway-order information "
                    f"for {pathway_id}"
                )
                continue

            reactions_list: List[
                Tuple[str, str, str]
            ] = []

            for reaction_id, records in reaction_order.items():

                if not records:
                    reactions_list.append(
                        (
                            reaction_id,
                            "",
                            "",
                        )
                    )
                    continue

                for record in records:

                    next_steps = record.get(
                        "next_step",
                        [],
                    )

                    previous_steps = record.get(
                        "previous",
                        [],
                    )

                    if next_steps and previous_steps:

                        for next_id in next_steps:
                            for previous_id in previous_steps:
                                reactions_list.append(
                                    (
                                        reaction_id,
                                        next_id,
                                        previous_id,
                                    )
                                )

                    elif next_steps:

                        for next_id in next_steps:
                            reactions_list.append(
                                (
                                    reaction_id,
                                    next_id,
                                    "",
                                )
                            )

                    elif previous_steps:

                        for previous_id in previous_steps:
                            reactions_list.append(
                                (
                                    reaction_id,
                                    "",
                                    previous_id,
                                )
                            )

                    else:

                        reactions_list.append(
                            (
                                reaction_id,
                                "",
                                "",
                            )
                        )

            (
                graph,
                starting_nodes,
                ending_nodes,
                edge_rows,
                node_rows,
            ) = PathwayGraphFunctions.make_pathway_graph(
                reactions_list,
                pathway_id,
            )


            pd.DataFrame(
                edge_rows
            ).to_csv(
                os.path.join(
                    pathway_output_dir,
                    "graph_edges.csv",
                ),
                index=False,
            )

            pd.DataFrame(
                node_rows
            ).to_csv(
                os.path.join(
                    pathway_output_dir,
                    "graph_nodes.csv",
                ),
                index=False,
            )

            simple_paths: List[List[str]] = []

            for start_node in starting_nodes:
                for end_node in ending_nodes:

                    if start_node == end_node:
                        continue

                    if (
                        start_node not in graph
                        or end_node not in graph
                    ):
                        continue

                    try:
                        paths = nx.all_simple_paths(
                            graph,
                            source=start_node,
                            target=end_node,
                        )

                        simple_paths.extend(
                            list(paths)
                        )

                    except nx.NetworkXNoPath:
                        continue

            paths_to_process = (
                PathwayGraphFunctions
                .remove_duplicates_order(
                    simple_paths
                )
            )

            filtered_paths = (
                PathwayGraphFunctions
                .find_non_subsequences(
                    paths_to_process
                )
            )

            if not filtered_paths:
                print(
                    f"[INFO] No linear start-to-end paths "
                    f"for {pathway_id}"
                )
                continue

            ordered_paths_rows = []

            for path_id, path in enumerate(
                filtered_paths
            ):
                for step, reaction_id in enumerate(path):

                    ordered_paths_rows.append(
                        {
                            "pathway_id": pathway_id,
                            "path_id": path_id,
                            "step": step,
                            "reaction_id": reaction_id,
                        }
                    )

            pd.DataFrame(
                ordered_paths_rows
            ).to_csv(
                ordered_paths_file,
                index=False,
            )

    def run(self) -> bool:
        """
        Execute the fully local Reactome analysis workflow.
        """
        target_info = (
            self.resolve_target_information()
        )

        contexts = (
            self.build_local_event_contexts()
        )

        if not contexts["all"]:
            print(
                f"[INFO] No local Reactome events found "
                f"for {self.uniprot_ac}."
            )
            return False

        legacy_data = (
            self.build_legacy_input(contexts)
        )

        return self.build_and_write_result(
            legacy_data,
            target_info,
        )
