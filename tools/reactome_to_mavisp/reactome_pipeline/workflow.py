import os
from collections import defaultdict
from typing import Any, Dict, List, Set, Tuple


import pandas as pd
import networkx as nx

try:
    import reactome2py as rc
    from reactome2py import analysis, content, utils
except ImportError:
    rc = None  # type: ignore

try:
    import pybiopax
except ImportError:
    pybiopax = None  # type: ignore

from reactome_pipeline.graph_utils import PathwayGraphFunctions
from reactome_pipeline.data_processing import DataProcessingFunctions
from reactome_pipeline.reactome_analysis import ReactomeAnalysisFunctions
from reactome_pipeline.uniprot_utils import UniprotFunctions

class ReactomeScript:
    """Encapsulate the original script's workflow in an object oriented form."""

    def __init__(
        self,
        uniprot_ac: str,
        skip_pathway_order: bool = False,
        output_dir: str = "."
    ) -> None:
        self.uniprot_ac = uniprot_ac
        self.skip_pathway_order = skip_pathway_order
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

    def output_path(self, *parts: str) -> str:
        return os.path.join(self.output_dir, *parts)

    def retrieve_reactome_pathways(self) -> List[Dict[str, Any]]:
        return ReactomeAnalysisFunctions.safe_reactome_call(
            rc.content.mapping,
            id=self.uniprot_ac,
            resource="UniProt",
            species="9606",
            by="pathways",
            default=[]
        )
    def resolve_target_information(self) -> Dict[str, Any]:
        """
        Retrieve target-related information.

        search_fireworks is used only as an optional extra check.
        If it fails, the workflow continues using:
        1. UniProt-derived protein name.
        2. Reactome complexes associated with the UniProt accession.
        3. Direct UniProt accession matching inside BioPAX reaction models.
        """
        target_protein_stId: List[str] = []
        target_protein_stId_other: List[str] = []

        # Primary target name source: UniProt.
        uniprot_target_name = UniprotFunctions.uniprot_ac_to_protein_name(
            self.uniprot_ac
        )

        # Optional target name source: Reactome search_fireworks.
        reactome_target_name = None

        target_protein_info = ReactomeAnalysisFunctions.safe_reactome_call(
            rc.content.search_fireworks,
            query=self.uniprot_ac,
            default={"entries": []}
        )

        entries = target_protein_info.get("entries", []) if target_protein_info else []

        if not entries:
            print(
                f"[WARNING] search_fireworks failed or returned no entries for "
                f"{self.uniprot_ac}. Continuing without Reactome protein-form IDs."
            )

        for entry in entries:
            if entry.get("stId"):
                target_protein_stId.append(entry["stId"])

            if reactome_target_name is None:
                reactome_target_name = (
                    entry.get("displayName")
                    or entry.get("name")
                )

        # Final target name priority:
        # 1. UniProt protein name
        # 2. Reactome display/name
        # 3. UniProt accession fallback
        target_name = (
            uniprot_target_name
            or reactome_target_name
            or self.uniprot_ac
        )

        # Optional extra check: alternative Reactome forms of the target protein.
        for prot_id in target_protein_stId:
            other_forms = ReactomeAnalysisFunctions.safe_reactome_call(
                rc.content.entity_other_form,
                id=prot_id,
                default=[]
            )

            if other_forms:
                for form in other_forms:
                    if form.get("stId"):
                        target_protein_stId_other.append(form["stId"])

        all_target_protein_stId = list(
            set(target_protein_stId + target_protein_stId_other)
        )

        # Main useful information: complexes associated with the UniProt accession.
        target_protein_complexes_names: List[str] = []

        complexes_list = ReactomeAnalysisFunctions.safe_reactome_call(
            rc.content.entities_complexes,
            id=self.uniprot_ac,
            resource="UniProt",
            default=[]
        )

        for complex_entry in complexes_list:
            if complex_entry.get("stId"):
                target_protein_complexes_names.append(complex_entry["stId"])

        return {
            "target_name": target_name,
            "target_protein_stId": target_protein_stId,
            "target_protein_stId_other": target_protein_stId_other,
            "all_target_protein_stId": all_target_protein_stId,
            "target_protein_complexes_names": target_protein_complexes_names,
        }
    def write_pathway_ordering(self, reactome_pathways: List[Dict[str, Any]]) -> None:
        """
        Build and write reaction ordering information for each Reactome pathway.
        """
        reactions_order: Dict[str, Dict[str, List[Dict[str, List[str]]]]] = {}

        for i in reactome_pathways:
            biopax_model = ReactomeAnalysisFunctions.safe_model_from_reactome_cached(i['stId'])

            if biopax_model is None:
                print(f"[WARNING] Skipping pathway order for {i['stId']}: BioPAX unavailable.")
                continue

            reaction_order = ReactomeAnalysisFunctions.get_pathway_order(biopax_model)
            reactions_order[i['stId']] = reaction_order

        for hrsa, reac in reactions_order.items():
            reactions_lists: List[Tuple[str, str, str]] = []

            for key, value_list in reac.items():
                for value in value_list:
                    next_steps = value.get('next_step', [])
                    previous_steps = value.get('previous', [])

                    if next_steps and previous_steps:
                        for next_step in next_steps:
                            for next_item in next_step:
                                for previous_step in previous_steps:
                                    for prev_item in previous_step:
                                        reactions_lists.append((key, next_item, prev_item))

                    if not next_steps and previous_steps:
                        for previous_step in previous_steps:
                            for prev_item in previous_step:
                                reactions_lists.append((key, "", prev_item))

                    if next_steps and not previous_steps:
                        for next_step in next_steps:
                            for next_item in next_step:
                                reactions_lists.append((key, next_item, ""))

            G, starting_nodes, ending_nodes, edge_rows, node_rows = PathwayGraphFunctions.make_pathway_graph(
                reactions_lists,
                hrsa
            )

            pathway_output_dir = self.output_path("pathways_order", hrsa)
            os.makedirs(pathway_output_dir, exist_ok=True)

            pd.DataFrame(edge_rows).to_csv(
                os.path.join(pathway_output_dir, "graph_edges.csv"),
                index=False
            )

            pd.DataFrame(node_rows).to_csv(
                os.path.join(pathway_output_dir, "graph_nodes.csv"),
                index=False
            )

            simple_paths: List[List[str]] = []

            for start_node in starting_nodes:
                for end_node in ending_nodes:
                    if start_node != end_node and start_node in G.nodes and end_node in G.nodes:
                        paths = list(nx.all_simple_paths(G, source=start_node, target=end_node))
                        simple_paths.extend(paths)

            paths_to_process = PathwayGraphFunctions.remove_duplicates_order(simple_paths)
            filtered_paths = PathwayGraphFunctions.find_non_subsequences(paths_to_process)

            if filtered_paths:
                ordered_paths_rows = []

                for path_id, path in enumerate(filtered_paths):
                    for step, reaction in enumerate(path):
                        ordered_paths_rows.append({
                            "pathway_id": hrsa,
                            "path_id": path_id,
                            "step": step,
                            "reaction": reaction
                        })

                pd.DataFrame(ordered_paths_rows).to_csv(
                    os.path.join(pathway_output_dir, "ordered_paths.csv"),
                    index=False
                )

    def collect_candidate_reactions(
        self,
        reactome_pathways: List[Dict[str, Any]],
        skipped_reactions: List[Dict[str, str]]
    ) -> Tuple[List[Dict[int, Any]], List[str]]:
        """
        Retrieve candidate reactions from the lowest-level Reactome pathways.

        Returns:
            pathways_with_reactions:
                Pathway hierarchy dictionaries with an added 'reactions' field.
            pathways_to_not_analyze:
                Lowest-level pathway names used to avoid treating pathways as reactions.
        """
        hierarchy_to_low_level = ReactomeAnalysisFunctions.top_level_pathways(reactome_pathways)

        # ------------------------------ Filter steps ----------------------------- #

        # Lowest-level pathway names are excluded later to avoid treating pathway
        # containers as actual reactions.

        pathways_to_not_analyze: List[str] = []
        for dic in hierarchy_to_low_level:
            last_pathway = list(dic.keys())[-1]
            pathways_to_not_analyze.append(dic[last_pathway]['name'])

        # Retrieve and keep only true reaction events from each lowest-level pathway.

        pathways_with_reactions: List[Dict[int, Any]] = []
        for dic in hierarchy_to_low_level:
            last_pathway = list(dic.keys())[-1]
            std_id = dic[last_pathway]['id']
            reactions = ReactomeAnalysisFunctions.safe_reactome_call(
                rc.content.pathway_contained_event,
                id=std_id,
                default=[]
            )
            filtered_reactions: List[Dict[str, Any]] = []
            for reaction in reactions:

                if not isinstance(reaction, int) and reaction['displayName'] not in pathways_to_not_analyze:
                    filtered_reactions.append(reaction)
                elif isinstance(reaction, int):
                    rhsa = "R-HSA-" + str(reaction)
                    reaction = ReactomeAnalysisFunctions.safe_reactome_call(
                        rc.content.query_id,
                        id=rhsa,
                        enhanced=False,
                        default=None
                    )
                    if reaction is None:
                        print(f"[WARNING] Skipping reaction {rhsa}: Reactome returned no usable data.")
                        skipped_reactions.append({
                            "uniprot_ac": self.uniprot_ac,
                            "reaction_id": rhsa,
                            "pathway_id": std_id,
                            "pathway_name": dic[last_pathway].get("name", ""),
                            "reason": "query_id_returned_none_after_retries"
                        })

                        continue
                    if reaction.get('displayName') not in pathways_to_not_analyze and reaction.get('className') == "Reaction":
                        filtered_reactions.append(reaction)
            dic['reactions'] = filtered_reactions
            pathways_with_reactions.append(dic)
        return pathways_with_reactions, pathways_to_not_analyze

    def is_target_reaction(
        self,
        reaction: Dict[str, Any],
        target_info: Dict[str, Any]
    ) -> bool:
        """
        Return True if the reaction contains the target protein, one of its
        Reactome protein forms, or one of its associated complexes.
        """
        rsa = reaction["stId"]

        reaction_cache_key = (self.uniprot_ac, rsa)

        if reaction_cache_key in ReactomeAnalysisFunctions.reaction_target_cache:
            return ReactomeAnalysisFunctions.reaction_target_cache[
                reaction_cache_key
            ]

        is_target_reaction = False

        # Retrieve physical entities participating in the reaction.
        components = ReactomeAnalysisFunctions.safe_reactome_call(
            rc.content.participants_physical_entities,
            rsa,
            default=[]
        )

        components_name: List[str] = []
        for component in components:
            if component and component.get("stId"):
                components_name.append(component["stId"])

        # Detect reactions containing any complex associated with the target protein.
        if any(
            target in component_id
            for component_id in components_name
            for target in target_info["target_protein_complexes_names"]
        ):
            is_target_reaction = True

        # Detect reactions containing any Reactome form of the target protein.
        reaction_proteins = ReactomeAnalysisFunctions.safe_model_from_reactome_cached(rsa)

        if reaction_proteins is not None:
            protein_stid: List[str] = []

            for protein in reaction_proteins.get_objects_by_type(pybiopax.biopax.Protein):
                for protein_ref in getattr(protein, "xref", []):
                    protein_ref_id = getattr(protein_ref, "id", "")

                    if protein_ref_id.startswith("R-HSA"):
                        protein_stid.append(protein_ref_id)

            if any(
                target in protein_id
                for protein_id in protein_stid
                for target in target_info["all_target_protein_stId"]
            ):
                is_target_reaction = True

            if ReactomeAnalysisFunctions.biopax_model_contains_uniprot(
                reaction_proteins,
                self.uniprot_ac
            ):
                is_target_reaction = True

        ReactomeAnalysisFunctions.reaction_target_cache[
            reaction_cache_key
        ] = is_target_reaction

        return is_target_reaction



    def extract_target_reaction_information(
        self,
        pathways_with_reactions: List[Dict[int, Any]],
        pathways_to_not_analyze: List[str],
        target_info: Dict[str, Any],
        skipped_reactions: List[Dict[str, str]]
    ) -> List[Dict[Any, Any]]:
        """
        Filter candidate reactions for the target protein and extract BioPAX-level
        reaction information.
        """
        pathways_with_reactions_information: List[Dict[Any, Any]] = []
        for pathway in pathways_with_reactions:
            pathways_with_reaction_information: Dict[Any, Any] = {}
            pathways_keys = [key for key in pathway.keys() if key != 'reactions']
            reactions_to_analyze: List[Dict[str, Any]] = []
            unique_reactions = {
                reaction["stId"]: reaction
                for reaction in pathway["reactions"]
                if reaction and reaction.get("stId")
            }

            reactions = list(unique_reactions.values())

            for reaction in reactions:
                if self.is_target_reaction(
                    reaction,
                    target_info
                ):
                    reactions_to_analyze.append(reaction)

            # Extract BioPAX-level annotations for target reactions.

            reactions_information: List[Dict[str, Any]] = []

            for reaction in reactions_to_analyze:
                reaction_stID = reaction["stId"]

                reaction_info_cache_key = (
                    reaction_stID,
                    tuple(sorted(pathways_to_not_analyze))
                )

                if reaction_info_cache_key in ReactomeAnalysisFunctions.reaction_information_cache:
                    reactions_information.extend(
                        ReactomeAnalysisFunctions.reaction_information_cache[reaction_info_cache_key]
                    )
                    continue

                reaction_biopax_model = ReactomeAnalysisFunctions.safe_model_from_reactome_cached(
                    reaction_stID
                )

                if reaction_biopax_model is None:
                    print(f"[WARNING] Skipping reaction {reaction_stID}: BioPAX unavailable.")
                    skipped_reactions.append({
                        "uniprot_ac": self.uniprot_ac,
                        "reaction_id": reaction_stID,
                        "pathway_id": pathway[pathways_keys[-1]].get("id", ""),
                        "pathway_name": pathway[pathways_keys[-1]].get("name", ""),
                        "reason": "biopax_unavailable_after_retries"
                    })
                    ReactomeAnalysisFunctions.reaction_information_cache[reaction_info_cache_key] = []
                    continue

                reaction_information_list = ReactomeAnalysisFunctions.build_reaction_information(
                    reaction,
                    reaction_biopax_model,
                    pathways_to_not_analyze
                )

                ReactomeAnalysisFunctions.reaction_information_cache[reaction_info_cache_key] = reaction_information_list
                reactions_information.extend(reaction_information_list)
            for key in pathways_keys:
                pathways_with_reaction_information[key] = pathway[key]
            pathways_with_reaction_information["reaction"] = reactions_information
            pathways_with_reactions_information.append(pathways_with_reaction_information)
        return pathways_with_reactions_information

    def remove_duplicate_pathway_reactions(
        self,
        result_df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Remove duplicated reaction rows caused by overlapping Reactome pathway
        hierarchies.

        When the same reaction appears in multiple nested pathways, the row
        associated with the shorter, less specific pathway path is removed.
        """
        grouped = result_df.groupby("reaction_name")
        paths_to_remove: List[List[str]] = []
        for reaction in grouped.groups.keys():
            reaction_df = grouped.get_group(reaction)
            # list of R-HSA id of the lowest-level pathway
            lowest_pathway_id = list(set(reaction_df["lowest_pathway_id"].to_list()))
            # list of R-HSA id of reaction
            reaction_id = list(set(reaction_df["reaction_id"].to_list()))
            # find the duplicated reactions
            reaction_dict: Dict[str, Set[str]] = defaultdict(set)
            if len(lowest_pathway_id) != len(reaction_id):
                # create a dictionary with the reaction name and the R-HSA id of 
                # the lowest pathway that contains the same reaction
                for q in zip(reaction_df["lowest_pathway_id"].to_list(), reaction_df["reaction_id"].to_list()):
                    reaction_dict[q[1]].add(q[0])
                reaction_with_multiple_pathways = {k: list(v) for k, v in reaction_dict.items() if len(v) > 1}

                # for each double reaction and R-HSA id of the lowest pathway to compare:
                # extract the df and build a string with all the R-HSA ids of the pathways from that reaction,
                # put the two strings to be compared in a list 

                for reaction_id_val, pathways in reaction_with_multiple_pathways.items():
                    paths_to_check: List[str] = []
                    for pathway in pathways:
                        df_for_comparison = result_df[(result_df["reaction_name"] == reaction) & (result_df["reaction_id"] == reaction_id_val) & (result_df["lowest_pathway_id"] == pathway)]
                        df_for_comparison = df_for_comparison.astype(str)
                        pathway_cols = [col for col in result_df.columns if "pathway" in col and "id" in col]
                        pathway_cols.append("reaction_id")
                        for _, row in df_for_comparison.iterrows():
                            combined_pathway = "_".join(row[pathway_cols].dropna())
                            paths_to_check.append(combined_pathway)
                    paths_to_remove.append(list(set(paths_to_check)))

        # Identify the shorter duplicated pathway path and mark its
        # lowest_pathway_id/reaction_id pair for removal.

        lines_to_remove: List[Tuple[str, str]] = []
        for paths in paths_to_remove:
            compare_pathways = [path.split("_") for path in paths]
            processed_pathways: List[List[str]] = []
            for pathway in compare_pathways:
                pathway_no_nan: List[str] = []
                for reaction_elem in pathway:
                    if reaction_elem != "nan":
                        pathway_no_nan.append(reaction_elem)
                processed_pathways.append(pathway_no_nan)
            to_remove = min(processed_pathways, key=len)
            lines_to_remove.append((to_remove[-2], to_remove[-1]))

        result_df_filtered = result_df[~result_df.apply(lambda row: (row["lowest_pathway_id"], row["reaction_id"]) in set(lines_to_remove), axis=1)]

        return result_df_filtered

    def order_reactions_by_pathway_files(
        self,
        result_df_filtered: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Reorder reactions using pathway ordering files when available.

        If no ordering file exists for a pathway, rows are kept in their original
        order and marked as unordered.
        """
        ordered_result_df_filtered = pd.DataFrame()

        for stid in set(result_df_filtered["lowest_pathway_id"].dropna().tolist()):

            ordered_paths_file = self.output_path(
                "pathways_order",
                stid,
                "ordered_paths.csv"
            )

            if os.path.exists(ordered_paths_file):
                ordered_paths_df = pd.read_csv(ordered_paths_file)

                df_sorted = DataProcessingFunctions.reorder_dataframe_from_ordered_paths(
                    result_df_filtered,
                    stid,
                    ordered_paths_df
                )

            else:
                df_sorted = result_df_filtered[
                    result_df_filtered["lowest_pathway_id"] == stid
                ].copy()
                df_sorted["ordered"] = False

            ordered_result_df_filtered = pd.concat(
                [ordered_result_df_filtered, df_sorted],
                ignore_index=True
            )
        return ordered_result_df_filtered


    def build_and_write_result(
        self,
        pathways_with_reactions_information: List[Dict[Any, Any]],
        target_info: Dict[str, Any],
        skipped_reactions: List[Dict[str, str]]
    ) -> bool:
        """
        Build the final result DataFrame, clean duplicated pathway rows,
        order reactions when pathway ordering is available, and write outputs.
        """
        if skipped_reactions:
            pd.DataFrame(skipped_reactions).drop_duplicates().to_csv(
                self.output_path("skipped_reactions.csv"),
                index=False
            )
        result_df = DataProcessingFunctions.process_data(pathways_with_reactions_information)
        if result_df.empty:
            print(f"[INFO] No Reactome reactions found for {self.uniprot_ac}.")
            return False

        sorted_columns = sorted(result_df.columns, key=DataProcessingFunctions.column_sort_key)
        result_df = result_df[sorted_columns]

        # Propagate complex and stoichiometry annotations within each reaction group.
        if "complex_of" in result_df.columns:
            result_df = result_df.groupby("reaction_name", group_keys=False).apply(DataProcessingFunctions.fill_complex_info)

        # Remove rows representing protein families rather than individual proteins.
        result_df = result_df.loc[result_df['is_a_protein_family'] == False]
        if result_df.empty:
            print(
                f"[INFO] No Reactome reactions left for {self.uniprot_ac} "
                "after removing protein families."
            )
            return False
        
        result_df_filtered = self.remove_duplicate_pathway_reactions(result_df)
        
        ordered_result_df_filtered = self.order_reactions_by_pathway_files(
            result_df_filtered
        )

        if ordered_result_df_filtered.empty:
            print(f"[INFO] Final Reactome output is empty for {self.uniprot_ac}.")
            return False

        ordered_result_df_filtered.insert(
            0,
            "target_uniprot_ac",
            self.uniprot_ac
        )

        ordered_result_df_filtered.insert(
            1,
            "target_name",
            target_info["target_name"]
        )

        ordered_result_df_filtered.to_csv(
            self.output_path("result.csv"),
            sep=",",
            index=False
        )

        return True

    def run(self) -> bool:
        """Execute the full Reactome analysis workflow."""
        if rc is None or pybiopax is None:
            raise ImportError(
                "reactome2py and pybiopax must be installed to run this script"
            )

        skipped_reactions: List[Dict[str, str]] = []

        reactome_pathways = self.retrieve_reactome_pathways()

        if not self.skip_pathway_order:
            self.write_pathway_ordering(reactome_pathways)

        target_info = self.resolve_target_information()

        pathways_with_reactions, pathways_to_not_analyze = (
            self.collect_candidate_reactions(
                reactome_pathways,
                skipped_reactions
            )
        )

        pathways_with_reactions_information = (
            self.extract_target_reaction_information(
                pathways_with_reactions,
                pathways_to_not_analyze,
                target_info,
                skipped_reactions
            )
        )

        return self.build_and_write_result(
            pathways_with_reactions_information,
            target_info,
            skipped_reactions
        )
           