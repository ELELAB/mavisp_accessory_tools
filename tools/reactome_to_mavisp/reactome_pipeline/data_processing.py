import re
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd

try:
    import reactome2py as rc
    from reactome2py import analysis, content, utils
except ImportError:
    rc = None  # type: ignore

from reactome_pipeline.reactome_analysis import ReactomeAnalysisFunctions

class DataProcessingFunctions:
    """Functions for cleaning and manipulating DataFrames."""

    @staticmethod
    def replace_empty_with_nan(value: Any) -> Any:
        if isinstance(value, (set, list)) and len(value) == 0:
            return np.nan
        return value

    @staticmethod
    def process_data(data: List[Dict[str, Any]]) -> pd.DataFrame:
        """Flatten nested pathway/reaction annotations into a tabular DataFrame.

        Each row represents one protein entry within one Reactome reaction and includes
        pathway hierarchy, reaction metadata, protein features, complex membership,
        stoichiometry, regulatory information, and disease annotations when available.
        """
        rows: List[Dict[str, Any]] = []
        for entry in data:
            pathways_keys = [int(k) for k in entry.keys() if k != 'reaction']
            pathways_keys.sort()
            for reaction in entry['reaction']:
                proteins = reaction.get('proteins', [])
                for protein in proteins:
                    row: Dict[str, Any] = {}
                    row["protein"] = protein['display_name']
                    row["uniprot_ac"] = protein['uniprot_ac']
                    row["cellular_location"] = protein['cellular_location']
                    prot_stId = protein['stId']

                    # Protein sequence features and residue modifications

                    sequence_interval: List[str] = []
                    sequence_site: List[str] = []
                    modification_type: List[str] = []
                    if 'feature' in protein:
                        features = protein['feature']
                        for feature in features:

                            if 'modification_type' in feature:
                                mod_type_match = re.search(r'"(.*?)"', str(feature['modification_type']))
                                if mod_type_match:
                                    modification_type.append(mod_type_match.group(1))
                                elif mod_type_match is None:
                                    modification_type.append("could be altered by a mutation event")
                                else:
                                    raise ValueError("Unexpected modification_type annotation")
                            position = re.findall(r'\((\d+)\)', str(feature['feature_location']))
                            position = list(map(str, position))

                            if "SequenceInterval" in str(feature['feature_location']):
                                sequence_interval.append("-".join(position))
                            if "SequenceSite" in str(feature['feature_location']) and not "SequenceInterval" in str(feature['feature_location']):
                                sequence_site.append("-".join(position))
                    if sequence_interval:
                        row["SequenceInterval"] = "_".join(sequence_interval)
                    if sequence_site:
                        row["SequenceSite"] = "_".join(sequence_site)
                    if modification_type:
                        row["Modification_type"] = "_".join(modification_type)

                    # Protein family annotation
                    if not protein['member_physical_entity_of'] and protein['member_physical_entity']:
                        row['is_a_protein_family'] = True
                    else:
                        row['is_a_protein_family'] = False

                    # Complex membership and stoichiometry

                    # Get the complexes in which the protein participates
                    complex_of = protein['component_of']
                    # Get all the complexes of the reaction
                    reaction_complexes = reaction.get('complex', {})
                    # Annotate all the complexes in which the protein participates
                    complexes_per_protein: List[str] = []
                    # Annotate the stoichiometry of the protein in the corresponding complex
                    stoichiometry: List[str] = []
                    protein_complex_ids = {
                        prot_complex["stId"]
                        for prot_complex in complex_of
                        if prot_complex.get("stId")
                    }

                    prot_reaction_complexes = [
                        rcx for rcx in reaction_complexes
                        if protein_complex_ids.intersection(set(rcx.get("complexes_ids", [])))
                    ]
                    for prc in prot_reaction_complexes:
                        for entity in prc["complex_components"]:
                            entity_id = entity["stId"]

                            if prot_stId == entity_id:
                                complexes_per_protein.append(str(prc["complex_name"]))
                                stoichiometry.append(str(entity["stoc_coefficient"]))

                            else:
                                if entity_id not in ReactomeAnalysisFunctions.complex_entities_cache:
                                    ReactomeAnalysisFunctions.complex_entities_cache[entity_id] = (
                                        ReactomeAnalysisFunctions.safe_reactome_call(
                                            rc.content.entities_complex,
                                            entity_id,
                                            default=[]
                                        ) if rc else []
                                    )

                                components_in_complex = ReactomeAnalysisFunctions.complex_entities_cache[entity_id]

                                for component in components_in_complex:
                                    if component and prot_stId == component.get("stId"):
                                        if str(prc["complex_name"]) not in complexes_per_protein:
                                            complexes_per_protein.append(str(prc["complex_name"]))
                                            stoichiometry.append(str(entity["stoc_coefficient"]))
                    if complex_of:
                        row["complex_of"] = "_".join(complexes_per_protein)
                        row["stoichiometry"] = "_".join(stoichiometry)

                    # Parent protein family, when available
        
                    match = re.search(r"\((.*?)\)", str(protein['member_physical_entity_of']))
                    if match:
                        result = match.group(1)
                    else:
                        result = None

                    row['member_physical_entity_of'] = result

                    # Pathway hierarchy annotation

                    row['highest_pathway'] = entry[pathways_keys[0]].get('name')
                    row['highest_pathway_id'] = entry[pathways_keys[0]].get('id')
                    row['lowest_pathway'] = entry[pathways_keys[-1]].get('name')
                    row['lowest_pathway_id'] = entry[pathways_keys[-1]].get('id')
                    for key, index in zip(pathways_keys[1:-1], range(1, len(pathways_keys) - 1)):
                        row[f'pathway_{index}'] = entry[key].get('name')
                        row[f'pathway_{index}_id'] = entry[key].get('id')

                    # Reaction identifiers
                    row['reaction_id'] = reaction['stID']
                    row['reaction_name'] = reaction["Display_Name"]

                    # Biochemical reaction metadata
                    biochemical = reaction.get('biochemical', {})
                    if len(biochemical) > 1:
                        raise ValueError("Unexpected multiple biochemical entries")
                    for bio in biochemical:
                        for key, value in bio.items():
                            if isinstance(value, list):
                                str_value: List[str] = []
                                for val in value:
                                    str_value.append(str(val))
                                if str_value:
                                    row[f'reaction_{key}'] = " ".join(str_value)
                                else:
                                    row[f'reaction_{key}'] = []
                            else:
                                row[f'reaction_{key}'] = value

                    # Reaction regulation, such as activation or inhibition

                    control_information = reaction.get('control_information', [])

                    # Store the regulatory role of each controller, such as activation or inhibition.
                    
                    for reaction_control in control_information:
                        controllers: List[str] = []

                        for controller in reaction_control['controller']:
                            match = re.search(r"\((.*?)\)", str(controller))
                            if match:
                                controllers.append(match.group(1))

                        control_type = reaction_control['control_type']
                        row[f'Controller_of_reaction_{control_type}'] = "_".join(controllers)
                            
                    # Catalytic metadata, when available
                    catalytic = reaction.get('catalytic', {})
                    for key, value in catalytic.items():
                        if isinstance(value, list):
                            str_value: List[str] = []
                            for val in value:
                                str_value.append(str(val))
                            if str_value:
                                row[f'catalytic_{key}'] = " ".join(str_value)
                            else:
                                row[f'catalytic_{key}'] = []
                        else:
                            row[f'catalytic_{key}'] = value

                    # Disease annotation, when available
                    row["disease_name"] = reaction["disease"]
                    rows.append(row)
        df = pd.DataFrame(rows)
        df = df.applymap(DataProcessingFunctions.replace_empty_with_nan)
        df.drop_duplicates(inplace=True)
        return df

    @staticmethod
    def get_column_type(col: str) -> str:
        """Classify a DataFrame column into a semantic category."""

        # Pathways: Columns that start with "highest_pathway", "pathway_", or "lowest_pathway"
        if "pathway" in col:
            if "highest_pathway" in col:
                return "highest_pathway"
            elif "lowest_pathway" in col:
                return "lowest_pathway"
            elif col.startswith("pathway_"):
                return "pathway"
        if col.startswith("disease_name"):
            return "disease_name"

        # Reaction: Columns that start with "reaction"
        if col.startswith("reaction"):
            if "reaction_name" in col:
                return "reaction_name"
            elif "reaction_id" in col:
                return "reaction_id"
            elif "reaction_Left" in col:
                return "reaction_left"
            elif "reaction_Right" in col:
                return "reaction_right"
            elif "reaction_Conversion_Direction" in col:
                return "reaction_conversion_direction"
            elif "reaction_controller" in col:
                return "reaction_controller"
            return "reaction"
        # Catalytic: Columns that start with "catalytic"

        if col.startswith("catalytic"):
            return "catalytic"

        # Complex or other physical entities: Columns that start with "Complex", "SmallMolecule", or "Protein"
        if col.startswith("Complex") or col.startswith("SmallMolecule") or col.startswith("Protein"):
            return "Physical_entity"
        return "Physical_entity" # Default category for uncategorized columns

    @staticmethod
    def column_sort_key(col: str) -> Tuple[int, str]:
        """Return a key for sorting columns."""
        column_type = DataProcessingFunctions.get_column_type(col)
        # Assign priority for each column type
        type_priority = {
            "highest_pathway": 0,
            "disease_name": 1,
            "pathway": 2,
            "lowest_pathway": 3,
            "reaction_name": 4,
            "reaction_id": 5,
            "reaction_right": 7,
            "reaction_left": 6,
            "reaction_conversion_direction": 8,
            "reaction_controller": 9,
            "reaction": 10,
            "catalytic": 11,
            "complex": 12,
            "Physical_entity": 13,
        }
        # First sort by category type, then alphabetically within each category
        return (type_priority.get(column_type, 99), col)

    @staticmethod
    def reorder_dataframe_from_ordered_paths(
        df: pd.DataFrame,
        lowest_pathway_id: str,
        ordered_paths_df: pd.DataFrame
    ) -> pd.DataFrame:

        df_filtered = df[df["lowest_pathway_id"] == lowest_pathway_id].copy()

        if df_filtered.empty or ordered_paths_df.empty:
            return df_filtered

        reaction_set = set(df_filtered["reaction_name"].dropna().tolist())

        best_path_id = None
        best_overlap = set()

        for path_id, path_df in ordered_paths_df.groupby("path_id"):
            path_reactions = set(path_df["reaction"].dropna().tolist())
            overlap = reaction_set.intersection(path_reactions)

            if len(overlap) > len(best_overlap):
                best_overlap = overlap
                best_path_id = path_id

        if best_path_id is None or not best_overlap:
            df_filtered["ordered"] = False
            return df_filtered

        best_path_df = ordered_paths_df[
            ordered_paths_df["path_id"] == best_path_id
        ].sort_values("step")

        ordered_reactions = [
            r for r in best_path_df["reaction"].tolist()
            if r in reaction_set
        ]

        missing_reactions = [
            r for r in reaction_set
            if r not in ordered_reactions
        ]

        final_order = ordered_reactions + missing_reactions

        df_filtered["ordered"] = df_filtered["reaction_name"].isin(ordered_reactions)

        df_filtered["reaction_name"] = pd.Categorical(
            df_filtered["reaction_name"],
            categories=final_order,
            ordered=True
        )

        df_filtered = df_filtered.sort_values("reaction_name")

        return df_filtered

    @staticmethod
    def fill_complex_info(group: pd.DataFrame) -> pd.DataFrame:
        """Fill missing complex and stoichiometry values in grouped DataFrame."""
        for index, row in group.iterrows():
            if pd.isna(row["complex_of"]):
                matching_row = group[group["protein"] == row["member_physical_entity_of"]]
                if not matching_row.empty:
                    group.at[index, "complex_of"] = matching_row.iloc[0]["complex_of"]
                    group.at[index, "stoichiometry"] = matching_row.iloc[0]["stoichiometry"]
        return group