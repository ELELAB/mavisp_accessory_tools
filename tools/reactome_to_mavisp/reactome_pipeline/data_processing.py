import re
from typing import Any, Dict, List, Set, Tuple

import numpy as np
import pandas as pd


class DataProcessingFunctions:
    """Functions for cleaning and manipulating DataFrames."""

    @staticmethod
    def replace_empty_with_nan(value: Any) -> Any:
        if isinstance(value, (set, list, dict)) and len(value) == 0:
            return np.nan
        return value

    @staticmethod
    def _complex_record_id(complex_record: Dict[str, Any]) -> Any:
        """
        Return the stable ID of a legacy complex record.

        LegacyEventAdapter._descendant_complex_ids() always places the
        current complex itself first in ``complexes_ids``, followed by nested
        descendant complex IDs. This lets data processing reconstruct nested
        complex membership locally without querying Reactome.
        """
        explicit_id = complex_record.get("stId")
        if explicit_id:
            return explicit_id

        complex_ids = complex_record.get("complexes_ids", []) or []
        return complex_ids[0] if complex_ids else None

    @staticmethod
    def _build_complex_index(
        reaction_complexes: List[Dict[str, Any]],
    ) -> Dict[str, Dict[str, Any]]:
        """Index local legacy complex records by their own Reactome stable ID."""
        result: Dict[str, Dict[str, Any]] = {}

        for complex_record in reaction_complexes:
            complex_id = DataProcessingFunctions._complex_record_id(
                complex_record
            )
            if complex_id:
                result[str(complex_id)] = complex_record

        return result

    @staticmethod
    def _complex_contains_entity(
        complex_id: Any,
        target_stid: Any,
        complex_index: Dict[str, Dict[str, Any]],
        visited: Any = None,
    ) -> bool:
        """
        Return True when ``target_stid`` occurs below ``complex_id``.

        Only locally serialized Complex -> component/member relationships are
        traversed. Stoichiometry is not multiplied across nested levels,
        preserving the historical output semantics.
        """
        if not complex_id or not target_stid:
            return False

        complex_id = str(complex_id)
        target_stid = str(target_stid)

        if visited is None:
            visited = set()

        if complex_id in visited:
            return False

        visited = set(visited)
        visited.add(complex_id)

        complex_record = complex_index.get(complex_id)
        if not complex_record:
            return False

        for entity in complex_record.get("complex_components", []) or []:
            entity_id = entity.get("stId")

            if entity_id is None:
                continue

            entity_id = str(entity_id)

            if entity_id == target_stid:
                return True

            if entity_id in complex_index:
                if DataProcessingFunctions._complex_contains_entity(
                    entity_id,
                    target_stid,
                    complex_index,
                    visited,
                ):
                    return True

        return False

    @staticmethod
    def _flatten_unique(values: List[Any]) -> List[Any]:
        """Flatten nested list-like values, preserving order and uniqueness."""
        result: List[Any] = []
        seen = set()

        for value in values:
            items = (
                value
                if isinstance(value, (list, tuple, set))
                else [value]
            )

            for item in items:
                if item in (None, ""):
                    continue

                key = str(item)
                if key in seen:
                    continue

                seen.add(key)
                result.append(item)

        return result

    @staticmethod
    def _join_variant_field(
        variants: List[Dict[str, Any]],
        key: str,
    ) -> str:
        """Join one disease-variant field across matched variant records."""
        values = [
            variant.get(key)
            for variant in variants
            if variant.get(key) not in (None, "", [])
        ]

        flattened = DataProcessingFunctions._flatten_unique(values)

        return "_".join(str(value) for value in flattened)

    @staticmethod
    def _matched_disease_variants(
        reaction: Dict[str, Any],
        protein_stid: Any,
    ) -> List[Dict[str, Any]]:
        """
        Return only disease variants represented by this protein row.

        Matching uses the Reactome stable ID of the mutant physical entity.
        This prevents mutation metadata from being copied onto unrelated
        proteins participating in the same disease reaction.
        """
        if protein_stid in (None, ""):
            return []

        protein_stid = str(protein_stid)

        return [
            variant
            for variant in reaction.get(
                "disease_variants",
                [],
            ) or []
            if str(
                variant.get(
                    "variant_entity_stable_id",
                    "",
                )
            ) == protein_stid
        ]

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

                    # ``protein["component_of"]`` is populated by LegacyEventAdapter
                    # only from direct BioPAX ``component`` relationships.  Therefore
                    # it is the authoritative source for the actual complexes that
                    # directly contain this protein.  Do not recursively promote
                    # outer EntitySets/containers into ``complex_of``.
                    direct_complexes = protein.get("component_of", []) or []

                    # Complex records are still useful for recovering the distinct
                    # BioPAX relation in which a direct complex is a
                    # ``member_physical_entity`` of an EntitySet-like container.
                    reaction_complexes = reaction.get("complex", []) or []

                    # Index complex records by their own stable ID.
                    complex_index = DataProcessingFunctions._build_complex_index(
                        reaction_complexes
                    )

                    # Build child-complex -> EntitySet-like parent(s) mapping from
                    # explicit BioPAX member_physical_entity edges.  This relation is
                    # deliberately kept separate from component_of.
                    entity_set_parents: Dict[str, List[str]] = {}
                    for parent_complex in reaction_complexes:
                        parent_name = parent_complex.get("complex_name")
                        for child in parent_complex.get("complex_components", []) or []:
                            if child.get("relation") != "member_physical_entity":
                                continue

                            child_id = child.get("stId")
                            if child_id in (None, "") or parent_name in (None, ""):
                                continue

                            child_id = str(child_id)
                            entity_set_parents.setdefault(child_id, [])
                            if str(parent_name) not in entity_set_parents[child_id]:
                                entity_set_parents[child_id].append(str(parent_name))

                    complex_names: List[str] = []
                    complex_ids: List[str] = []
                    stoichiometry: List[str] = []
                    complex_entity_sets: List[str] = []

                    # Keep one entry per direct Complex stable ID.  Two different
                    # complexes may legitimately have the same display name, so the
                    # stable-ID column is retained to keep them distinguishable.
                    seen_complex_ids = set()

                    for direct_complex in direct_complexes:
                        complex_id = direct_complex.get("stId")
                        if complex_id in (None, ""):
                            continue

                        complex_id = str(complex_id)
                        if complex_id in seen_complex_ids:
                            continue
                        seen_complex_ids.add(complex_id)

                        complex_name = direct_complex.get("display_name")
                        if complex_name in (None, ""):
                            # Fallback to the serialized complex record if needed.
                            complex_record = complex_index.get(complex_id, {})
                            complex_name = complex_record.get("complex_name")

                        complex_names.append(
                            str(complex_name) if complex_name not in (None, "") else complex_id
                        )
                        complex_ids.append(complex_id)

                        coefficient = direct_complex.get("stoichiometry")
                        if coefficient is None:
                            stoichiometry.append("NA")
                        else:
                            stoichiometry.append(str(coefficient))

                        for entity_set_name in entity_set_parents.get(complex_id, []):
                            if entity_set_name not in complex_entity_sets:
                                complex_entity_sets.append(entity_set_name)

                    if complex_names:
                        row["complex_of"] = "_".join(complex_names)
                        row["complex_of_stid"] = "_".join(complex_ids)
                        row["stoichiometry"] = "_".join(stoichiometry)

                    if complex_entity_sets:
                        row["complex_entity_set"] = "_".join(complex_entity_sets)

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
                    row["catalytic_EC_Number"] = np.nan
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

                    # Event-context metadata.
                    row["event_context_type"] = reaction.get(
                        "event_context_type",
                        "NORMAL",
                    )

                    # --------------------------------------------------
                    # NORMAL EVENTS
                    # Preserve the historical behavior unchanged.
                    # --------------------------------------------------
                    if row["event_context_type"] != "DISEASE":

                        row["disease_name"] = ""

                        for field in (
                            "disease_cross_reference",
                            "disease_identifier",
                            "disease_reaction_id",
                            "disease_reaction_name",
                            "disease_pathway_id",
                            "disease_pathway_name",
                            "normal_reaction_id",
                            "normal_reaction_name",
                            "normal_pathway_id",
                            "normal_pathway_name",
                        ):
                            row[field] = ""

                        row["mutation"] = ""
                        row["variant"] = ""
                        row["variant_entity_stable_id"] = ""
                        row["variant_uniprot_ac"] = ""
                        row["variant_modification_class"] = ""
                        row["variant_modification_description"] = ""
                        row["functional_status"] = ""
                        row["variant_literature_pubmed"] = ""

                        rows.append(row)
                        continue

                    # --------------------------------------------------
                    # DISEASE EVENTS
                    # --------------------------------------------------

                    for field in (
                        "disease_reaction_id",
                        "disease_reaction_name",
                        "disease_pathway_id",
                        "disease_pathway_name",
                        "normal_reaction_id",
                        "normal_reaction_name",
                        "normal_pathway_id",
                        "normal_pathway_name",
                    ):
                        row[field] = reaction.get(
                            field,
                            "",
                        )

                    matched_variants = (
                        DataProcessingFunctions._matched_disease_variants(
                            reaction,
                            prot_stId,
                        )
                    )

                    # Protein participating in a DISEASE reaction but not
                    # corresponding to the mutant physical entity.
                    if not matched_variants:
                        row["disease_name"] = ""
                        row["disease_cross_reference"] = ""
                        row["disease_identifier"] = ""

                        row["mutation"] = ""
                        row["variant"] = ""
                        row["variant_entity_stable_id"] = ""
                        row["variant_uniprot_ac"] = ""
                        row["variant_modification_class"] = ""
                        row["variant_modification_description"] = ""
                        row["functional_status"] = ""
                        row["variant_literature_pubmed"] = ""

                        rows.append(row)
                        continue

                    # Exact mutant physical entity:
                    # one output row for each variant × disease association.
                    for variant in matched_variants:

                        diseases = (
                            variant.get("diseases", [])
                            or [""]
                        )

                        disease_identifiers = (
                            variant.get(
                                "disease_identifiers",
                                [],
                            )
                            or []
                        )

                        for disease_index, disease in enumerate(
                            diseases
                        ):
                            variant_row = row.copy()

                            variant_row["disease_name"] = disease

                            if disease_index < len(
                                disease_identifiers
                            ):
                                variant_row[
                                    "disease_identifier"
                                ] = disease_identifiers[
                                    disease_index
                                ]
                            else:
                                variant_row[
                                    "disease_identifier"
                                ] = ""

                            variant_row[
                                "disease_cross_reference"
                            ] = "_".join(
                                str(value)
                                for value in (
                                    variant.get(
                                        "cross_references",
                                        [],
                                    )
                                    or []
                                )
                            )

                            variant_row["variant"] = (
                                variant.get(
                                    "variant_display_name",
                                    "",
                                )
                            )

                            variant_display_name = str(
                                variant.get(
                                    "variant_display_name",
                                    "",
                                )
                            )

                            mutation_matches = re.findall(
                                r"\b[A-Z]\d+(?:[A-Z](?:fs\*\d+)?|fs\*\d+|del|ins\d+)",
                                variant_display_name,
                            )

                            variant_row["mutation"] = ";".join(
                                mutation_matches
                            )

                            variant_row[
                                "variant_entity_stable_id"
                            ] = variant.get(
                                "variant_entity_stable_id",
                                "",
                            )

                            variant_row[
                                "variant_uniprot_ac"
                            ] = variant.get(
                                "uniprot_ac",
                                "",
                            )

                            variant_row[
                                "variant_modification_class"
                            ] = variant.get(
                                "modification_class",
                                "",
                            )

                            variant_row[
                                "variant_modification_description"
                            ] = variant.get(
                                "modification_description",
                                "",
                            )

                            variant_row[
                                "functional_status"
                            ] = variant.get(
                                "functional_status",
                                "",
                            )

                            variant_row[
                                "variant_literature_pubmed"
                            ] = "_".join(
                                str(value)
                                for value in (
                                    variant.get(
                                        "literature_pubmed",
                                        [],
                                    )
                                    or []
                                )
                            )

                            rows.append(variant_row)
        df = pd.DataFrame(rows)
        df = df.applymap(DataProcessingFunctions.replace_empty_with_nan)

        # Some legacy fields intentionally remain Python lists (for example
        # UniProt accessions or multi-valued cellular locations). Pandas cannot
        # hash lists during drop_duplicates(), so build a temporary hashable
        # representation only for duplicate detection while preserving the
        # original values in the returned DataFrame.
        def _freeze_for_hashing(value: Any) -> Any:
            if isinstance(value, list):
                return tuple(
                    _freeze_for_hashing(item)
                    for item in value
                )
            if isinstance(value, set):
                return tuple(
                    sorted(
                        _freeze_for_hashing(item)
                        for item in value
                    )
                )
            if isinstance(value, dict):
                return tuple(
                    sorted(
                        (
                            key,
                            _freeze_for_hashing(item),
                        )
                        for key, item in value.items()
                    )
                )
            return value

        hashable_df = df.applymap(_freeze_for_hashing)
        df = df.loc[~hashable_df.duplicated()].copy()

        return df

    @staticmethod
    def get_column_type(col: str) -> str:
        """Classify a DataFrame column into a semantic category."""

        # Target / query
        if col in {
            "target_uniprot_ac",
            "target_name"
        }:
            return "target"

        # Pathway hierarchy
        if col == "highest_pathway" or col == "highest_pathway_id":
            return "highest_pathway"

        if col.startswith("pathway_"):
            return "pathway"

        if col == "lowest_pathway" or col == "lowest_pathway_id":
            return "lowest_pathway"

        # Reaction
        if col == "reaction_name":
            return "reaction_name"

        if col == "reaction_id":
            return "reaction_id"

        if col.startswith("reaction_"):
            return "reaction"

        if col.startswith("Controller_of_reaction_"):
            return "reaction_control"

        if col.startswith("catalytic_"):
            return "catalytic"

        # Protein / physical entity
        if col in {
            "protein",
            "uniprot_ac",
            "cellular_location",
            "SequenceInterval",
            "SequenceSite",
            "Modification_type",
        }:
            return "protein"

        # Complex / family
        if col in {
            "complex_of",
            "complex_of_stid",
            "complex_entity_set",
            "stoichiometry",
            "is_a_protein_family",
            "member_physical_entity_of",
        }:
            return "complex"

        # Event context
        if col == "event_context_type":
            return "event_context"

        # Disease
        if col in {
            "disease_name",
            "disease_cross_reference",
            "disease_identifier",
            "disease_pathway_id",
            "disease_pathway_name",
            "disease_reaction_id",
            "disease_reaction_name",
            "normal_pathway_id",
            "normal_pathway_name",
            "normal_reaction_id",
            "normal_reaction_name",
        }:
            return "disease"

        # Variant
        if col in {
            "mutation",
            "variant",
            "variant_entity_stable_id",
            "variant_uniprot_ac",
            "variant_modification_class",
            "variant_modification_description",
            "functional_status",
            "variant_literature_pubmed",
        }:
            return "variant"

        if col == "ordered":
            return "status"

        return "other"


    @staticmethod
    def column_sort_key(col: str) -> Tuple[int, int, str]:
        """Return a stable biologically meaningful ordering for output columns."""

        explicit_order = [
            # Target / query
            "target_uniprot_ac",
            "target_name",

            # Highest pathway
            "highest_pathway",
            "highest_pathway_id",

            # Lowest pathway is positioned dynamically after intermediate pathways

            # Reaction
            "reaction_name",
            "reaction_id",
            "reaction_Left",
            "reaction_Right",
            "reaction_Conversion_Direction",

            # Catalysis / regulation
            "reaction_EC_Number",
            "catalytic_EC_Number",
            "Controller_of_reaction_ACTIVATION",
            "Controller_of_reaction_INHIBITION",

            # Protein
            "protein",
            "uniprot_ac",
            "cellular_location",
            "SequenceInterval",
            "SequenceSite",
            "Modification_type",

            # Complex / family
            "complex_of",
            "complex_of_stid",
            "complex_entity_set",
            "stoichiometry",
            "is_a_protein_family",
            "member_physical_entity_of",

            # Context
            "event_context_type",

            # Disease
            "disease_name",
            "disease_cross_reference",
            "disease_identifier",
            "disease_pathway_id",
            "disease_pathway_name",
            "disease_reaction_id",
            "disease_reaction_name",
            "normal_pathway_id",
            "normal_pathway_name",
            "normal_reaction_id",
            "normal_reaction_name",

            # Variant
            "mutation",
            "variant",
            "variant_entity_stable_id",
            "variant_uniprot_ac",
            "variant_modification_class",
            "variant_modification_description",
            "functional_status",
            "variant_literature_pubmed",

            # Final status
            "ordered",
        ]

        explicit_index = {
            name: index
            for index, name in enumerate(explicit_order)
        }

        # Intermediate pathway columns need numeric ordering:
        # pathway_1, pathway_1_id, pathway_2, pathway_2_id, ...
        pathway_match = re.fullmatch(r"pathway_(\d+)(_id)?", col)

        if pathway_match:
            pathway_number = int(pathway_match.group(1))
            is_id = pathway_match.group(2) is not None

            # Place them after highest_pathway and before lowest_pathway.
            return (
                1,
                pathway_number * 2 + (1 if is_id else 0),
                col,
            )

        if col == "lowest_pathway":
            return (2, 0, col)

        if col == "lowest_pathway_id":
            return (2, 1, col)

        if col in explicit_index:
            # +100 keeps these after pathway hierarchy.
            return (
                3,
                explicit_index[col],
                col,
            )

        # Any unexpected/new columns are kept at the end,
        # rather than disappearing or breaking the pipeline.
        return (
            99,
            0,
            col,
        )
    
    @staticmethod
    def reorder_dataframe_from_ordered_paths(
        df: pd.DataFrame,
        lowest_pathway_id: str,
        ordered_paths_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Reorder reactions according to locally generated ordered_paths.csv.

        Matching is performed using Reactome reaction stable IDs rather than
        reaction display names.
        """

        df_filtered = df[
            df["lowest_pathway_id"] == lowest_pathway_id
        ].copy()

        if df_filtered.empty or ordered_paths_df.empty:
            return df_filtered

        if (
            "reaction_id" not in df_filtered.columns
            or "reaction_id" not in ordered_paths_df.columns
        ):
            df_filtered["ordered"] = False
            return df_filtered

        reaction_set = set(
            df_filtered[
                "reaction_id"
            ].dropna().astype(str).tolist()
        )

        best_path_id = None
        best_overlap: Set[str] = set()

        for path_id, path_df in ordered_paths_df.groupby(
            "path_id"
        ):
            path_reactions = set(
                path_df[
                    "reaction_id"
                ].dropna().astype(str).tolist()
            )

            overlap = reaction_set.intersection(
                path_reactions
            )

            if len(overlap) > len(best_overlap):
                best_overlap = overlap
                best_path_id = path_id

        if best_path_id is None or not best_overlap:
            df_filtered["ordered"] = False
            return df_filtered

        best_path_df = (
            ordered_paths_df[
                ordered_paths_df[
                    "path_id"
                ]
                == best_path_id
            ]
            .sort_values("step")
        )

        ordered_reactions = [
            str(reaction_id)
            for reaction_id in best_path_df[
                "reaction_id"
            ].tolist()
            if (
                pd.notna(reaction_id)
                and str(reaction_id) in reaction_set
            )
        ]

        # Preserve deterministic order for reactions not represented
        # in the selected linear path.
        missing_reactions = sorted(
            reaction_id
            for reaction_id in reaction_set
            if reaction_id not in ordered_reactions
        )

        final_order = (
            ordered_reactions
            + missing_reactions
        )

        df_filtered["ordered"] = (
            df_filtered[
                "reaction_id"
            ]
            .astype(str)
            .isin(ordered_reactions)
        )

        df_filtered["_reaction_order"] = pd.Categorical(
            df_filtered[
                "reaction_id"
            ].astype(str),
            categories=final_order,
            ordered=True,
        )

        df_filtered = (
            df_filtered
            .sort_values(
                "_reaction_order",
                kind="stable",
            )
            .drop(
                columns=[
                    "_reaction_order"
                ]
            )
        )

        return df_filtered

    @staticmethod
    def fill_complex_info(group: pd.DataFrame) -> pd.DataFrame:
        """
        Propagate complex annotations from a protein-family/entity-set row to
        its member proteins when the member itself is not a direct complex
        component. Keep all complex-related columns aligned.
        """
        if "complex_of" not in group.columns:
            return group

        for index, row in group.iterrows():
            if pd.isna(row.get("complex_of")):
                matching_row = group[
                    group["protein"] == row.get("member_physical_entity_of")
                ]

                if matching_row.empty:
                    continue

                source = matching_row.iloc[0]

                for column in (
                    "complex_of",
                    "complex_of_stid",
                    "complex_entity_set",
                    "stoichiometry",
                ):
                    if column in group.columns:
                        group.at[index, column] = source.get(column)

        return group
