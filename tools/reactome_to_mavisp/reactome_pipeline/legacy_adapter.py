from __future__ import annotations

from collections import defaultdict
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple


class LegacyEventAdapter:
    """
    Convert the generic local EventContext representation into the nested
    structure expected by the legacy DataProcessingFunctions.process_data().

    The adapter never queries Reactome. Everything is reconstructed from the
    EventContext entity registry, pathway chains, event roles and controls.
    """

    @staticmethod
    def _unique(values: Iterable[Any]) -> List[Any]:
        result: List[Any] = []
        seen: Set[Any] = set()
        for value in values:
            if value in seen:
                continue
            seen.add(value)
            result.append(value)
        return result

    @staticmethod
    def _entity_id(entity: Dict[str, Any]) -> Optional[str]:
        """Prefer a Reactome stable ID; otherwise use the BioPAX UID."""
        for xref in entity.get("xref", []) or []:
            db = str(xref.get("db", "") or "")
            xid = str(xref.get("id", "") or "")
            if db.lower() == "reactome" and xid.startswith("R-HSA-"):
                return xid
        return entity.get("uid")

    @staticmethod
    def _entity_label(entity: Dict[str, Any]) -> str:
        entity_type = entity.get("entity_type") or "PhysicalEntity"
        name = entity.get("display_name") or entity.get("uid") or ""
        return f"{entity_type}({name})"

    @staticmethod
    def _location_value(entity: Dict[str, Any]) -> Any:
        location = entity.get("cellular_location")
        if not location:
            return []
        terms = location.get("terms", []) or []
        if len(terms) == 1:
            return terms[0]
        return terms

    @staticmethod
    def _legacy_feature_location(location: Optional[Dict[str, Any]]) -> Any:
        """Convert EventContext sequence locations to legacy-parseable strings."""
        if not location:
            return []

        location_type = location.get("type")

        if location_type == "SequenceSite":
            position = location.get("position")
            if position is None:
                return "SequenceSite"
            return f"SequenceSite({position})"

        if location_type == "SequenceInterval":
            begin = location.get("begin") or {}
            end = location.get("end") or {}
            positions = ""
            if begin.get("position") is not None:
                positions += f"({begin['position']})"
            if end.get("position") is not None:
                positions += f"({end['position']})"
            return f"SequenceInterval{positions}"

        return str(location)

    @staticmethod
    def _legacy_modification_type(values: List[str]) -> Any:
        # data_processing.py historically searches for text in double quotes.
        if not values:
            return []
        if len(values) == 1:
            return f'"{values[0]}"'
        return "_".join(f'"{value}"' for value in values)

    def _legacy_features(self, entity: Dict[str, Any]) -> List[Dict[str, Any]]:
        features: List[Dict[str, Any]] = []
        for feature in entity.get("features", []) or []:
            record: Dict[str, Any] = {
                "feature_location": self._legacy_feature_location(feature.get("location")),
            }
            modification_type = feature.get("modification_type", []) or []
            if modification_type:
                record["modification_type"] = self._legacy_modification_type(modification_type)
            features.append(record)
        return features

    @staticmethod
    def _parent_maps(
        registry: Dict[str, Dict[str, Any]],
    ) -> Tuple[
        Dict[str, List[Tuple[str, Dict[str, Any]]]],
        Dict[str, List[Tuple[str, Dict[str, Any]]]],
    ]:
        component_parents = defaultdict(list)
        member_parents = defaultdict(list)

        for parent_uid, parent in registry.items():
            for child in parent.get("children", []) or []:
                child_uid = child.get("entity_uid")
                if not child_uid:
                    continue
                if child.get("relation") == "component":
                    component_parents[child_uid].append((parent_uid, child))
                elif child.get("relation") == "member_physical_entity":
                    member_parents[child_uid].append((parent_uid, child))

        return component_parents, member_parents

    def _legacy_proteins(self, context: Dict[str, Any]) -> List[Dict[str, Any]]:
        registry = context["entities"]
        component_parents, member_parents = self._parent_maps(registry)
        proteins: List[Dict[str, Any]] = []

        for uid, entity in registry.items():
            if entity.get("entity_type") != "Protein":
                continue

            component_of = []
            for parent_uid, edge in component_parents.get(uid, []):
                parent = registry[parent_uid]
                component_of.append(
                    {
                        "stId": self._entity_id(parent),
                        "display_name": parent.get("display_name"),
                        "stoichiometry": edge.get("stoichiometry"),
                    }
                )

            member_of = [
                self._entity_label(registry[parent_uid])
                for parent_uid, _ in member_parents.get(uid, [])
            ]

            member_children = []
            for child in entity.get("children", []) or []:
                if child.get("relation") != "member_physical_entity":
                    continue
                child_uid = child.get("entity_uid")
                if child_uid in registry:
                    member_children.append(self._entity_label(registry[child_uid]))

            proteins.append(
                {
                    "display_name": entity.get("display_name"),
                    "uniprot_ac": entity.get("uniprot_accessions", []),
                    "cellular_location": self._location_value(entity),
                    "stId": self._entity_id(entity),
                    "feature": self._legacy_features(entity),
                    "member_physical_entity_of": member_of,
                    "member_physical_entity": member_children,
                    "component_of": component_of,
                }
            )

        return proteins

    def _descendant_complex_ids(
        self,
        uid: str,
        registry: Dict[str, Dict[str, Any]],
        visited: Optional[Set[str]] = None,
    ) -> List[str]:
        if visited is None:
            visited = set()
        if uid in visited:
            return []

        visited = set(visited)
        visited.add(uid)
        entity = registry[uid]
        result: List[str] = []

        if entity.get("entity_type") == "Complex":
            entity_id = self._entity_id(entity)
            if entity_id:
                result.append(entity_id)

        for child in entity.get("children", []) or []:
            child_uid = child.get("entity_uid")
            if child_uid not in registry:
                continue
            if registry[child_uid].get("entity_type") == "Complex":
                result.extend(self._descendant_complex_ids(child_uid, registry, visited))

        return self._unique(result)

    def _direct_complex_components(
        self,
        complex_entity: Dict[str, Any],
        registry: Dict[str, Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        result: List[Dict[str, Any]] = []
        for child in complex_entity.get("children", []) or []:
            child_uid = child.get("entity_uid")
            if child_uid not in registry:
                continue
            child_entity = registry[child_uid]
            result.append(
                {
                    "stId": self._entity_id(child_entity),
                    "display_name": child_entity.get("display_name"),
                    "entity_type": child_entity.get("entity_type"),
                    "relation": child.get("relation"),
                    "stoc_coefficient": child.get("stoichiometry"),
                }
            )
        return result

    def _legacy_complexes(self, context: Dict[str, Any]) -> List[Dict[str, Any]]:
        registry = context["entities"]
        result = []

        for uid, entity in registry.items():
            if entity.get("entity_type") != "Complex":
                continue
            result.append(
                {
                    "stId": self._entity_id(entity),
                    "complex_name": entity.get("display_name"),
                    "complexes_ids": self._descendant_complex_ids(uid, registry),
                    "complex_components": self._direct_complex_components(entity, registry),
                }
            )

        return result

    def _role_labels(self, context: Dict[str, Any], role: str) -> List[str]:
        registry = context["entities"]
        labels = []
        for participant in context.get("roles", {}).get(role, []):
            uid = participant.get("entity_uid")
            if uid in registry:
                labels.append(self._entity_label(registry[uid]))
        return labels

    def _legacy_biochemical(self, context: Dict[str, Any]) -> List[Dict[str, Any]]:
        roles = context.get("roles", {})
        event_specific = context.get("event_specific", {})
        record: Dict[str, Any] = {}

        if "LEFT" in roles or "RIGHT" in roles:
            record["Left"] = self._role_labels(context, "LEFT")
            record["Right"] = self._role_labels(context, "RIGHT")
            record["Conversion_Direction"] = event_specific.get("conversion_direction")
            ec_numbers = event_specific.get("e_c_number", [])

            if ec_numbers:
                record["EC_Number"] = ec_numbers
        else:
            # Do not invent left/right semantics for TemplateReaction or generic Interaction.
            for role in ("TEMPLATE", "PRODUCT", "PARTICIPANT"):
                if role in roles:
                    record[role.title()] = self._role_labels(context, role)

        return [record] if record else []

    def _legacy_controls(self, context: Dict[str, Any]) -> List[Dict[str, Any]]:
        registry = context["entities"]
        result = []

        for control in context.get("controls", []) or []:
            controllers = []
            controller_ids = []
            for uid in control.get("controller_uids", []) or []:
                entity = registry.get(uid)
                if not entity:
                    continue
                controllers.append(self._entity_label(entity))
                controller_ids.append(self._entity_id(entity))

            result.append(
                {
                    "controller": controllers,
                    "controller_stid": controller_ids,
                    "control_type": control.get("control_type"),
                    "control_class": control.get("control_class"),
                }
            )

        return result

    @staticmethod
    def _legacy_catalytic(context: Dict[str, Any]) -> Dict[str, Any]:
        event_specific = context.get("event_specific", {})
        result: Dict[str, Any] = {}

        for key in ("delta_g", "delta_h", "delta_s", "k_eq"):
            value = event_specific.get(key)
            if value not in (None, [], ""):
                result[key] = value

        return result

    @staticmethod
    def _join_values(values: Iterable[Any]) -> str:
        """
        Join non-empty values while preserving order and removing duplicates.
        """
        seen = set()
        result = []

        for value in values or []:
            if value in (None, ""):
                continue

            value = str(value)

            if value in seen:
                continue

            seen.add(value)
            result.append(value)

        return "_".join(result)

    @classmethod
    def _collect_variant_values(
        cls,
        variants: Iterable[Dict[str, Any]],
        key: str,
    ) -> List[Any]:
        """
        Collect one field across disease-variant records.

        Scalar values are kept as scalars; list-valued fields are flattened.
        Order is preserved and duplicates are removed.
        """
        values: List[Any] = []

        for variant in variants or []:
            value = variant.get(key)

            if value in (None, "", []):
                continue

            if isinstance(value, (list, tuple, set)):
                values.extend(value)
            else:
                values.append(value)

        return cls._unique(values)

    

    @classmethod
    def _legacy_disease_event_metadata(
        cls,
        context: Dict[str, Any],
    ) -> Dict[str, Any]:
        """
        Build reaction-level disease metadata.

        These values describe the disease event itself and are therefore kept
        once on the reaction. Mutation-specific metadata remain in the raw
        ``disease_variants`` list and are matched to the corresponding protein
        later by DataProcessingFunctions using the Reactome stable entity ID.
        """
        variants = context.get("disease_variants", []) or []

        result: Dict[str, Any] = {
            "event_context_type": context.get(
                "context_type",
                "NORMAL",
            ),
            "query_uniprot_ac": context.get(
                "query_uniprot_ac"
            ),
            "disease_cross_reference": "",
            "disease_identifier": "",
            "disease_reaction_id": "",
            "disease_reaction_name": "",
            "disease_pathway_id": "",
            "disease_pathway_name": "",
            "normal_reaction_id": "",
            "normal_reaction_name": "",
            "normal_pathway_id": "",
            "normal_pathway_name": "",
            "disease_variants": variants,
        }

        if not variants:
            return result

        field_map = {
            "disease_cross_reference": "cross_references",
            "disease_identifier": "disease_identifiers",
            "disease_reaction_id": "disease_reaction_ids",
            "disease_reaction_name": "disease_reaction_names",
            "disease_pathway_id": "disease_pathway_ids",
            "disease_pathway_name": "disease_pathway_names",
            "normal_reaction_id": "normal_reaction_ids",
            "normal_reaction_name": "normal_reaction_names",
            "normal_pathway_id": "normal_pathway_ids",
            "normal_pathway_name": "normal_pathway_names",
        }

        for output_key, variant_key in field_map.items():
            result[output_key] = cls._join_values(
                cls._collect_variant_values(
                    variants,
                    variant_key,
                )
            )

        return result

    def _legacy_reaction(
        self,
        context: Dict[str, Any],
    ) -> Dict[str, Any]:
        """
        Convert one EventContext into exactly one legacy reaction record.

        A DISEASE EventContext is *not* duplicated once per variant. Instead,
        all disease variants remain attached to this single reaction and are
        later assigned to their matching mutant protein by stable entity ID.
        """
        reaction = {
            "stID": context.get("stable_id"),
            "Display_Name": context.get("display_name"),
            "proteins": self._legacy_proteins(context),
            "complex": self._legacy_complexes(context),
            "biochemical": self._legacy_biochemical(context),
            "control_information": self._legacy_controls(context),
            "catalytic": self._legacy_catalytic(context),
        }

        reaction.update(
            self._legacy_disease_event_metadata(context)
        )

        return reaction

    def adapt(
        self,
        context: Dict[str, Any],
    ) -> List[Dict[Any, Any]]:
        """
        Emit one legacy entry for each pathway chain of one EventContext.

        Each entry contains one reaction record. Disease variants are stored
        inside that reaction and are not expanded here.
        """
        reaction = self._legacy_reaction(context)
        chains = context.get("pathway_chains", []) or []

        if not chains:
            return [{"reaction": [reaction]}]

        result: List[Dict[Any, Any]] = []

        for chain in chains:
            entry: Dict[Any, Any] = {}

            for index, pathway in enumerate(chain):
                entry[index] = {
                    "name": pathway.get("display_name"),
                    "id": pathway.get("stable_id"),
                }

            entry["reaction"] = [reaction]
            result.append(entry)

        return result

    def adapt_many(
        self,
        contexts: Iterable[Dict[str, Any]],
    ) -> List[Dict[Any, Any]]:
        result: List[Dict[Any, Any]] = []

        for context in contexts:
            result.extend(self.adapt(context))

        return result
