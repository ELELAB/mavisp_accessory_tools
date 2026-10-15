from __future__ import annotations

from collections import defaultdict
from typing import Any, Dict, List, Optional, Set


class EventContextBuilder:
    """
    Generic BioPAX Level 3 event-context builder.

    The builder follows the BioPAX class hierarchy:
    - Conversion family -> LEFT / RIGHT + participant stoichiometry
    - TemplateReaction -> TEMPLATE / PRODUCT
    - other Interaction subclasses -> PARTICIPANT
    - Control subclasses -> controller / controlled / controlType

    Physical entities are stored once in an entity registry. Event roles and
    controls reference entities by BioPAX UID, preserving identity even when
    the same Protein/Complex appears in multiple biological roles.
    """

    def __init__(self, reactome_db) -> None:
        self.db = reactome_db

    @staticmethod
    def _as_list(value: Any) -> List[Any]:
        if value is None:
            return []
        if isinstance(value, (list, tuple, set)):
            return list(value)
        return [value]

    @staticmethod
    def _class_names(obj: Any) -> Set[str]:
        return {cls.__name__ for cls in type(obj).mro()}

    @staticmethod
    def _uid(obj: Any) -> Optional[str]:
        return getattr(obj, "uid", None)

    @staticmethod
    def _display_name(obj: Any) -> Optional[str]:
        return getattr(obj, "display_name", None)

    @staticmethod
    def _terms(obj: Any) -> List[str]:
        if obj is None:
            return []
        return [str(v) for v in (getattr(obj, "term", []) or [])]

    def _extract_xrefs(self, obj: Any) -> List[Dict[str, Any]]:
        result: List[Dict[str, Any]] = []

        for xref in self._as_list(getattr(obj, "xref", [])):
            result.append(
                {
                    "uid": self._uid(xref),
                    "type": type(xref).__name__,
                    "db": getattr(xref, "db", None),
                    "id": getattr(xref, "id", None),
                    "title": getattr(xref, "title", None),
                    "year": getattr(xref, "year", None),
                    "url": self._as_list(getattr(xref, "url", [])),
                }
            )

        return result

    def _extract_ec_numbers(self, obj: Any) -> List[str]:
        """
        Extract EC numbers explicitly attached to a BioPAX object.

        This intentionally does not infer EC numbers from GO terms or from the
        controlled BiochemicalReaction. Only direct EC annotations on the object
        itself are returned.
        """
        result: List[str] = []

        # Future-proofing: use a direct BioPAX/pybiopax attribute if present.
        for value in self._as_list(
            getattr(obj, "e_c_number", [])
        ):
            value = str(value).strip()
            if value and value not in result:
                result.append(value)

        # Also inspect explicit xrefs whose database is EC.
        for xref in self._as_list(
            getattr(obj, "xref", [])
        ):
            db = getattr(xref, "db", None)
            identifier = getattr(xref, "id", None)

            if db is None or identifier is None:
                continue

            db = str(db).strip().lower()
            identifier = str(identifier).strip()

            if db in {
                "ec",
                "enzyme commission",
                "enzyme nomenclature",
            }:
                identifier = re.sub(
                    r"^EC[:\s]*",
                    "",
                    identifier,
                    flags=re.IGNORECASE,
                )

                if identifier and identifier not in result:
                    result.append(identifier)

        return result

    def _extract_data_sources(self, obj: Any) -> List[Dict[str, Any]]:
        result = []

        for source in self._as_list(getattr(obj, "data_source", [])):
            result.append(
                {
                    "uid": self._uid(source),
                    "type": type(source).__name__,
                    "display_name": self._display_name(source),
                    "name": self._as_list(getattr(source, "name", [])),
                    "xref": self._extract_xrefs(source),
                }
            )

        return result

    def _extract_evidence(self, obj: Any) -> List[Dict[str, Any]]:
        result = []

        for evidence in self._as_list(getattr(obj, "evidence", [])):
            result.append(
                {
                    "uid": self._uid(evidence),
                    "type": type(evidence).__name__,
                    "xref": self._extract_xrefs(evidence),
                    "comments": self._as_list(
                        getattr(evidence, "comment", [])
                    ),
                }
            )

        return result

    def _extract_annotations(self, obj: Any) -> Dict[str, Any]:
        return {
            "xref": self._extract_xrefs(obj),
            "comments": [
                str(v)
                for v in self._as_list(getattr(obj, "comment", []))
            ],
            "data_source": self._extract_data_sources(obj),
            "evidence": self._extract_evidence(obj),
            "interaction_type": [
                self._terms(v)
                for v in self._as_list(
                    getattr(obj, "interaction_type", [])
                )
            ],
        }

    def _extract_sequence_location(
        self,
        location: Any,
    ) -> Optional[Dict[str, Any]]:
        if location is None:
            return None

        result: Dict[str, Any] = {
            "uid": self._uid(location),
            "type": type(location).__name__,
        }

        if hasattr(location, "sequence_position"):
            result["position"] = getattr(
                location,
                "sequence_position",
                None,
            )
            result["position_status"] = getattr(
                location,
                "position_status",
                None,
            )

        begin = getattr(
            location,
            "sequence_interval_begin",
            None,
        )
        end = getattr(
            location,
            "sequence_interval_end",
            None,
        )

        if begin is not None or end is not None:
            result["begin"] = (
                self._extract_sequence_location(begin)
                if begin is not None
                else None
            )
            result["end"] = (
                self._extract_sequence_location(end)
                if end is not None
                else None
            )

        return result

    def _extract_features(self, entity: Any) -> List[Dict[str, Any]]:
        features: List[Dict[str, Any]] = []
        all_features = []

        for feature in self._as_list(getattr(entity, "feature", [])):
            all_features.append(("feature", feature))

        for feature in self._as_list(getattr(entity, "not_feature", [])):
            all_features.append(("not_feature", feature))

        for relation, feature in all_features:
            modification_type = getattr(
                feature,
                "modification_type",
                None,
            )

            features.append(
                {
                    "uid": self._uid(feature),
                    "type": type(feature).__name__,
                    "relation": relation,
                    "modification_type": self._terms(
                        modification_type
                    ),
                    "location": self._extract_sequence_location(
                        getattr(
                            feature,
                            "feature_location",
                            None,
                        )
                    ),
                    "xref": self._extract_xrefs(feature),
                    "comments": self._as_list(
                        getattr(feature, "comment", [])
                    ),
                }
            )

        return features

    def _extract_entity_reference(
        self,
        entity: Any,
    ) -> Optional[Dict[str, Any]]:
        ref = getattr(entity, "entity_reference", None)

        if ref is None:
            return None

        return {
            "uid": self._uid(ref),
            "type": type(ref).__name__,
            "display_name": self._display_name(ref),
            "xref": self._extract_xrefs(ref),
            "sequence": getattr(ref, "sequence", None),
            "organism_uid": self._uid(
                getattr(ref, "organism", None)
            ),
            "comments": self._as_list(
                getattr(ref, "comment", [])
            ),
        }

    @staticmethod
    def _uniprot_accessions(
        reference: Optional[Dict[str, Any]],
    ) -> List[str]:
        if not reference:
            return []

        values = []

        for xref in reference.get("xref", []):
            if str(xref.get("db", "")).lower() == "uniprot":
                value = xref.get("id")
                if value:
                    values.append(str(value))

        return list(dict.fromkeys(values))

    @staticmethod
    def _component_stoichiometry(parent: Any, child: Any) -> Any:
        for stoich in (
            getattr(parent, "component_stoichiometry", [])
            or []
        ):
            if getattr(stoich, "physical_entity", None) is child:
                return getattr(
                    stoich,
                    "stoichiometric_coefficient",
                    None,
                )

        return None

    def _register_entity(
        self,
        entity: Any,
        registry: Dict[str, Dict[str, Any]],
        active: Optional[Set[str]] = None,
    ) -> Optional[str]:
        if entity is None:
            return None

        uid = self._uid(entity)
        if uid is None:
            uid = f"anonymous:{id(entity)}"

        if uid in registry:
            return uid

        if active is None:
            active = set()

        if uid in active:
            return uid

        next_active = set(active)
        next_active.add(uid)

        reference = self._extract_entity_reference(entity)
        location = getattr(entity, "cellular_location", None)

        record: Dict[str, Any] = {
            "uid": uid,
            "entity_type": type(entity).__name__,
            "class_hierarchy": sorted(self._class_names(entity)),
            "display_name": self._display_name(entity),
            "cellular_location": {
                "uid": self._uid(location),
                "terms": self._terms(location),
            } if location is not None else None,
            "entity_reference": reference,
            "uniprot_accessions": self._uniprot_accessions(reference),
            "features": self._extract_features(entity),
            "xref": self._extract_xrefs(entity),
            "comments": self._as_list(
                getattr(entity, "comment", [])
            ),
            "children": [],
        }

        registry[uid] = record

        for component in self._as_list(
            getattr(entity, "component", [])
        ):
            child_uid = self._register_entity(
                component,
                registry,
                next_active,
            )

            record["children"].append(
                {
                    "relation": "component",
                    "entity_uid": child_uid,
                    "stoichiometry": self._component_stoichiometry(
                        entity,
                        component,
                    ),
                }
            )

        for member in self._as_list(
            getattr(entity, "member_physical_entity", [])
        ):
            child_uid = self._register_entity(
                member,
                registry,
                next_active,
            )

            record["children"].append(
                {
                    "relation": "member_physical_entity",
                    "entity_uid": child_uid,
                    "stoichiometry": None,
                }
            )

        return uid

    @staticmethod
    def _participant_stoichiometry(event: Any) -> Dict[Any, Any]:
        result = {}

        for stoich in (
            getattr(event, "participant_stoichiometry", [])
            or []
        ):
            entity = getattr(stoich, "physical_entity", None)

            if entity is not None:
                result[entity] = getattr(
                    stoich,
                    "stoichiometric_coefficient",
                    None,
                )

        return result

    def _extract_roles(
        self,
        event: Any,
        registry: Dict[str, Dict[str, Any]],
    ) -> Dict[str, List[Dict[str, Any]]]:
        class_names = self._class_names(event)
        roles: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
        stoich_map = self._participant_stoichiometry(event)

        if "Conversion" in class_names:
            for role, attr in (
                ("LEFT", "left"),
                ("RIGHT", "right"),
            ):
                for entity in self._as_list(
                    getattr(event, attr, [])
                ):
                    entity_uid = self._register_entity(
                        entity,
                        registry,
                    )

                    roles[role].append(
                        {
                            "entity_uid": entity_uid,
                            "stoichiometry": stoich_map.get(entity),
                        }
                    )

        elif "TemplateReaction" in class_names:
            for role, attr in (
                ("TEMPLATE", "template"),
                ("PRODUCT", "product"),
            ):
                for entity in self._as_list(
                    getattr(event, attr, [])
                ):
                    entity_uid = self._register_entity(
                        entity,
                        registry,
                    )

                    roles[role].append(
                        {
                            "entity_uid": entity_uid,
                            "stoichiometry": None,
                        }
                    )

        else:
            for entity in self._as_list(
                getattr(event, "participant", [])
            ):
                entity_uid = self._register_entity(
                    entity,
                    registry,
                )

                roles["PARTICIPANT"].append(
                    {
                        "entity_uid": entity_uid,
                        "stoichiometry": None,
                    }
                )

        return dict(roles)

    def _extract_control(
        self,
        control: Any,
        registry: Dict[str, Dict[str, Any]],
        visited: Optional[Set[str]] = None,
    ) -> Dict[str, Any]:
        if visited is None:
            visited = set()

        uid = self._uid(control)
        key = uid or f"anonymous:{id(control)}"

        if key in visited:
            return {
                "uid": uid,
                "control_class": type(control).__name__,
                "cycle_reference": True,
            }

        next_visited = set(visited)
        next_visited.add(key)

        controllers = []

        for controller in self._as_list(
            getattr(control, "controller", [])
        ):
            entity_uid = self._register_entity(
                controller,
                registry,
            )
            controllers.append(entity_uid)

        controlled = getattr(control, "controlled", None)

        nested_controls = [
            self._extract_control(
                nested,
                registry,
                next_visited,
            )
            for nested in self._as_list(
                getattr(control, "_controlled_of", [])
            )
        ]

        return {
            "uid": uid,
            "control_class": type(control).__name__,
            "class_hierarchy": sorted(
                self._class_names(control)
            ),
            "display_name": self._display_name(control),
            "control_type": getattr(
                control,
                "control_type",
                None,
            ),
            "controller_uids": controllers,
            "controlled_uid": self._uid(controlled),
            "controlled_type": (
                type(controlled).__name__
                if controlled is not None
                else None
            ),
            "xref": self._extract_xrefs(control),
            "comments": self._as_list(
                getattr(control, "comment", [])
            ),
            "controlled_by": nested_controls,
        }

    def _extract_controls(
        self,
        event: Any,
        registry: Dict[str, Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        return [
            self._extract_control(
                control,
                registry,
            )
            for control in self._as_list(
                getattr(event, "_controlled_of", [])
            )
        ]

    def _event_specific_properties(
        self,
        event: Any,
    ) -> Dict[str, Any]:
        class_names = self._class_names(event)
        result: Dict[str, Any] = {}

        if "Conversion" in class_names:
            result.update(
                {
                    "conversion_direction": getattr(
                        event,
                        "conversion_direction",
                        None,
                    ),
                    "spontaneous": getattr(
                        event,
                        "spontaneous",
                        None,
                    ),
                }
            )

        if "BiochemicalReaction" in class_names:
            for output_key, attr in (
                ("delta_g", "delta_g"),
                ("delta_h", "delta_h"),
                ("delta_s", "delta_s"),
                ("k_eq", "k_e_q"),
                ("e_c_number", "e_c_number"),
            ):
                result[output_key] = self._as_list(
                    getattr(event, attr, [])
                )

        if "TemplateReaction" in class_names:
            result["template_direction"] = getattr(
                event,
                "template_direction",
                None,
            )

        if hasattr(event, "interaction_score"):
            result["interaction_score"] = self._as_list(
                getattr(event, "interaction_score", [])
            )

        if hasattr(event, "phenotype"):
            phenotype = getattr(event, "phenotype", None)
            result["phenotype"] = (
                {
                    "uid": self._uid(phenotype),
                    "type": type(phenotype).__name__,
                    "terms": self._terms(phenotype),
                }
                if phenotype is not None
                else None
            )

        return result

    def _pathway_context(
        self,
        event_stid: str,
    ) -> List[List[Dict[str, Any]]]:
        self.db.build_pathway_index()

        result = []

        for chain in self.db.get_event_pathway_chains(
            event_stid
        ):
            result.append(
                [
                    {
                        "stable_id": pathway_stid,
                        "uid": self._uid(
                            self.db.pathway_index[pathway_stid]
                        ),
                        "display_name": self._display_name(
                            self.db.pathway_index[pathway_stid]
                        ),
                    }
                    for pathway_stid in chain
                ]
            )

        return result

    def build(self, event_stid: str) -> Dict[str, Any]:
        event = self.db.get_event(event_stid)

        if event is None:
            raise KeyError(
                f"Reactome event not found in local BioPAX: "
                f"{event_stid}"
            )

        registry: Dict[str, Dict[str, Any]] = {}

        return {
            "stable_id": event_stid,
            "uid": self._uid(event),
            "event_type": type(event).__name__,
            "class_hierarchy": sorted(self._class_names(event)),
            "display_name": self._display_name(event),
            "pathway_chains": self._pathway_context(event_stid),
            "roles": self._extract_roles(event, registry),
            "entities": registry,
            "controls": self._extract_controls(event, registry),
            "annotations": self._extract_annotations(event),
            "event_specific": self._event_specific_properties(event),
            "disease": self.db.get_event_disease_context(event_stid),
        }

