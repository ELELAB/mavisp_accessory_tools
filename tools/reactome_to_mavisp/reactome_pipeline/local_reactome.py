from __future__ import annotations

from collections import defaultdict
from typing import Any, Dict, List, Set

import pandas as pd


class LocalReactomeDatabase:
    """
    Local Reactome/BioPAX data-access layer.

    Responsibilities:
    - load UniProt2Reactome_PE_Reactions.txt
    - load Homo_sapiens.owl once
    - dynamically index Reactome BioPAX events by R-HSA stable ID
    - index pathways and reconstruct pathway ancestry
    - expose BioPAX objects to higher-level EventContext builders

    Biological interpretation of an event is intentionally kept out of this
    class and belongs in event_context.py.
    """

    def __init__(
        self,
        reaction_map_file: str,
        biopax_file: str,
    ) -> None:
        self.reaction_map_file = reaction_map_file
        self.biopax_file = biopax_file

        self.uniprot_to_reactions: Dict[str, Set[str]] = defaultdict(set)

        self.model = None
        self.event_index: Dict[str, Any] = {}
        self.pathway_index: Dict[str, Any] = {}
        self.pathway_parent_index: Dict[str, Set[str]] = {}

        self._load_reaction_mapping()

    def _load_reaction_mapping(self) -> None:
        print(
            f"[INFO] Loading Reactome reaction mapping: "
            f"{self.reaction_map_file}"
        )

        df = pd.read_csv(
            self.reaction_map_file,
            sep="\t",
            header=None,
            dtype=str,
        )

        if df.shape[1] < 8:
            raise ValueError(
                "Unexpected UniProt2Reactome_PE_Reactions.txt format: "
                f"expected at least 8 columns, found {df.shape[1]}"
            )

        df = df[df[7] == "Homo sapiens"].copy()
        df = df[df[3].str.startswith("R-HSA-", na=False)]

        for uniprot_ac, reaction_id in zip(df[0], df[3]):
            self.uniprot_to_reactions[uniprot_ac].add(reaction_id)

        print(
            f"[INFO] Loaded reaction mappings for "
            f"{len(self.uniprot_to_reactions)} human UniProt accessions"
        )

    def get_reactions_for_uniprot(self, uniprot_ac: str) -> List[str]:
        return sorted(self.uniprot_to_reactions.get(uniprot_ac, set()))

    def has_uniprot(self, uniprot_ac: str) -> bool:
        return bool(self.uniprot_to_reactions.get(uniprot_ac))

    def load_biopax_model(self) -> None:
        if self.model is not None:
            return

        print(f"[INFO] Loading BioPAX model: {self.biopax_file}")

        import pybiopax

        self.model = pybiopax.api.model_from_owl_file(self.biopax_file)

        print(
            f"[INFO] BioPAX model loaded: "
            f"{len(self.model.objects)} objects"
        )

    @staticmethod
    def _class_names(obj: Any) -> Set[str]:
        return {cls.__name__ for cls in type(obj).mro()}

    @staticmethod
    def _reactome_stable_ids(obj: Any) -> List[str]:
        stable_ids: List[str] = []

        for xref in getattr(obj, "xref", []) or []:
            db = getattr(xref, "db", "")
            xref_id = getattr(xref, "id", "")

            if (
                isinstance(xref_id, str)
                and xref_id.startswith("R-HSA-")
                and db
                and "reactome" in str(db).lower()
            ):
                stable_ids.append(xref_id)

        return list(dict.fromkeys(stable_ids))

    def build_event_index(self) -> None:
        """
        Dynamically index BioPAX events.

        Any R-HSA object that is an Interaction (or subclass) but not a
        Control (or subclass) is considered a primary event.
        """

        self.load_biopax_model()

        if self.event_index:
            return

        print("[INFO] Building generic BioPAX event index...")

        counts: Dict[str, int] = defaultdict(int)

        for obj in self.model.objects.values():
            class_names = self._class_names(obj)

            if "Interaction" not in class_names:
                continue

            if "Control" in class_names:
                continue

            for stable_id in self._reactome_stable_ids(obj):
                self.event_index[stable_id] = obj
                counts[type(obj).__name__] += 1

        print(f"[INFO] Indexed {len(self.event_index)} Reactome events")

        for event_type, count in sorted(counts.items()):
            print(f"[INFO]   {event_type}: {count}")

    def get_event(self, event_stid: str):
        if not self.event_index:
            self.build_event_index()
        return self.event_index.get(event_stid)

    def get_detected_event_types(self) -> Dict[str, int]:
        if not self.event_index:
            self.build_event_index()

        counts: Dict[str, int] = defaultdict(int)
        for event in self.event_index.values():
            counts[type(event).__name__] += 1

        return dict(sorted(counts.items()))

    def build_pathway_index(self) -> None:
        self.load_biopax_model()

        if self.pathway_index:
            return

        import pybiopax

        print("[INFO] Building BioPAX pathway index...")

        pathways = self.model.get_objects_by_type(
            pybiopax.biopax.Pathway
        )

        for pathway in pathways:
            for stable_id in self._reactome_stable_ids(pathway):
                self.pathway_index[stable_id] = pathway

        print(
            f"[INFO] Indexed "
            f"{len(self.pathway_index)} Reactome pathways"
        )

    def get_pathway(self, pathway_stid: str):
        if not self.pathway_index:
            self.build_pathway_index()
        return self.pathway_index.get(pathway_stid)

    def get_pathway_reaction_order(
        self,
        pathway_stid: str,
    ) -> Dict[str, List[Dict[str, List[str]]]]:
        """
        Return local reaction-order relationships for one Reactome pathway.

        The ordering is reconstructed from BioPAX PathwayStep objects using:

            reaction.step_process_of
            PathwayStep.next_step
            PathwayStep.next_step_of
            PathwayStep.step_process

        Only PathwaySteps belonging to the requested pathway are considered.

        Reaction stable IDs are used as graph nodes.

        Returns
        -------
        dict
            Example:

            {
                "R-HSA-AAAA": [
                    {
                        "next_step": ["R-HSA-BBBB"],
                        "previous": ["R-HSA-CCCC"],
                    }
                ]
            }
        """
        pathway = self.get_pathway(pathway_stid)

        if pathway is None:
            return {}

        # Make sure we know which stable IDs correspond to primary events.
        self.build_event_index()

        # PathwayStep objects explicitly belonging to this pathway.
        pathway_steps = set(
            getattr(
                pathway,
                "pathway_order",
                [],
            )
            or []
        )

        result: Dict[
            str,
            List[Dict[str, List[str]]]
        ] = {}

        # Direct components of the pathway can include reactions and,
        # in some cases, nested pathways. We only want primary events.
        pathway_components = (
            getattr(
                pathway,
                "pathway_component",
                [],
            )
            or []
        )

        for event in pathway_components:

            event_stids = [
                stable_id
                for stable_id in self._reactome_stable_ids(
                    event
                )
                if stable_id in self.event_index
            ]

            if not event_stids:
                continue

            # Normally there is one Reactome stable ID per event.
            event_stid = event_stids[0]

            step_records: List[
                Dict[str, List[str]]
            ] = []

            for step in (
                getattr(
                    event,
                    "step_process_of",
                    [],
                )
                or []
            ):

                # A reaction can occur in more than one pathway.
                # Restrict the relationship to PathwaySteps of the
                # requested pathway.
                if (
                    pathway_steps
                    and step not in pathway_steps
                ):
                    continue

                record: Dict[
                    str,
                    List[str]
                ] = {}

                # -------------------------------------------------
                # Next reactions
                # -------------------------------------------------

                next_reactions: List[str] = []

                for next_step in (
                    getattr(
                        step,
                        "next_step",
                        [],
                    )
                    or []
                ):

                    if (
                        pathway_steps
                        and next_step not in pathway_steps
                    ):
                        continue

                    for process in (
                        getattr(
                            next_step,
                            "step_process",
                            [],
                        )
                        or []
                    ):

                        for process_stid in (
                            self._reactome_stable_ids(
                                process
                            )
                        ):

                            if (
                                process_stid
                                not in self.event_index
                            ):
                                continue

                            if (
                                process_stid
                                not in next_reactions
                            ):
                                next_reactions.append(
                                    process_stid
                                )

                if next_reactions:
                    record[
                        "next_step"
                    ] = next_reactions

                # -------------------------------------------------
                # Previous reactions
                # -------------------------------------------------

                previous_reactions: List[str] = []

                for previous_step in (
                    getattr(
                        step,
                        "next_step_of",
                        [],
                    )
                    or []
                ):

                    if (
                        pathway_steps
                        and previous_step
                        not in pathway_steps
                    ):
                        continue

                    for process in (
                        getattr(
                            previous_step,
                            "step_process",
                            [],
                        )
                        or []
                    ):

                        for process_stid in (
                            self._reactome_stable_ids(
                                process
                            )
                        ):

                            if (
                                process_stid
                                not in self.event_index
                            ):
                                continue

                            if (
                                process_stid
                                not in previous_reactions
                            ):
                                previous_reactions.append(
                                    process_stid
                                )

                if previous_reactions:
                    record[
                        "previous"
                    ] = previous_reactions

                step_records.append(
                    record
                )

            result[event_stid] = (
                step_records
            )

        return result

    def build_pathway_parent_index(self) -> None:
        self.build_pathway_index()

        if self.pathway_parent_index:
            return

        print("[INFO] Building pathway parent index...")

        object_to_stid = {
            pathway: stid
            for stid, pathway in self.pathway_index.items()
        }

        for parent_stid, parent_pathway in self.pathway_index.items():
            for component in getattr(
                parent_pathway,
                "pathway_component",
                [],
            ) or []:
                if component in object_to_stid:
                    child_stid = object_to_stid[component]
                    self.pathway_parent_index.setdefault(
                        child_stid,
                        set(),
                    ).add(parent_stid)

        print(
            f"[INFO] Indexed parent relationships for "
            f"{len(self.pathway_parent_index)} pathways"
        )

    def get_parent_pathways(self, pathway_stid: str) -> List[str]:
        if not self.pathway_parent_index:
            self.build_pathway_parent_index()

        return sorted(
            self.pathway_parent_index.get(pathway_stid, set())
        )

    def get_pathway_ancestry(
        self,
        pathway_stid: str,
    ) -> List[List[str]]:
        if not self.pathway_parent_index:
            self.build_pathway_parent_index()

        chains: List[List[str]] = []

        def _walk(
            current_stid: str,
            current_chain: List[str],
            visited: Set[str],
        ) -> None:
            if current_stid in visited:
                return

            next_visited = set(visited)
            next_visited.add(current_stid)

            parents = self.pathway_parent_index.get(
                current_stid,
                set(),
            )

            if not parents:
                chains.append(list(reversed(current_chain)))
                return

            for parent_stid in sorted(parents):
                _walk(
                    parent_stid,
                    current_chain + [parent_stid],
                    next_visited,
                )

        _walk(pathway_stid, [pathway_stid], set())
        return chains

    def get_direct_pathways_for_event(
        self,
        event_stid: str,
    ) -> List[str]:
        self.build_event_index()
        self.build_pathway_index()

        event = self.event_index.get(event_stid)

        if event is None:
            return []

        object_to_stid = {
            pathway: stid
            for stid, pathway in self.pathway_index.items()
        }

        direct_objects = getattr(
            event,
            "_pathway_component_of",
            None,
        )

        if direct_objects:
            return sorted(
                {
                    object_to_stid[pathway]
                    for pathway in direct_objects
                    if pathway in object_to_stid
                }
            )

        direct_pathways: List[str] = []

        for pathway_stid, pathway in self.pathway_index.items():
            if event in (
                getattr(pathway, "pathway_component", []) or []
            ):
                direct_pathways.append(pathway_stid)

        return sorted(set(direct_pathways))

    def get_event_pathway_chains(
        self,
        event_stid: str,
    ) -> List[List[str]]:
        chains: List[List[str]] = []

        for pathway_stid in self.get_direct_pathways_for_event(
            event_stid
        ):
            chains.extend(
                self.get_pathway_ancestry(pathway_stid)
            )

        seen = set()
        unique: List[List[str]] = []

        for chain in chains:
            key = tuple(chain)
            if key not in seen:
                seen.add(key)
                unique.append(chain)

        return unique

    def build_target_pathway_structure(
        self,
        uniprot_ac: str,
    ) -> List[Dict[Any, Any]]:
        reaction_ids = self.get_reactions_for_uniprot(uniprot_ac)
        chain_to_events: Dict[tuple, Set[str]] = {}

        for event_stid in reaction_ids:
            for chain in self.get_event_pathway_chains(event_stid):
                chain_to_events.setdefault(
                    tuple(chain),
                    set(),
                ).add(event_stid)

        result: List[Dict[Any, Any]] = []

        for chain, event_ids in chain_to_events.items():
            entry: Dict[Any, Any] = {}

            for level, pathway_stid in enumerate(chain):
                pathway = self.pathway_index[pathway_stid]
                entry[level] = {
                    "id": pathway_stid,
                    "uid": getattr(pathway, "uid", None),
                    "name": getattr(
                        pathway,
                        "display_name",
                        None,
                    ),
                }

            entry["reactions"] = sorted(event_ids)
            result.append(entry)

        return result
    def get_event_disease_context(
        self,
        event_stid: str,
        disease_root_stid: str = "R-HSA-1643685",
    ) -> dict:
        """
        Return disease-pathway context for a Reactome event.

        Disease membership is inferred from pathway ancestry:
        an event is considered disease-associated if at least one
        pathway chain contains the Reactome top-level Disease pathway
        R-HSA-1643685.

        Parameters
        ----------
        event_stid : str
            Reactome stable ID of the event.
        disease_root_stid : str
            Reactome stable ID of the top-level Disease pathway.

        Returns
        -------
        dict
            Dictionary with:
            - is_in_disease
            - disease_root
            - disease_chains
            - disease_pathways
            - lowest_disease_pathways
            - mondo
        """

        chains = self.get_event_pathway_chains(event_stid)

        disease_chains = []

        for chain in chains:
            if disease_root_stid not in chain:
                continue

            root_index = chain.index(disease_root_stid)

            # Keep only the Disease branch,
            # starting from top-level Disease.
            disease_branch = chain[root_index:]

            disease_chains.append(disease_branch)

        if not disease_chains:
            return {
                "is_in_disease": False,
                "disease_root": None,
                "disease_chains": [],
                "disease_pathways": [],
                "lowest_disease_pathways": [],
                "mondo": [],
            }

        disease_root = self.get_pathway(disease_root_stid)

        root_info = {
            "stable_id": disease_root_stid,
            "uid": getattr(disease_root, "uid", None),
            "display_name": getattr(disease_root, "display_name", None),
        }

        # Collect unique disease pathways while preserving order.
        seen_pathways = set()
        disease_pathways = []

        for chain in disease_chains:
            for pathway_stid in chain:

                if pathway_stid in seen_pathways:
                    continue

                pathway = self.get_pathway(pathway_stid)

                disease_pathways.append(
                    {
                        "stable_id": pathway_stid,
                        "uid": getattr(pathway, "uid", None),
                        "display_name": getattr(
                            pathway,
                            "display_name",
                            None,
                        ),
                    }
                )

                seen_pathways.add(pathway_stid)

        # One event may belong to multiple disease branches.
        # Therefore preserve ALL lowest disease pathways.
        lowest_seen = set()
        lowest_disease_pathways = []

        for chain in disease_chains:

            if not chain:
                continue

            lowest_stid = chain[-1]

            if lowest_stid in lowest_seen:
                continue

            pathway = self.get_pathway(lowest_stid)

            lowest_disease_pathways.append(
                {
                    "stable_id": lowest_stid,
                    "uid": getattr(pathway, "uid", None),
                    "display_name": getattr(
                        pathway,
                        "display_name",
                        None,
                    ),
                }
            )

            lowest_seen.add(lowest_stid)

        # MONDO is intentionally empty here.
        # In the tested Reactome BioPAX release no MONDO xrefs
        # were serialized on pathways or events.
        mondo = []

        return {
            "is_in_disease": True,
            "disease_root": root_info,
            "disease_chains": disease_chains,
            "disease_pathways": disease_pathways,
            "lowest_disease_pathways": lowest_disease_pathways,
            "mondo": mondo,
        }