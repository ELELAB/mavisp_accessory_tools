from __future__ import annotations

from collections import defaultdict
from typing import Any, Dict, Iterable, List, Optional, Set

import pandas as pd


class DiseaseVariantIndex:
    """
    Local index over Reactome disease_variant_ewas_mapping.tsv.

    The TSV is used to discover disease-specific Reactome reactions for a
    UniProt accession and to attach variant/disease metadata that are not
    represented as generic pathway metadata in the BioPAX EventContext.

    One disease reaction may contain multiple variants of the same protein.
    Therefore the index stores a list of variant records per reaction.
    """

    REACTION_COL = (
        "entityWithAccessionedSequence_reactionLikeEvent_stable_id"
    )
    REACTION_NAME_COL = (
        "entityWithAccessionedSequence_reactionLikeEvent_displayName"
    )
    PATHWAY_COL = (
        "entityWithAccessionedSequence_pathway_stable_id"
    )
    PATHWAY_NAME_COL = (
        "entityWithAccessionedSequence_pathway_displayName"
    )
    NORMAL_REACTION_COL = (
        "entityWithAccessionedSequence_reactionLikeEvent_normalReaction_stable_id"
    )
    NORMAL_REACTION_NAME_COL = (
        "entityWithAccessionedSequence_reactionLikeEvent_normalReaction_displayName"
    )
    NORMAL_PATHWAY_COL = (
        "entityWithAccessionedSequence_pathway_normalPathway_stable_id"
    )
    NORMAL_PATHWAY_NAME_COL = (
        "entityWithAccessionedSequence_pathway_normalPathway_displayName"
    )
    FUNCTIONAL_STATUS_COL = (
        "entityWithAccessionedSequence_reactionLikeEvent_"
        "entityFunctionalStatus_functionalStatus_functionalStatusType_displayName"
    )
    PUBMED_COL = (
        "reactionLikeEvent_literatureReference_pubMedIdentifier"
    )

    def __init__(self, mapping_file: str) -> None:
        self.mapping_file = mapping_file

        # UniProt -> disease reaction IDs
        self.uniprot_to_reactions: Dict[str, Set[str]] = defaultdict(set)

        # disease reaction ID -> variant metadata records
        self.reaction_to_variants: Dict[
            str, List[Dict[str, Any]]
        ] = defaultdict(list)

        # UniProt -> all variant records
        self.uniprot_to_variants: Dict[
            str, List[Dict[str, Any]]
        ] = defaultdict(list)

        self._load()

    @staticmethod
    def _clean(value: Any) -> Optional[str]:
        if value is None or pd.isna(value):
            return None

        value = str(value).strip()
        return value if value else None

    @classmethod
    def _split_pipe(cls, value: Any) -> List[str]:
        value = cls._clean(value)

        if value is None:
            return []

        return [
            item.strip()
            for item in value.split("|")
            if item.strip()
        ]

    @staticmethod
    def _uniprot_from_reference_entity_id(
        value: Optional[str],
    ) -> Optional[str]:
        if not value:
            return None

        if value.lower().startswith("uniprot:"):
            return value.split(":", 1)[1]

        return None

    def _load(self) -> None:
        print(
            f"[INFO] Loading Reactome disease variant mapping: "
            f"{self.mapping_file}"
        )

        df = pd.read_csv(
            self.mapping_file,
            sep="\t",
            dtype=str,
            keep_default_na=False,
        )

        required = {
            "Genename",
            "displayName",
            "stable_id",
            "referenceEntity_name",
            "referenceEntity_id",
            "hasModifiedResidue_displayName",
            "modifiedResidue_class",
            "cross_reference",
            "disease",
            "disease_identifier",
            self.REACTION_COL,
            self.REACTION_NAME_COL,
            self.FUNCTIONAL_STATUS_COL,
            self.PATHWAY_COL,
            self.PATHWAY_NAME_COL,
            self.NORMAL_REACTION_COL,
            self.NORMAL_REACTION_NAME_COL,
            self.NORMAL_PATHWAY_COL,
            self.NORMAL_PATHWAY_NAME_COL,
        }

        missing = sorted(required.difference(df.columns))

        if missing:
            raise ValueError(
                "Unexpected disease_variant_ewas_mapping.tsv format. "
                f"Missing columns: {missing}"
            )

        for _, row in df.iterrows():
            reference_entity_id = self._clean(
                row.get("referenceEntity_id")
            )
            uniprot_ac = self._uniprot_from_reference_entity_id(
                reference_entity_id
            )

            if not uniprot_ac:
                continue

            reaction_ids = self._split_pipe(
                row.get(self.REACTION_COL)
            )

            if not reaction_ids:
                continue

            record: Dict[str, Any] = {
                "gene": self._clean(row.get("Genename")),
                "variant_display_name": self._clean(
                    row.get("displayName")
                ),
                "variant_entity_stable_id": self._clean(
                    row.get("stable_id")
                ),
                "reference_entity_name": self._clean(
                    row.get("referenceEntity_name")
                ),
                "reference_entity_id": reference_entity_id,
                "uniprot_ac": uniprot_ac,
                "modification_description": self._clean(
                    row.get("hasModifiedResidue_displayName")
                ),
                "modification_class": self._clean(
                    row.get("modifiedResidue_class")
                ),
                "cross_references": self._split_pipe(
                    row.get("cross_reference")
                ),
                "diseases": self._split_pipe(
                    row.get("disease")
                ),
                "disease_identifiers": self._split_pipe(
                    row.get("disease_identifier")
                ),
                "functional_status": self._clean(
                    row.get(self.FUNCTIONAL_STATUS_COL)
                ),
                "literature_pubmed": self._split_pipe(
                    row.get(self.PUBMED_COL)
                ),
                "disease_reaction_ids": reaction_ids,
                "disease_reaction_names": self._split_pipe(
                    row.get(self.REACTION_NAME_COL)
                ),
                "disease_pathway_ids": self._split_pipe(
                    row.get(self.PATHWAY_COL)
                ),
                "disease_pathway_names": self._split_pipe(
                    row.get(self.PATHWAY_NAME_COL)
                ),
                "normal_reaction_ids": self._split_pipe(
                    row.get(self.NORMAL_REACTION_COL)
                ),
                "normal_reaction_names": self._split_pipe(
                    row.get(self.NORMAL_REACTION_NAME_COL)
                ),
                "normal_pathway_ids": self._split_pipe(
                    row.get(self.NORMAL_PATHWAY_COL)
                ),
                "normal_pathway_names": self._split_pipe(
                    row.get(self.NORMAL_PATHWAY_NAME_COL)
                ),
            }

            self.uniprot_to_variants[uniprot_ac].append(record)

            for reaction_id in reaction_ids:
                if not reaction_id.startswith("R-HSA-"):
                    continue

                self.uniprot_to_reactions[uniprot_ac].add(
                    reaction_id
                )
                self.reaction_to_variants[reaction_id].append(
                    record
                )

        print(
            f"[INFO] Loaded disease mappings for "
            f"{len(self.uniprot_to_reactions)} UniProt accessions"
        )
        print(
            f"[INFO] Indexed "
            f"{len(self.reaction_to_variants)} disease reactions"
        )

    def has_uniprot(self, uniprot_ac: str) -> bool:
        return bool(
            self.uniprot_to_reactions.get(uniprot_ac)
        )

    def get_reactions_for_uniprot(
        self,
        uniprot_ac: str,
    ) -> List[str]:
        return sorted(
            self.uniprot_to_reactions.get(
                uniprot_ac,
                set(),
            )
        )

    def get_variants_for_uniprot(
        self,
        uniprot_ac: str,
    ) -> List[Dict[str, Any]]:
        return list(
            self.uniprot_to_variants.get(
                uniprot_ac,
                [],
            )
        )

    def get_variants_for_reaction(
        self,
        reaction_stid: str,
        uniprot_ac: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        records = self.reaction_to_variants.get(
            reaction_stid,
            [],
        )

        if uniprot_ac is None:
            return list(records)

        return [
            record
            for record in records
            if record.get("uniprot_ac") == uniprot_ac
        ]

