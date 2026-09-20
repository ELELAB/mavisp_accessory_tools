from __future__ import annotations

from typing import Any, Dict, List, Optional

from reactome_pipeline.event_context import EventContextBuilder
from reactome_pipeline.disease_variants import DiseaseVariantIndex


class ReactomeEventContextFactory:
    """
    Build independent NORMAL and DISEASE EventContexts for one UniProt target.

    Both branches use the same EventContextBuilder and the same local BioPAX
    model. The difference is event discovery:

    NORMAL:
        UniProt2Reactome_PE_Reactions.txt

    DISEASE:
        disease_variant_ewas_mapping.tsv

    Disease EventContexts receive additional disease_variants metadata from
    the TSV.
    """

    def __init__(
        self,
        reactome_db,
        disease_variant_index: DiseaseVariantIndex,
    ) -> None:
        self.db = reactome_db
        self.disease_index = disease_variant_index
        self.builder = EventContextBuilder(reactome_db)

    def build_normal_contexts(
        self,
        uniprot_ac: str,
        exclude_disease_variant_events: bool = True,
    ) -> List[Dict[str, Any]]:
        reaction_ids = set(
            self.db.get_reactions_for_uniprot(
                uniprot_ac
            )
        )

        if exclude_disease_variant_events:
            reaction_ids.difference_update(
                self.disease_index.get_reactions_for_uniprot(
                    uniprot_ac
                )
            )

        contexts = []

        for reaction_id in sorted(reaction_ids):
            context = self.builder.build(reaction_id)
            context["context_type"] = "NORMAL"
            context["query_uniprot_ac"] = uniprot_ac
            context["disease_variants"] = []
            contexts.append(context)

        return contexts

    def build_disease_contexts(
        self,
        uniprot_ac: str,
    ) -> List[Dict[str, Any]]:
        contexts = []

        for reaction_id in (
            self.disease_index.get_reactions_for_uniprot(
                uniprot_ac
            )
        ):
            context = self.builder.build(reaction_id)

            context["context_type"] = "DISEASE"
            context["query_uniprot_ac"] = uniprot_ac
            context["disease_variants"] = (
                self.disease_index.get_variants_for_reaction(
                    reaction_id,
                    uniprot_ac=uniprot_ac,
                )
            )

            contexts.append(context)

        return contexts

    def build_all_contexts(
        self,
        uniprot_ac: str,
    ) -> Dict[str, List[Dict[str, Any]]]:
        normal = self.build_normal_contexts(
            uniprot_ac
        )
        disease = self.build_disease_contexts(
            uniprot_ac
        )

        return {
            "normal": normal,
            "disease": disease,
            "all": normal + disease,
        }

