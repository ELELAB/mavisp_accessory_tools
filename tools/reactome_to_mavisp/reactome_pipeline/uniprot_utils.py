import re
from typing import List, Optional

import pandas as pd
import requests

class UniprotFunctions:
    """Helper functions for UniProt lookups."""

    @staticmethod
    def uniprot_gene_to_uniprot_ac(gene_name: str) -> Optional[str]:
        gene_name = re.sub(r"\(.*", "", str(gene_name)).strip()

        params = {
            "query": f"gene:{gene_name}",
            "fields": "accession,gene_names,protein_name,organism_name",
            "format": "json"
        }

        try:
            response = requests.get(
                "https://rest.uniprot.org/uniprotkb/search",
                params=params,
                timeout=20
            )

            if response.status_code != 200:
                return None

            data = response.json()

            for entry in data.get("results", []):
                if (
                    entry.get("organism", {}).get("scientificName") == "Homo sapiens"
                    and entry.get("entryType") == "UniProtKB reviewed (Swiss-Prot)"
                ):
                    return entry.get("primaryAccession")

            return None

        except requests.exceptions.RequestException:
            return None

    @staticmethod
    def uniprot_id_to_uniprot_ac(uniprot_id: str) -> Optional[str]:
        uniprot_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        uniprot_response = requests.get(uniprot_url)
        if uniprot_response.status_code == 200:
            gene_information = uniprot_response.json()
            if gene_information['organism']['scientificName'] == 'Homo sapiens' and gene_information['entryType'] == "UniProtKB reviewed (Swiss-Prot)":
                return gene_information['primaryAccession']
            return None
        else:
            return None
    @staticmethod
    def uniprot_ac_to_protein_name(uniprot_ac: str) -> Optional[str]:
        uniprot_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_ac}.json"

        try:
            response = requests.get(uniprot_url, timeout=20)

            if response.status_code != 200:
                return None

            data = response.json()

            protein_description = data.get("proteinDescription", {})

            recommended_name = protein_description.get("recommendedName", {})
            full_name = recommended_name.get("fullName", {})

            protein_name = full_name.get("value")

            if protein_name:
                return protein_name

            submission_names = protein_description.get("submissionNames", [])
            if submission_names:
                return submission_names[0].get("fullName", {}).get("value")

            return None

        except requests.exceptions.RequestException:
            return None

    @staticmethod
    def process_gene_names(gene_names: str) -> str:
        if pd.isna(gene_names) or not isinstance(gene_names, str):
            return "None_uniprot_ac"
        genes = gene_names.split('_')
        accessions: List[str] = []
        for gene in genes:
            single_proteins = gene.split("-")
            for protein in single_proteins:
                try:
                    accession = UniprotFunctions.uniprot_gene_to_uniprot_ac(protein.strip().upper())
                    if accession and accession not in accessions:
                        accessions.append(str(accession))
                except Exception:
                    pass
        return "_".join(accessions)