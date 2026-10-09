#!/usr/bin/env python3
"""Filter spliceAI_lookup CSVs by the NM transcript encoded by an input NP.

Run in data collection, before copying outputs into MAVISp. Requires pandas.
"""
import argparse
import ast
import json
import re
from pathlib import Path
from functools import lru_cache
from urllib.parse import urlencode
from urllib.request import urlopen
import xml.etree.ElementTree as ET
import pandas as pd

@lru_cache(maxsize=None)
def _transcript_from_protein(protein_id):
    """Resolve the RefSeq transcript from the protein's NCBI coded_by field.

    Cache within this process so simple/ensemble modes reuse the lookup.
    Do not guess when the record has no unique NM transcript.
    """
    if not re.fullmatch(r'NP_\d+(?:\.\d+)?', protein_id):
        raise ValueError(f'Expected a RefSeq NP accession, got {protein_id!r}')
    query = urlencode({
        'db': 'protein', 'id': protein_id,
        'rettype': 'gp', 'retmode': 'xml', 'tool': 'mavisp_splicing',
    })
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + query
    try:
        with urlopen(url, timeout=30) as response:
            root = ET.fromstring(response.read())
    except Exception as error:
        raise ValueError(
            f'Cannot retrieve NCBI protein record for {protein_id}: {error}'
        ) from error
    transcripts = set()
    for qualifier in root.findall('.//GBQualifier'):
        if qualifier.findtext('GBQualifier_name') == 'coded_by':
            coded_by = qualifier.findtext('GBQualifier_value', default='')
            transcripts.update(re.findall(r'NM_\d+(?:\.\d+)?', coded_by))
    if len(transcripts) != 1:
        raise ValueError(
            f'{protein_id}: expected one NM transcript in NCBI coded_by; '
            f'found {sorted(transcripts)}'
        )
    transcript = transcripts.pop()
    return transcript



def filter_outputs(protein_id, input_dir, output_dir):
    transcript = _transcript_from_protein(protein_id)
    prepared = []
    for filename in ('pangolin_output.csv', 'spliceai_output.csv'):
        source = Path(input_dir) / filename
        data = pd.read_csv(source, dtype=str)
        required = {'ref_seq_id', 'Mutation', 'variant_coordinate', 'Δ_type', 'Δ_score'}
        if not required.issubset(data.columns):
            raise ValueError(f'{source}: missing columns {sorted(required - set(data.columns))}')
        def matches(value):
            if pd.isna(value):
                return False
            refs = ast.literal_eval(value)
            if not isinstance(refs, list):
                raise ValueError(f'{source}: ref_seq_id must contain a list: {value!r}')
            return any(
                ref.split('.')[0] == transcript.split('.')[0]
                for ref in refs
                if isinstance(ref, str)
            )
        filtered = data.loc[data['ref_seq_id'].map(matches)].copy()
        if filtered.empty:
            raise ValueError(f'{source}: no rows for {transcript}; check RefSeq mappings')
        prepared.append((filename, filtered, len(data)))
    output = Path(output_dir)
    if output.resolve() == Path(input_dir).resolve():
        raise ValueError('Use a separate output directory to preserve the raw predictions')
    output.mkdir(parents=True, exist_ok=True)
    provenance = {'refseq_protein_id': protein_id, 'refseq_transcript_id': transcript,
                  'mapping_source': 'NCBI protein record coded_by', 'files': {}}
    for filename, filtered, total in prepared:
        filtered.to_csv(output / filename, index=False)
        provenance['files'][filename] = {'input_rows': total, 'retained_rows': len(filtered)}
        print(f'{filename}: {len(filtered)}/{total} rows retained for {transcript}')
    (output / 'splicing_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    return provenance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-r', required=True, dest="refseq_protein", help='NP accession, preferably versioned')
    parser.add_argument('-i', required=True, dest='input_dir')
    parser.add_argument('-o', required=True, dest='output_dir')
    args = parser.parse_args()
    try:
        filter_outputs(args.refseq_protein.strip(), args.input_dir, args.output_dir)
    except Exception as error:
        parser.exit(1, f'Error: {error}\n')


if __name__ == '__main__':
    main()

