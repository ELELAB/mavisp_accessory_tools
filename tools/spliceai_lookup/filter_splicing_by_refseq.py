#!/usr/bin/env python3
"""Filter spliceAI_lookup CSVs by the NM transcript encoded by an input NP.

Run in data collection, before copying outputs into MAVISp.
Requires pandas and internet access for the NCBI lookup.
"""
import argparse
import ast
import json
import re
import sys
from functools import lru_cache
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
import xml.etree.ElementTree as ET

import pandas as pd


@lru_cache(maxsize=None)
def _transcript_from_protein(protein_id):
    """Resolve one RefSeq NM transcript from the NCBI coded_by field."""
    if not re.fullmatch(r'NP_\d+(?:\.\d+)?', protein_id):
        raise ValueError(
            f'Expected a RefSeq NP accession, got {protein_id!r}'
        )

    query = urlencode({
        'db': 'protein',
        'id': protein_id,
        'rettype': 'gp',
        'retmode': 'xml',
        'tool': 'mavisp_splicing',
    })
    url = (
        'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'
        + query
    )

    try:
        with urlopen(url, timeout=30) as response:
            root = ET.fromstring(response.read())
    except Exception as error:
        raise ValueError(
            f'Cannot retrieve NCBI protein record for {protein_id}: '
            f'{error}'
        ) from error

    transcripts = set()
    for qualifier in root.findall('.//GBQualifier'):
        if qualifier.findtext('GBQualifier_name') == 'coded_by':
            coded_by = qualifier.findtext(
                'GBQualifier_value', default=''
            )
            transcripts.update(
                re.findall(r'NM_\d+(?:\.\d+)?', coded_by)
            )

    if len(transcripts) != 1:
        raise ValueError(
            f'{protein_id}: expected one NM transcript in NCBI coded_by; '
            f'found {sorted(transcripts)}'
        )

    return transcripts.pop()


def filter_outputs(protein_id, input_dir, output_dir):
    input_path = Path(input_dir)
    output = Path(output_dir)

    # Reject reused results before fetching NCBI records.
    if output.resolve() == input_path.resolve():
        raise ValueError(
            'Use a separate output directory to preserve the raw predictions'
        )

    if output.exists():
        if not output.is_dir():
            raise ValueError(f'{output}: output path is not a directory')
        if any(output.iterdir()):
            raise ValueError(
                f'{output}: output directory is not empty; '
                'use a new or empty directory'
            )

    transcript = _transcript_from_protein(protein_id)

    # Match the accession regardless of its version suffix.
    target_accession = transcript.split('.')[0]
    prepared = []

    for filename in ('pangolin_output.csv', 'spliceai_output.csv'):
        source = input_path / filename
        data = pd.read_csv(source, dtype=str)

        required = {
            'ref_seq_id',
            'Mutation',
            'variant_coordinate',
            'Δ_type',
            'Δ_score',
        }
        missing = required - set(data.columns)
        if missing:
            raise ValueError(
                f'{source}: missing columns {sorted(missing)}'
            )

        def matches(value):
            if pd.isna(value):
                return False

            try:
                refs = ast.literal_eval(value)
            except (ValueError, SyntaxError) as error:
                raise ValueError(
                    f'{source}: invalid ref_seq_id list: {value!r}'
                ) from error

            if not isinstance(refs, list):
                raise ValueError(
                    f'{source}: ref_seq_id must contain a list: {value!r}'
                )

            return any(
                isinstance(ref, str)
                and ref.split('.')[0] == target_accession
                for ref in refs
            )

        mask = data['ref_seq_id'].map(matches).astype(bool)
        filtered = data.loc[mask].copy()
        warning = None

        if filtered.empty:
            warning = (
                f'{source}: no rows for {transcript}; '
                'check RefSeq mappings. '
                'The filtered CSV will contain only column headers.'
            )
            print(f'Warning: {warning}', file=sys.stderr)

        prepared.append((filename, filtered, len(data), warning))

    if all(filtered.empty for _, filtered, _, _ in prepared):
        raise ValueError(
            f'No rows for {transcript} in either predictor; '
            'check RefSeq mappings. No outputs were written.'
        )

    output.mkdir(parents=True, exist_ok=True)

    provenance = {
        'refseq_protein_id': protein_id,
        'refseq_transcript_id': transcript,
        'mapping_source': 'NCBI protein record coded_by',
        'files': {},
    }

    for filename, filtered, total, warning in prepared:
        filtered.to_csv(output / filename, index=False)

        details = {
            'input_rows': total,
            'retained_rows': len(filtered),
        }
        if warning is not None:
            details['warning'] = warning

        provenance['files'][filename] = details

        print(
            f'{filename}: {len(filtered)}/{total} '
            f'rows retained for {transcript}'
        )

    (output / 'splicing_provenance.json').write_text(
        json.dumps(provenance, indent=2) + '\n',
        encoding='utf-8',
    )

    return provenance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '-r',
        required=True,
        dest='refseq_protein',
        help='RefSeq NP accession, preferably versioned',
    )
    parser.add_argument(
        '-i',
        required=True,
        dest='input_dir',
        help='Directory containing the two predictor CSVs',
    )
    parser.add_argument(
        '-o',
        required=True,
        dest='output_dir',
        help='New or empty output directory',
    )
    args = parser.parse_args()

    try:
        filter_outputs(
            args.refseq_protein.strip(),
            args.input_dir,
            args.output_dir,
        )
    except Exception as error:
        parser.exit(1, f'Error: {error}\n')


if __name__ == '__main__':
    main()