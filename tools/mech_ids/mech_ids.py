import os
import re
import argparse
import logging as log
import pandas as pd
from functools import lru_cache
from pathlib import Path


mechanistic_code = {
    r'(?i)destabilizing': 1,
    r'(?i)pathogenic': 1,
    r'(?i)damaging': 1,
    r'(?i)mixed_effects': 1,
    r'(?i)stabilizing': 1,
    r'(?i)neutral': 0,
    r'(?i)benign': 0,
    r'(?i)uncertain': 0,
    r'(?i)ambiguous': 0,
    r'(?i)true': 1,
    r'(?i)false': 0}

effects = {
    "Stability classification": "STABILITY",
    "PTM effect in stability": "PTM STABILITY",
    "AlloSigMA": "LONG RANGE",
    r"Local Int\. classification": "LOCAL INTERACTION",
    "Functional sites": "FUNCTIONAL SITE",
    "EFoldMine - part of early folding region": "EARLY FOLDING",
    "Predicted de-novo disulfide bridge" : "SS-BOND GAIN",
    "Loss of disulfide bridge" : "SS-BOND LOSS",
    "Phosphorylation - loss of function" : "PHOSPHO LOSS",
    "Mutation predicted to add new phosphorylation site": "DENOVO PHOSPHO",
    "Classification of change in folding free energy with phosphorylation" : "DENOVO PHOSPHO",
    "Classification of change in binding free energy with phosphorylation" : "DENOVO PHOSPHO"}

@lru_cache(maxsize=None)
def load_mode_table(protein, mode, db):
    fp = Path(db) / mode / "dataset_tables" / f"{protein}-{mode}.csv"

    if not fp.exists():
        log.warning("%s missing in %s", protein, mode)
        return None

    df = pd.read_csv(fp)

    if "HGVSp" not in df.columns:
        log.warning("HGVSp column missing in %s", fp)
        return None

    df["_variant"] = df["HGVSp"].apply(normalize_variant)
    return df

def normalize_variant(hgvsp):
    if pd.isna(hgvsp):
        return None
    return str(hgvsp).split(":")[-1]

def max_numeric(row, regex):
    values = pd.to_numeric(
        row.filter(regex=regex),
        errors="coerce",)
    return values.max()

def has_value(row, regex, value=1):
    values = pd.to_numeric(
        row.filter(regex=regex),
        errors="coerce",)
    return (values == value).any()

def detect_mechanisms(row, d_threshold):
    # Replace mechanistic codes 
    temp = row.replace(mechanistic_code, regex=True)

    mechs = set()

    for regex, mechanism in effects.items():
        if has_value(temp, regex):
            mechs.add(mechanism)

    if max_numeric(temp, "damaging conformations") > 0:
        mechs.add("STABILITY CONFORMATION DEPENDENT")

    demask = pd.to_numeric(
        row.get("DeMaSk delta fitness"),
        errors="coerce")

    # Only use LoF/GoF if no specific mechanism was found
    if not mechs and pd.notna(demask):
        if demask <= -d_threshold:
            mechs.add("LoF")
        elif demask >= d_threshold:
            mechs.add("GoF")

    return "_".join(sorted(mechs)) if mechs else "NONE"

def process_mode(protein, variant, mode, db, d_threshold):
    df = load_mode_table(protein, mode, db)

    if df is None:
        return "NA"

    hit = df[df["_variant"] == variant]

    if hit.empty:
        log.warning("%s missing for %s in %s", variant, protein, mode)
        return "NA"

    return detect_mechanisms(hit.iloc[0], d_threshold)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help='CSV file containing protein and variants in HGVSp format')
    parser.add_argument("-m", "--mavisp_db", required=True, help='MAVISp DB folder containing simple/ensemble mode directories')
    parser.add_argument("-o", "--output", default='mechanistic_indicators.csv', help='Output file name')
    parser.add_argument("-d", "--demask_threshold", type=float, default=0.25, help='Demask threshold to define LoF/GoF, else default 0.25')
    args = parser.parse_args()

    log.basicConfig(level=log.INFO, format="%(levelname)s - %(message)s")
    
    modes = [d.name for d in Path(args.mavisp_db).iterdir() if d.is_dir() and d.name.endswith("_mode")]

    input = pd.read_csv(args.input)
    out_rows = []

    for r in input.itertuples(index=False):
        row = {"protein": r.protein, 
                "variant": r.variant}
        found = False
        
        for mode in modes:
            res = process_mode(r.protein, r.variant, mode, args.mavisp_db, args.demask_threshold)
            row[mode] = res

            if res != "NA":
                found = True

        if not found:
            log.warning("%s %s missing from all modes", r.protein, r.variant)
        
        out_rows.append(row)

    pd.DataFrame(out_rows).to_csv(args.output, index=False)

if __name__ == "__main__":
    main()
