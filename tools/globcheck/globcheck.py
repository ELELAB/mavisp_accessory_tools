import os
import sys
import argparse
import warnings
import pandas as pd
import logging as log
import MDAnalysis as mda
from Bio.PDB import PDBParser
from cath_alphaflow.models.domains import GeneralDomainID, ChoppingPdbResLabel
from cath_alphaflow.commands.measure_globularity import (
    calculate_normed_radius_of_gyration,calculate_packing_density)

# Basic logging configuration
log.basicConfig(level=log.INFO, \
                format='%(levelname)s - %(message)s')

# Thresholds
asphericity_threshold = 0.1
rg_threshold = 0.356
pd_threshold = 10.333

def validate_residue_range(pdb_file, trimming_range):
    """Validates that the provided residue trimming range is 
    within the actual residue range present in the PDB structure.

    Args:
        pdb_file (str): Path to the input PDB
        trimming_range (str): Residue range in format "start-end" ("1-122")
    """
    if not os.path.isfile(pdb_file):
        log.error(f"❌ PDB file '{pdb_file}' does not exist.")
        sys.exit(1)

    start_res, end_res = map(int, trimming_range.split('-'))
    u = mda.Universe(pdb_file)
    residue_ids = [res.resid for res in u.residues]
    min_res, max_res = min(residue_ids), max(residue_ids)

    if not (min_res <= start_res <= max_res) or not (min_res <= end_res <= max_res):
        log.error(f"❌ Error: Trimming range {trimming_range} is outside residue range of the PDB {min_res}-{max_res}")
        sys.exit(1)
    log.info(f"✅ Residue range {trimming_range} is valid.")

def calculate_asphericity(pdb_file):
    """Calculates the asphericity of the structure using MDAnalysis.

    Asphericity is a measure of how non-spherical a structure is,
    with 0 being a perfect sphere.

    Args:
        pdb_file (str):Path to the input PDB

    Returns:
        float: Asphericity value for the structure
    """
    u = mda.Universe(pdb_file)
    asphericity = u.atoms.asphericity(pbc=False)
    log.info(f"✅ Asphericity = {asphericity:.3f}")

    return asphericity

def calculate_globularity_metrics(pdb_file, trimming_range):
    """Calculates globularity metrics for a trimmed domain using CATH AlphaFlow.

    The metrics include:
    - Normalized radius of gyration (Rg)
    - Packing density (PD)

    Args:
        pdb_file (str): Path to the input PDB
        trimming_range (str): Residue range in format "start-end" ("1-122")

    Returns:
        tuple: A tuple (normed_rg, packing_density), both floats
    """
    chopping = ChoppingPdbResLabel.from_str(f"{trimming_range}")

    model_id = os.path.splitext(os.path.basename(pdb_file))[0]
    domain_id = GeneralDomainID(raw_id=model_id, chopping=chopping)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(model_id, pdb_file)

    normed_rg = calculate_normed_radius_of_gyration(domain_id, structure, volume_resolution=5)
    packing_density = calculate_packing_density(domain_id, structure, include_all_atoms=False, distance_cutoff=5)
    log.info(f"✅ Normalized Rg = {normed_rg:.3f}, Packing Density = {packing_density:.3f}")

    return normed_rg, packing_density

def process_results(model_id, asphericity, normed_rg, packing_density, output_file):
    """Processes results to finalise an output
    file in CSV format. Adding column asserting
    Pass/Fail of thresholds.

    Args:
        model_id (str): Identifier for the model
        asphericity (flot): Asphericity value
        normed_rg (float): Normalised rg value
        packing_density (float): Packing Density
        output_file (str): Path to output file 
    """
    result = {
        "model_id": model_id,
        "asphericity": asphericity,
        "normalized_radius_gyration": normed_rg,
        "packing_density": packing_density,}

    df = pd.DataFrame([result])

    df['satisfies_thresholds'] = (
        (df['asphericity'] < asphericity_threshold) &
        (df['normalized_radius_gyration'] < rg_threshold) &
        (df['packing_density'] > pd_threshold))

    df.to_csv(output_file, index=False, float_format="%.3f")
    log.info(f"✅ Final results written to {output_file}")

    row = df.iloc[0]
    status = "✅ PASS" if row['satisfies_thresholds'] else "❌ FAIL"
    log.info(f"{status}: {model_id} | Asp = {row['asphericity']:.3f}, Rg = {row['normalized_radius_gyration']:.3f}, PD = {row['packing_density']:.3f}")

def main():
    warnings.filterwarnings("ignore", category=UserWarning)
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    parser = argparse.ArgumentParser(description="Run globularity metrics using CATH AlphaFlow and MDAnalysis.")
    parser.add_argument("-p", "--pdb", required=True, help="Input trimmed PDB file")
    parser.add_argument("-r", "--range", required=True, help="Residue trimming range (e.g., 160-214)")
    parser.add_argument("-o", "--output", help="Output CSV file")
    args = parser.parse_args()

    model_id = os.path.splitext(os.path.basename(args.pdb))[0]

    output_file = args.output or f"{model_id}_globularity_recap.csv"

    try:
        validate_residue_range(args.pdb, args.range)
        asphericity = calculate_asphericity(args.pdb)
        normed_rg, packing_density = calculate_globularity_metrics(args.pdb, args.range)
        process_results(model_id, asphericity, normed_rg, packing_density, output_file)

    except Exception as e:
        log.error(f"An error occured {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
