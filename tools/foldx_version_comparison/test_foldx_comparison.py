#!/usr/bin/env python3
"""
Test script for FoldX Version Comparison Tool

Usage:
    python test_foldx_comparison.py
"""

import sys
import pandas as pd
from pathlib import Path
from foldx_version_comparison import FoldXVersionComparison


FOLDX5_DIR = "/data/user/shared_projects/mavisp_ensemble_sim_length/foldx5.1_evaluation/foldx5_csv_repo/28112025_foldx5_candidates"
FOLDX51_DIR = "/data/user/shared_projects/mavisp_ensemble_sim_length/foldx5.1_evaluation/foldx5.1_csv_repo/1712025_folx5.1_candidates"
TEST_OUTPUT_DIR = "/tmp/test_foldx_comparison"


def test_initialization():
    """Test 1: Tool initialization"""
    print("\n" + "="*70)
    print("TEST 1: Tool Initialization")
    print("="*70)

    try:
        comparison = FoldXVersionComparison(
            version1_dir=FOLDX5_DIR,
            version2_dir=FOLDX51_DIR,
            output_dir=TEST_OUTPUT_DIR,
            mode='simple',
            verbose=False
        )
        print("Tool initialized successfully")
        print(f"  Mode: {comparison.mode}")
        print(f"  Threshold: {comparison.ddg_threshold} kcal/mol")
        return True
    except Exception as e:
        print(f"Initialization failed: {e}")
        return False


def test_file_discovery():
    """Test 2: Finding protein files"""
    print("\n" + "="*70)
    print("TEST 2: File Discovery")
    print("="*70)

    try:
        comparison = FoldXVersionComparison(
            version1_dir=FOLDX5_DIR,
            version2_dir=FOLDX51_DIR,
            output_dir=TEST_OUTPUT_DIR,
            mode='simple',
            verbose=False
        )

        proteins = comparison.find_protein_files()

        if len(proteins) == 0:
            print("No proteins found")
            return False

        print(f"Found {len(proteins)} proteins")
        print(f"  Examples: {', '.join(list(proteins.keys())[:3])}")

        expected_proteins = ['ABI1', 'ACVR1B', 'ADCK2']
        found_expected = [p for p in expected_proteins if p in proteins]
        print(f"  Expected proteins found: {len(found_expected)}/{len(expected_proteins)}")

        return len(proteins) > 0

    except Exception as e:
        print(f"File discovery failed: {e}")
        return False


def test_data_loading():
    """Test 3: Loading and processing protein data"""
    print("\n" + "="*70)
    print("TEST 3: Data Loading")
    print("="*70)

    try:
        comparison = FoldXVersionComparison(
            version1_dir=FOLDX5_DIR,
            version2_dir=FOLDX51_DIR,
            output_dir=TEST_OUTPUT_DIR,
            mode='simple',
            verbose=False
        )

        proteins = comparison.find_protein_files()

        if len(proteins) == 0:
            print("No proteins to test")
            return False

        test_protein = 'ABI1' if 'ABI1' in proteins else list(proteins.keys())[0]
        print(f"  Testing with protein: {test_protein}")

        files = proteins[test_protein]
        df = comparison.load_protein_data(test_protein, files['version1'], files['version2'])

        if df is None:
            print("Failed to load protein data")
            return False

        print(f"Data loaded successfully")
        print(f"  Mutations: {len(df)}")

        required_cols = ['Mutation', 'FoldX5_ddG', 'FoldX5.1_ddG']
        has_required = all(col in df.columns for col in required_cols)

        if not has_required:
            print(f"Missing required columns")
            return False

        print(f"All required columns present")

        null_v1 = df['FoldX5_ddG'].isna().sum()
        null_v2 = df['FoldX5.1_ddG'].isna().sum()

        if null_v1 > 0 or null_v2 > 0:
            print(f"Data contains null values")
            return False

        print(f"No null values in data")
        return True

    except Exception as e:
        print(f"Data loading failed: {e}")
        return False


def test_single_protein():
    """Test 4: Complete single protein comparison"""
    print("\n" + "="*70)
    print("TEST 4: Single Protein Comparison")
    print("="*70)

    try:
        comparison = FoldXVersionComparison(
            version1_dir=FOLDX5_DIR,
            version2_dir=FOLDX51_DIR,
            output_dir=TEST_OUTPUT_DIR,
            mode='simple',
            verbose=False
        )

        proteins = comparison.find_protein_files()
        test_protein = 'ABI1' if 'ABI1' in proteins else list(proteins.keys())[0]
        files = proteins[test_protein]

        print(f"  Testing protein: {test_protein}")

        result = comparison.compare_protein(test_protein, files['version1'], files['version2'])

        if result is None:
            print("Comparison returned None")
            return False

        print(f"Comparison completed")
        print(f"  Mutations: {result['n_mutations']}")
        print(f"  Pearson r: {result['pearson_r']:.3f}")

        output_dir = Path(TEST_OUTPUT_DIR) / test_protein
        scatter_file = output_dir / f'{test_protein}_scatter.png'
        confusion_file = output_dir / f'{test_protein}_confusion_matrix.png'
        data_file = output_dir / f'{test_protein}_comparison_data.csv'

        if not (scatter_file.exists() and confusion_file.exists() and data_file.exists()):
            print("Some output files missing")
            return False

        print(f"All output files created")
        return True

    except Exception as e:
        print(f"Single protein comparison failed: {e}")
        return False


def main():
    """Run all tests"""
    print("="*70)
    print("FoldX Version Comparison Tool - Test Suite")
    print("="*70)

    foldx5_exists = Path(FOLDX5_DIR).exists()
    foldx51_exists = Path(FOLDX51_DIR).exists()

    if not foldx5_exists:
        print(f" FoldX5 directory not found: {FOLDX5_DIR}")
        return 1

    if not foldx51_exists:
        print(f"FoldX5.1 directory not found: {FOLDX51_DIR}")
        return 1

    print(f"Both directories found")

    tests = [
        ("Initialization", test_initialization),
        ("File Discovery", test_file_discovery),
        ("Data Loading", test_data_loading),
        ("Single Protein", test_single_protein),
    ]

    results = []
    for test_name, test_func in tests:
        try:
            success = test_func()
            results.append((test_name, success))
        except Exception as e:
            print(f"Test '{test_name}' crashed: {e}")
            results.append((test_name, False))

    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)

    for test_name, success in results:
        status = "PASS" if success else "FAIL"
        print(f"{status}: {test_name}")

    passed = sum(1 for _, success in results if success)
    total = len(results)

    print(f"\nPassed: {passed}/{total}")

    if passed == total:
        print("\n All tests passed!")
        print(f"\nTest output: {TEST_OUTPUT_DIR}")
        return 0
    else:
        print(f"\n {total - passed} test(s) failed")
        return 1


if __name__ == '__main__':
    sys.exit(main())
