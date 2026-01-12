#!/usr/bin/env python3
"""
FoldX Version Comparison Tool
Compares protein stability predictions between different FoldX versions.
Generates correlations, scatter plots, confusion matrices, and performance metrics.

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import warnings
warnings.filterwarnings('ignore')

__version__ = "1.0.0"


class FoldXVersionComparison:

    def __init__(self, version1_dir, version2_dir, output_dir,
                 version1_name="FoldX5", version2_name="FoldX5.1",
                 mode='auto', ddg_threshold=1.0, verbose=True):
        self.version1_dir = Path(version1_dir)
        self.version2_dir = Path(version2_dir)
        self.output_dir = Path(output_dir)
        self.version1_name = version1_name
        self.version2_name = version2_name
        self.mode = mode
        self.ddg_threshold = ddg_threshold
        self.verbose = verbose

        self.output_dir.mkdir(parents=True, exist_ok=True)

        if self.mode == 'auto':
            self.mode = self._detect_mode()

        if self.verbose:
            print(f"FoldX Version Comparison Tool v{__version__}")
            print(f"  Mode: {self.mode}")
            print(f"  {self.version1_name}: {self.version1_dir}")
            print(f"  {self.version2_name}: {self.version2_dir}")
            print(f"  Output: {self.output_dir}")
            print(f"  Threshold: {self.ddg_threshold} kcal/mol")

    # File discovery

    def _detect_mode(self):
        if (self.version1_dir / 'simple_mode').exists():
            return 'simple'
        elif (self.version1_dir / 'ensemble_mode').exists():
            return 'ensemble'
        csv_files = list(self.version1_dir.glob('**/*.csv'))
        if csv_files:
            if 'simple_mode' in csv_files[0].name:
                return 'simple'
            elif 'ensemble_mode' in csv_files[0].name:
                return 'ensemble'
        return 'simple'

    def _get_column_mapping(self):
        if self.mode == 'simple':
            return {
                'mutation': 'Mutation',
                'version1_stability': f'Stability ({self.version1_name}, alphafold, kcal/mol)',
                'version2_stability': f'Stability ({self.version2_name}, alphafold, kcal/mol)',
            }
        else:
            return {
                'mutation': 'Mutation',
                'version1_stability': f'Stability ({self.version1_name}, kcal/mol) [md]',
                'version2_stability': f'Stability ({self.version2_name}, kcal/mol) [md]',
            }

    def _find_stability_column(self, df, preferred_names):
        for name in preferred_names:
            if name in df.columns:
                return name
        for col in df.columns:
            col_lower = col.lower()
            if 'stability' in col_lower and 'foldx' in col_lower:
                return col
        return None

    def find_protein_files(self):
        proteins = {}

        if self.mode == 'simple':
            subdir = 'simple_mode/dataset_tables'
        else:
            subdir = 'ensemble_mode/dataset_tables'

        version1_subdir = self.version1_dir / subdir
        version2_subdir = self.version2_dir / subdir

        if not version1_subdir.exists():
            version1_subdir = self.version1_dir
        if not version2_subdir.exists():
            version2_subdir = self.version2_dir

        mode_suffix = f"-{self.mode}_mode.csv"
        version1_files = list(version1_subdir.glob(f"*{mode_suffix}"))

        for v1_file in version1_files:
            protein_name = v1_file.name.replace(mode_suffix, '')
            v2_file = version2_subdir / v1_file.name
            if v2_file.exists():
                proteins[protein_name] = {
                    'version1': v1_file,
                    'version2': v2_file
                }

        if self.verbose:
            print(f"\nFound {len(proteins)} proteins in both versions:")
            for protein in sorted(proteins.keys()):
                print(f"  - {protein}")

        return proteins

    # Load data

    def load_protein_data(self, protein_name, v1_file, v2_file):
        try:
            df1 = pd.read_csv(v1_file)
            df2 = pd.read_csv(v2_file)

            col_map = self._get_column_mapping()

            v1_stability_options = [
                col_map['version1_stability'],
                'Stability (FoldX5, alphafold, kcal/mol)',
                'Stability (FoldX, alphafold, kcal/mol)',
            ]
            v2_stability_options = [
                col_map['version2_stability'],
                'Stability (FoldX5.1, alphafold, kcal/mol)',
                'Stability (FoldX, alphafold, kcal/mol)',
            ]

            v1_col = self._find_stability_column(df1, v1_stability_options)
            v2_col = self._find_stability_column(df2, v2_stability_options)

            if v1_col is None or v2_col is None:
                if self.verbose:
                    print(f"  Warning: Stability column not found for {protein_name}")
                return None

            df1_subset = df1[[col_map['mutation'], v1_col]].copy()
            df2_subset = df2[[col_map['mutation'], v2_col]].copy()

            df1_subset.columns = ['Mutation', f'{self.version1_name}_ddG']
            df2_subset.columns = ['Mutation', f'{self.version2_name}_ddG']

            merged = df1_subset.merge(df2_subset, on='Mutation', how='inner')
            merged = merged.dropna()

            if len(merged) == 0:
                if self.verbose:
                    print(f"  Warning: No overlapping mutations for {protein_name}")
                return None

            merged['Protein'] = protein_name
            return merged

        except Exception as e:
            if self.verbose:
                print(f"  Error loading {protein_name}: {e}")
            return None

    def classify_mutations(self, df):
        df = df.copy()
        v1_col = f'{self.version1_name}_ddG'
        v2_col = f'{self.version2_name}_ddG'

        df[f'{self.version1_name}_class'] = pd.cut(
            df[v1_col],
            bins=[-np.inf, -self.ddg_threshold, self.ddg_threshold, np.inf],
            labels=['Stabilizing', 'Neutral', 'Destabilizing']
        )

        df[f'{self.version2_name}_class'] = pd.cut(
            df[v2_col],
            bins=[-np.inf, -self.ddg_threshold, self.ddg_threshold, np.inf],
            labels=['Stabilizing', 'Neutral', 'Destabilizing']
        )

        return df

    # Analysis

    def calculate_correlation(self, df):
        """Calculate Pearson correlation using numpy"""
        if len(df) < 2:
            return np.nan, np.nan

        v1_col = f'{self.version1_name}_ddG'
        v2_col = f'{self.version2_name}_ddG'

        x = df[v1_col].values
        y = df[v2_col].values

        r = np.corrcoef(x, y)[0, 1]

        # Approximate p-value
        n = len(x)
        t_stat = r * np.sqrt(n - 2) / np.sqrt(1 - r**2)
        from math import exp
        p_value = 2 * (1 - 0.5 * (1 + np.sign(t_stat) * np.sqrt(1 - exp(-2.77 * abs(t_stat) / np.sqrt(n - 2)))))

        return r, p_value

    def _confusion_matrix(self, y_true, y_pred, labels):
        """Confusion matrix implementation"""
        n_labels = len(labels)
        cm = np.zeros((n_labels, n_labels), dtype=int)

        label_to_idx = {label: idx for idx, label in enumerate(labels)}

        for true, pred in zip(y_true, y_pred):
            true_idx = label_to_idx[true]
            pred_idx = label_to_idx[pred]
            cm[true_idx, pred_idx] += 1

        return cm

    def _cohen_kappa(self, y_true, y_pred):
        """Cohen's kappa implementation"""
        y_true = np.array(y_true)
        y_pred = np.array(y_pred)

        # Observed agreement
        po = np.mean(y_true == y_pred)

        # Expected agreement
        n = len(y_true)
        unique_labels = np.unique(np.concatenate([y_true, y_pred]))
        pe = 0
        for label in unique_labels:
            true_count = np.sum(y_true == label)
            pred_count = np.sum(y_pred == label)
            pe += (true_count * pred_count) / (n * n)

        # Cohen's kappa
        if pe == 1:
            return 1.0
        kappa = (po - pe) / (1 - pe)
        return kappa

    def calculate_metrics(self, df):
        v1_class = f'{self.version1_name}_class'
        v2_class = f'{self.version2_name}_class'

        y_true = df[v1_class].values
        y_pred = df[v2_class].values

        metrics = {}
        metrics['accuracy'] = np.mean(y_true == y_pred)
        metrics['cohen_kappa'] = self._cohen_kappa(y_true, y_pred)

        labels = ['Stabilizing', 'Neutral', 'Destabilizing']
        cm = self._confusion_matrix(y_true, y_pred, labels)


        for i, label in enumerate(labels):
            label_lower = label.lower()

            # True positives, false positives, false negatives
            tp = cm[i, i]
            fp = cm[:, i].sum() - tp
            fn = cm[i, :].sum() - tp


            precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0
            f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

            metrics[f'{label_lower}_precision'] = precision
            metrics[f'{label_lower}_recall'] = recall
            metrics[f'{label_lower}_f1'] = f1


        for label in labels:
            label_lower = label.lower()
            metrics[f'n_{label_lower}_{self.version1_name.lower()}'] = np.sum(y_true == label)
            metrics[f'n_{label_lower}_{self.version2_name.lower()}'] = np.sum(y_pred == label)

        return metrics

    # Visualization

    def plot_scatter(self, df, protein_name, correlation, output_file):
        fig, ax = plt.subplots(figsize=(10, 10))

        v1_col = f'{self.version1_name}_ddG'
        v2_col = f'{self.version2_name}_ddG'

        ax.scatter(df[v1_col], df[v2_col], alpha=0.5, s=30, edgecolors='none')

        min_val = min(df[v1_col].min(), df[v2_col].min())
        max_val = max(df[v1_col].max(), df[v2_col].max())
        ax.plot([min_val, max_val], [min_val, max_val], 'r--', 
                alpha=0.5, linewidth=2, label='Perfect agreement')

        ax.axhline(y=self.ddg_threshold, color='gray', linestyle=':', alpha=0.5)
        ax.axhline(y=-self.ddg_threshold, color='gray', linestyle=':', alpha=0.5)
        ax.axvline(x=self.ddg_threshold, color='gray', linestyle=':', alpha=0.5)
        ax.axvline(x=-self.ddg_threshold, color='gray', linestyle=':', alpha=0.5)

        ax.set_xlabel(f'{self.version1_name} ΔΔG (kcal/mol)', fontsize=12)
        ax.set_ylabel(f'{self.version2_name} ΔΔG (kcal/mol)', fontsize=12)
        ax.set_title(f'{protein_name} - {self.version1_name} vs {self.version2_name}\n'
                    f'Pearson r = {correlation:.3f} (n={len(df)})', 
                    fontsize=14, fontweight='bold')

        ax.set_aspect('equal', adjustable='box')
        ax.grid(True, alpha=0.3)
        ax.legend()

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

    def plot_confusion_matrix(self, df, protein_name, output_file):
        v1_class = f'{self.version1_name}_class'
        v2_class = f'{self.version2_name}_class'

        labels = ['Stabilizing', 'Neutral', 'Destabilizing']
        y_true = df[v1_class].values
        y_pred = df[v2_class].values

        cm = self._confusion_matrix(y_true, y_pred, labels)
        cm_normalized = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]

        fig, ax = plt.subplots(figsize=(10, 8))

        im = ax.imshow(cm_normalized, cmap='Blues', vmin=0, vmax=1, aspect='auto')

        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Proportion', fontsize=11)

        ax.set_xticks(np.arange(len(labels)))
        ax.set_yticks(np.arange(len(labels)))
        ax.set_xticklabels(labels)
        ax.set_yticklabels(labels)

        for i in range(len(labels)):
            for j in range(len(labels)):
                text = ax.text(j, i, f'{cm_normalized[i, j]:.2f}',
                             ha="center", va="center",
                             color="black" if cm_normalized[i, j] < 0.5 else "white",
                             fontsize=12)

        ax.set_xlabel(f'{self.version2_name} Classification', fontsize=12)
        ax.set_ylabel(f'{self.version1_name} Classification', fontsize=12)
        ax.set_title(f'{protein_name} - Classification Agreement\n(n={len(df)} mutations)',
                    fontsize=14, fontweight='bold')

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

    def create_overall_plots(self, summary_df):
        if len(summary_df) == 0:
            return

        fig, axes = plt.subplots(2, 2, figsize=(14, 12))

        axes[0, 0].hist(summary_df['pearson_r'], bins=20, edgecolor='black', alpha=0.7)
        axes[0, 0].axvline(summary_df['pearson_r'].mean(), color='red',
                          linestyle='--', linewidth=2, 
                          label=f'Mean: {summary_df["pearson_r"].mean():.3f}')
        axes[0, 0].set_xlabel('Pearson Correlation (r)', fontsize=11)
        axes[0, 0].set_ylabel('Number of Proteins', fontsize=11)
        axes[0, 0].set_title('Distribution of Correlations', fontsize=12, fontweight='bold')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)

        axes[0, 1].scatter(summary_df['pearson_r'], summary_df['accuracy'], alpha=0.6, s=80)
        axes[0, 1].set_xlabel('Pearson Correlation (r)', fontsize=11)
        axes[0, 1].set_ylabel('Classification Accuracy', fontsize=11)
        axes[0, 1].set_title('Accuracy vs Correlation', fontsize=12, fontweight='bold')
        axes[0, 1].grid(True, alpha=0.3)

        axes[1, 0].hist(summary_df['cohen_kappa'], bins=20, edgecolor='black', 
                       alpha=0.7, color='green')
        axes[1, 0].axvline(summary_df['cohen_kappa'].mean(), color='red',
                          linestyle='--', linewidth=2,
                          label=f'Mean: {summary_df["cohen_kappa"].mean():.3f}')
        axes[1, 0].set_xlabel("Cohen's Kappa", fontsize=11)
        axes[1, 0].set_ylabel('Number of Proteins', fontsize=11)
        axes[1, 0].set_title('Distribution of Agreement', fontsize=12, fontweight='bold')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)

        axes[1, 1].scatter(summary_df['n_mutations'], summary_df['pearson_r'], 
                          alpha=0.6, s=80, color='purple')
        axes[1, 1].set_xlabel('Number of Mutations', fontsize=11)
        axes[1, 1].set_ylabel('Pearson Correlation (r)', fontsize=11)
        axes[1, 1].set_title('Sample Size vs Correlation', fontsize=12, fontweight='bold')
        axes[1, 1].grid(True, alpha=0.3)

        plt.tight_layout()
        overall_file = self.output_dir / 'overall_comparison_plots.png'
        plt.savefig(overall_file, dpi=300, bbox_inches='tight')
        plt.close()

        if self.verbose:
            print(f"\nOverall plots saved: {overall_file}")

    # Main

    def compare_protein(self, protein_name, v1_file, v2_file):
        if self.verbose:
            print(f"\nProcessing {protein_name}...")

        df = self.load_protein_data(protein_name, v1_file, v2_file)
        if df is None or len(df) == 0:
            return None

        if self.verbose:
            print(f"  {len(df)} mutations loaded")

        df = self.classify_mutations(df)
        r, p = self.calculate_correlation(df)

        if self.verbose:
            print(f"  Pearson r={r:.3f}, p={p:.2e}")

        protein_output_dir = self.output_dir / protein_name
        protein_output_dir.mkdir(exist_ok=True)

        scatter_file = protein_output_dir / f'{protein_name}_scatter.png'
        self.plot_scatter(df, protein_name, r, scatter_file)

        confusion_file = protein_output_dir / f'{protein_name}_confusion_matrix.png'
        self.plot_confusion_matrix(df, protein_name, confusion_file)

        metrics = self.calculate_metrics(df)

        if self.verbose:
            print(f"  Accuracy: {metrics['accuracy']:.3f}, Kappa: {metrics['cohen_kappa']:.3f}")

        data_file = protein_output_dir / f'{protein_name}_comparison_data.csv'
        df.to_csv(data_file, index=False)

        results = {
            'protein': protein_name,
            'n_mutations': len(df),
            'pearson_r': r,
            'pearson_p': p,
            **metrics
        }

        return results

    def compare_all_proteins(self):
        proteins = self.find_protein_files()

        if len(proteins) == 0:
            print("\nNo matching proteins found!")
            return pd.DataFrame()

        all_results = []
        for protein_name, files in proteins.items():
            result = self.compare_protein(protein_name, files['version1'], files['version2'])
            if result is not None:
                all_results.append(result)

        if len(all_results) == 0:
            print("\nNo successful comparisons!")
            return pd.DataFrame()

        summary_df = pd.DataFrame(all_results)
        summary_file = self.output_dir / 'comparison_summary.csv'
        summary_df.to_csv(summary_file, index=False)

        if self.verbose:
            print(f"\n{'='*70}")
            print(f"Summary: {summary_file}")
            print(f"{'='*70}")
            print(f"Proteins compared: {len(summary_df)}")
            print(f"Mean Pearson r: {summary_df['pearson_r'].mean():.3f} ± {summary_df['pearson_r'].std():.3f}")
            print(f"Mean accuracy: {summary_df['accuracy'].mean():.3f} ± {summary_df['accuracy'].std():.3f}")
            print(f"Mean Cohen's kappa: {summary_df['cohen_kappa'].mean():.3f} ± {summary_df['cohen_kappa'].std():.3f}")

        return summary_df


def main():
    parser = argparse.ArgumentParser(
        description='Compare protein stability predictions between FoldX versions'
    )

    parser.add_argument('--version1-dir', type=str, required=True,
                       help='Directory with first version CSV files')
    parser.add_argument('--version2-dir', type=str, required=True,
                       help='Directory with second version CSV files')
    parser.add_argument('--output-dir', type=str, required=True,
                       help='Output directory')
    parser.add_argument('--version1-name', type=str, default='FoldX5',
                       help='Display name for first version (default: FoldX5)')
    parser.add_argument('--version2-name', type=str, default='FoldX5.1',
                       help='Display name for second version (default: FoldX5.1)')
    parser.add_argument('--mode', type=str, default='auto',
                       choices=['auto', 'simple', 'ensemble'],
                       help='Comparison mode (default: auto)')
    parser.add_argument('--threshold', type=float, default=1.0,
                       help='Stability threshold in kcal/mol (default: 1.0)')
    parser.add_argument('--quiet', action='store_true',
                       help='Suppress progress messages')

    args = parser.parse_args()

    comparison = FoldXVersionComparison(
        version1_dir=args.version1_dir,
        version2_dir=args.version2_dir,
        output_dir=args.output_dir,
        version1_name=args.version1_name,
        version2_name=args.version2_name,
        mode=args.mode,
        ddg_threshold=args.threshold,
        verbose=not args.quiet
    )

    summary_df = comparison.compare_all_proteins()

    if len(summary_df) > 0:
        comparison.create_overall_plots(summary_df)

    print("\n" + "="*70)
    print("Comparison completed")
    print("="*70)


if __name__ == '__main__':
    main()
