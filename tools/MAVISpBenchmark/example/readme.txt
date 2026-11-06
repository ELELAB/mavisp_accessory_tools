Updated MLH1 Simple mode example (2026).
Input: MLH1-simple_mode.csv (11,324 variants), with current MAVISp column names.
Run all analyses and plots: bash run.sh
For a fast metrics-only check: python ../MAVISpBenchmark.py -S -y config.yaml --no-plots -rr exclude_residues.csv
The YAML specifies GEMME/DeMaSk thresholds and includes direct single-column comparisons through explicit priority strategies as well as MAVISp voting/weighted/priority combinations.
Outputs: output/simple/common_variants.csv, performance_summary.csv, performance_comparison.png and individual <comparison>/<strategy>/ folders.
