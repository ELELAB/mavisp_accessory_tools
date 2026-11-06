#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
module load python
python ../MAVISpBenchmark.py -S -y config.yaml -rr exclude_residues.csv -sp
