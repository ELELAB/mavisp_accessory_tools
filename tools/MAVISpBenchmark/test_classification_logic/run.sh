#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
module load python
python generate_test.py
python ../MAVISpBenchmark.py -S -y config.yaml --no-plots

