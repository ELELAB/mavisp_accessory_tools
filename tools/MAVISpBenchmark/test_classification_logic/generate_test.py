"""Generate deterministic synthetic data for priority, weighted, voting, and common-MCC tests."""
from pathlib import Path
import itertools
import pandas as pd

ST = "Stability classification, AFDB, (RaSP, FoldX)"
LOC = "Local Int. classification (PMS2_AFmulti)"
EXP = "Experimental data (DHFR-PCA yeast, Stability (MAVEdb))"
CLASS = "Experimental data classification (DHFR-PCA yeast, Stability (MAVEdb))"
AA = "ACDEFGHIKLMNPQRSTVWY"
rows = []
for i, (st, loc, target, demask, gemme) in enumerate(itertools.product(("Damaging", "Neutral"), repeat=5)):
    rows.append({"Mutation": f"{AA[i % 20]}{i + 1}{AA[(i + 3) % 20]}", ST: st, LOC: loc,
                 EXP: -0.6 if target == "Damaging" else 0.3, CLASS: target,
                 "DeMaSk delta fitness": -0.4 if demask == "Damaging" else -0.1,
                 "GEMME Score": -4.0 if gemme == "Damaging" else -1.0,
                 "Relative Side Chain Solvent Accessibility in wild-type": 20.0})
for j, missing in enumerate(("DeMaSk delta fitness", "GEMME Score")):
    row = {"Mutation": f"{AA[j]}{101+j}{AA[j+5]}", ST: "Uncertain", LOC: "Uncertain",
           EXP: -0.5, CLASS: "Damaging", "DeMaSk delta fitness": -0.4, "GEMME Score": -4.0,
           "Relative Side Chain Solvent Accessibility in wild-type": 20.0}
    row[missing] = float("nan")
    rows.append(row)

pd.DataFrame(rows).to_csv(Path(__file__).with_name("test-simple_mode.csv"), index=False)
print(f"Generated {len(rows)} deterministic test variants (32 complete, 2 with missing VEP score).")
