#!/usr/bin/env python3
import argparse, csv, sys
import numpy as np, pandas as pd
from sklearn.mixture import GaussianMixture

def args():
    p = argparse.ArgumentParser(
        description="ProteoCast-style 3-class GMM on GEMME scores (neutral/mild/impactful) per protein file."
    )
    p.add_argument("-i","--input",required=True, help="MAVISp CSV for a SINGLE protein")
    p.add_argument("-o","--output",required=True, help="Output CSV (same + GMM columns)")
    p.add_argument("--gemme-col", default="GEMME Score",
                   help='Column with GEMME scores; should match normPred_evolCombi (default: "GEMME Score")')
    p.add_argument("--random-state", type=int, default=42)
    p.add_argument("--max-iter", type=int, default=500)
    p.add_argument("--tol", type=float, default=1e-6)
    return p.parse_args()

def main():
    a = args()
    df = pd.read_csv(a.input, dtype=str, keep_default_na=False)
    if a.gemme_col not in df.columns:
        sys.exit(f'ERROR: column "{a.gemme_col}" not found')

    x = pd.to_numeric(df[a.gemme_col].replace({"": np.nan}), errors="coerce")
    valid = x.notna()
    vals = x[valid].to_numpy().reshape(-1,1)

    # Per ProteoCast: 3-component GMM on GEMME distribution
    if vals.shape[0] < 30:
        sys.exit("ERROR: need >=30 valid scores to fit a 3-component GMM reliably.")

    gmm = GaussianMixture(n_components=3, covariance_type="full",
                          random_state=a.random_state, max_iter=a.max_iter, tol=a.tol)
    gmm.fit(vals)

    # Map components to classes by mean order (more negative GEMME => more impactful)
    means = gmm.means_.flatten()
    order = np.argsort(means)  # ascending: most negative -> most positive
    comp_to_class = {order[0]: "impactful", order[1]: "mild", order[2]: "neutral"}

    post = np.full((len(df),3), np.nan)
    post[valid.values] = gmm.predict_proba(vals)

    # Reorder posteriors into fixed class order
    p_imp = pd.Series(post[:, order[0]], index=df.index)
    p_mild = pd.Series(post[:, order[1]], index=df.index)
    p_neut = pd.Series(post[:, order[2]], index=df.index)

    # Class = argmax posterior
    hard_comp = np.full(len(df), np.nan)
    hard_comp[valid.values] = post[valid.values].argmax(axis=1)
    hard_class = pd.Series(hard_comp, index=df.index).map(comp_to_class)

    # Append columns (names mirror ProteoCast concepts)
    df["ProteoCast_GMM_p_impactful"] = p_imp.round(6)
    df["ProteoCast_GMM_p_mild"]       = p_mild.round(6)
    df["ProteoCast_GMM_p_neutral"]    = p_neut.round(6)
    df["ProteoCast_GMM_class"]        = hard_class
    df["ProteoCast_GMM_means"]        = ";".join(f"{m:.6g}" for m in means)
    df["ProteoCast_GMM_weights"]      = ";".join(f"{w:.6g}" for w in gmm.weights_)

    df.to_csv(a.output, index=False, quoting=csv.QUOTE_MINIMAL)
    print(f"Wrote {a.output}")

if __name__ == "__main__":
    main()

