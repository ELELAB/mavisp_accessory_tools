#!/usr/bin/env python3
import argparse, csv, os, re, shlex, subprocess, sys, tempfile, time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import requests

UNI_BASE = "https://rest.uniprot.org/uniprotkb/"

# Compartments that permit stable disulfide formation
DISULFIDE_POSITIVE = {
    "endoplasmic reticulum", "er", "er lumen",
    "golgi", "golgi apparatus", "golgi lumen",
    "secreted", "extracellular", "extracellular region",
    "cell membrane (extracellular domain)", 
    "mitochondrial intermembrane space", "ims",
}

# Strongly reducing compartments (disulfides not applicable)
DISULFIDE_NEGATIVE = {
    "cytosol", "cytoplasm", "cytoplasmic", "nucleus", "nuclear", "lysosome", "vacuole",
    "mitochondrial matrix", "peroxisome", "mitochondrion matrix", "chloroplast stroma", "plastid stroma",
}

_norm_re = re.compile(r"[^a-z0-9 ]+")

def norm_loc(s: str) -> str:
    s = s.strip().lower()
    s = s.replace("membrane", " membrane")
    s = _norm_re.sub("", s)
    s = re.sub(r"\s+", " ", s)
    return s

def decide_disulfide_capability(locs: List[str]) -> Tuple[bool, str]:
    if not locs:
        return False, "no localization available"
    pos = {l for l in locs if any(p in l for p in DISULFIDE_POSITIVE)}
    neg = {l for l in locs if any(n in l for n in DISULFIDE_NEGATIVE)}
    if pos:
        return True, f"allows disulfides (evidence: {sorted(list(pos))[:3]})"
    if neg and not pos:
        return False, f"reducing compartments (evidence: {sorted(list(neg))[:3]})"
    return False, "ambiguous/unsupported compartments"

def fetch_uniprot_record(acc: str, timeout: int = 25) -> Optional[Dict]:
    r = requests.get(f"{UNI_BASE}{acc}.json", timeout=timeout)
    return r.json() if r.status_code == 200 else None

def parse_uniprot_localizations(rec: Dict) -> Tuple[List[str], Optional[bool]]:
    locs: List[str] = []
    reviewed = None
    try:
        reviewed = rec.get("entryType", "").lower() == "reviewed"
    except Exception:
        reviewed = None

    for c in rec.get("comments", []):
        if c.get("commentType") == "SUBCELLULAR LOCATION":
            for sl in c.get("subcellularLocations", []):
                for f in (sl.get("location"), sl.get("topology"), sl.get("orientation")):
                    if f and isinstance(f, dict):
                        v = f.get("value")
                        if v:
                            locs.append(norm_loc(v))

    for kw in rec.get("keywords", []):
        v = kw.get("value")
        if v:
            vv = norm_loc(v)
            if vv in {"secreted", "extracellular region"}:
                locs.append(vv)

    for xref in rec.get("uniProtKBCrossReferences", []):
        if xref.get("database") == "GO":
            props = {p["key"]: p.get("value", "") for p in xref.get("properties", []) if isinstance(p, dict) and "key" in p}
            if props.get("GoTerm") == "C":
                desc = props.get("GoTermName") or props.get("Term", "")
                if desc:
                    locs.append(norm_loc(desc))

    seen, uniq = set(), []
    for l in locs:
        if l not in seen:
            uniq.append(l); seen.add(l)
    return uniq, reviewed

def get_uniprot_sequence(rec: Dict) -> Optional[str]:
    seq = rec.get("sequence", {}).get("value")
    return seq if isinstance(seq, str) and seq else None

def run_deeploc2_server(fasta: Path, module_load: str, workdir: Path) -> Tuple[Optional[str], Optional[str]]:
    """
    Runs:
      module load <module_load>
      mkdir -p tmp/model tmp/transformer out/
      export TORCH_HOME=tmp/model ; export TRANSFORMERS_CACHE=tmp/transformer
      deeploc2 -f <fasta> -o out/ -m Accurate -p
      rm -rf tmp
    Returns (top_label, raw_summary)
    """
    bash_script = f"""
set -euo pipefail
module load {module_load}
mkdir -p tmp/model tmp/transformer
mkdir -p out
export TORCH_HOME=$(realpath tmp/model)
export TRANSFORMERS_CACHE=$(realpath tmp/transformer)
deeploc2 -f {shlex.quote(str(fasta))} -o out/ -m Accurate -p
rm -rf tmp
"""
    try:
        proc = subprocess.run(["bash", "-lc", bash_script], cwd=str(workdir),
                              capture_output=True, text=True, timeout=7200)
    except Exception as e:
        return None, f"deeploc2 server run failed: {e}"
    out_dir = workdir / "out"
    return parse_deeploc_outputs(out_dir, proc.stdout)

def parse_deeploc_outputs(out_location: Path, stdout_fallback: str) -> Tuple[Optional[str], Optional[str]]:
    files = []
    if out_location.is_file():
        files = [out_location]
    elif out_location.is_dir():
        files = sorted([p for p in out_location.iterdir() if p.suffix.lower() in {".tsv", ".csv"}],
                       key=lambda p: p.stat().st_mtime, reverse=True)
    if files:
        try:
            path = files[0]
            lines = [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]
            if not lines:
                return None, None
            header = re.split(r"\t|,", lines[0])
            data = re.split(r"\t|,", lines[1]) if len(lines) > 1 else []
            loc_idx = prob_idx = None
            for i, h in enumerate(header):
                hl = h.lower()
                if "loc" in hl and loc_idx is None: loc_idx = i
                if any(p in hl for p in ("prob", "score")) and prob_idx is None: prob_idx = i
            if data and loc_idx is not None:
                top_loc = data[loc_idx].strip()
                preds = []
                for ln in lines[1:6]:
                    parts = re.split(r"\t|,", ln)
                    loc = parts[loc_idx].strip() if loc_idx < len(parts) else ""
                    pr  = parts[prob_idx].strip() if (prob_idx is not None and prob_idx < len(parts)) else ""
                    if loc: preds.append(f"{loc}:{pr}")
                raw = "; ".join(preds) if preds else None
                return top_loc, raw
        except Exception:
            pass
    stdout = (stdout_fallback or "").strip()
    if stdout:
        m = re.search(r"(cytosol|nucleus|nuclear|extracellular|secreted|er|golgi|mitochondrial intermembrane space|mitochondrial matrix|peroxisome)", stdout, re.I)
        top_loc = m.group(1) if m else None
        return top_loc, (stdout[:300] + ("..." if len(stdout) > 300 else ""))
    return None, None

def main():
    ap = argparse.ArgumentParser(
        description="Single-protein localization: UniProt + DeepLoc2 (server flow) -> disulfide plausibility.")
    ap.add_argument("--id", required=True, help="UniProt accession (e.g., Q9UHQ9)")
    ap.add_argument("--fasta", default=None, help="Path to FASTA; if omitted, tries <ID>.fasta or <ID>.fa in CWD; otherwise fetches sequence from UniProt.")
    ap.add_argument("--deeploc", choices=["never", "auto", "always"], default="auto",
                    help="Run DeepLoc2: never, only if UniProt lacks localization (auto), or always.")
    ap.add_argument("--module-load", default="python/3.10/modulefile",
                    help="Module passed to 'module load' for DeepLoc2 server run (default: python/3.10/modulefile)")
    ap.add_argument("-o", "--out", default="localization_summary.csv", help="Output CSV path")
    ap.add_argument("--timeout", type=int, default=25, help="Timeout for UniProt requests (seconds)")
    args = ap.parse_args()

    uid = args.id.strip()
    rec = fetch_uniprot_record(uid, timeout=args.timeout)
    if not rec:
        print(f"[WARN] UniProt fetch failed for {uid}", file=sys.stderr)
    uni_locs, reviewed = (parse_uniprot_localizations(rec) if rec else ([], None))
    seq = get_uniprot_sequence(rec) if rec else None

    basis = "none"
    locs_for_decision: List[str] = []
    deeploc_top: Optional[str] = None
    deeploc_raw: Optional[str] = None

    if uni_locs:
        basis = "uniprot"
        locs_for_decision.extend(uni_locs)

    # Should we run DeepLoc2?
    need_deeploc = (args.deeploc != "never") and (args.deeploc == "always" or not uni_locs)
    fasta_path: Optional[Path] = None

    if need_deeploc:
        if args.fasta:
            fasta_path = Path(args.fasta)
        else:
            for cand in [Path(f"{uid}.fasta"), Path(f"{uid}.fa")]:
                if cand.exists():
                    fasta_path = cand; break
        if fasta_path is None:
            # write temp FASTA if we have the UniProt sequence
            if not seq:
                print(f"[WARN] No sequence available for {uid}; skipping DeepLoc2", file=sys.stderr)
            else:
                tmpdir = Path(tempfile.mkdtemp(prefix=f"deeploc_{uid}_"))
                fasta_path = tmpdir / f"{uid}.fasta"
                with open(fasta_path, "w") as fh:
                    fh.write(f">{uid}\n")
                    for i in range(0, len(seq), 60):
                        fh.write(seq[i:i+60] + "\n")

        if fasta_path and fasta_path.exists():
            dl_top, dl_raw = run_deeploc2_server(fasta_path, args.module_load, Path.cwd())
            deeploc_top, deeploc_raw = dl_top, dl_raw
            if dl_top:
                basis = "both" if basis == "uniprot" else "deeploc"
                locs_for_decision.append(norm_loc(dl_top))

    capable, reason = decide_disulfide_capability(locs_for_decision)

    # Write one-row CSV
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "uniprot_id", "uniprot_localizations", "uniprot_reviewed",
        "deeploc_top_location", "deeploc_all",
        "disulfide_capable", "basis", "reason",
    ]
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerow({
            "uniprot_id": uid,
            "uniprot_localizations": "; ".join(uni_locs) if uni_locs else "",
            "uniprot_reviewed": "" if reviewed is None else str(bool(reviewed)),
            "deeploc_top_location": deeploc_top or "",
            "deeploc_all": deeploc_raw or "",
            "disulfide_capable": str(bool(capable)),
            "basis": basis,
            "reason": reason,
        })
    print(f"Done. Wrote {out_path}")

if __name__ == "__main__":
    main()

