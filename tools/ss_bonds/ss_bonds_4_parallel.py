#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MAVISp Disulfide Module — WT Analysis + De Novo Prediction (parallel + verbose)
===============================================================================

Features
- WT mode: detect native disulfides using strict/loose geometry; optional coarse rotamer scan.
- De novo mode: scan variant folders; implicit-Cys folder names (e.g., IA1124) supported.
- Parallelization: --jobs for concurrent per-variant scanning.
- Ignores WT PDB: only reads mutant files named like *_X_Y.pdb (e.g., ..._1_0.pdb).
- Majority vote: uses ONLY models where the specified site is actually a cysteine.

Examples
WT mode:
  python disulfide_module.py wt WT.pdb \
    --native-mode LOOSE --ss-loose-min 1.8 --ss-loose-max 3.5 --no-angle \
    --rotamer-scan --rot-angles "-60,60,180" --rot-threshold 10.0 \
    --dump-pairs all_pairs.csv --out-csv native_disulfides.csv \
    --variants-file variants_noheader.txt --chain A \
    --out-mutlist mutlist_disulfide_disruptions.tsv

De novo (parallel, verbose):
  python disulfide_module.py denovo /path/to/variants_root \
    --out denovo_disulfide_summary.csv --details denovo_disulfide_details.csv \
    --majority 0.5 --ss-min 1.95 --ss-max 2.25 --pot-min 2.25 --pot-max 3.5 \
    --jobs 8 --verbose
"""

import os
import re
import math
import time
import argparse
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser

# ---------------------------- geometry helpers ---------------------------- #

@dataclass
class AtomRef:
    chain: str
    resname: str
    resnum: int
    icode: str
    atomname: str
    coord: np.ndarray

def dist(a, b):
    return float(np.linalg.norm(a - b))

def angle(a, b, c):
    ba, bc = a - b, c - b
    cosang = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc) + 1e-12)
    return float(np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0))))

def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = -1.0*(p1-p0), p2-p1, p3-p2
    b1 /= np.linalg.norm(b1) + 1e-12
    v, w = b0 - np.dot(b0,b1)*b1, b2 - np.dot(b2,b1)*b1
    x, y = np.dot(v,w), np.dot(np.cross(b1,v),w)
    return float(np.degrees(np.arctan2(y,x)))

def _rotate_point_around_axis(point: np.ndarray, axis_point: np.ndarray, axis_dir: np.ndarray, angle_deg: float) -> np.ndarray:
    """Rigid rotation of a point around an axis defined by (axis_point, axis_dir)."""
    p = point - axis_point
    u = axis_dir / (np.linalg.norm(axis_dir) + 1e-12)
    theta = math.radians(angle_deg)
    p_rot = p * math.cos(theta) + np.cross(u, p) * math.sin(theta) + u * np.dot(u, p) * (1 - math.cos(theta))
    return axis_point + p_rot

# ---------------------------- PDB utilities ------------------------------ #

def read_cysteine_sg_atoms(pdb_path: str):
    """Return list of SG atoms and a residue map for cysteines (incl. CYX/CYM/CYD/CSE)."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", pdb_path)
    sg_list, resmap = [], {}
    for model in structure:
        for chain in model:
            ch = chain.id.strip() or 'A'
            for residue in chain:
                resname = residue.get_resname().strip()
                hetflag, resseq, icode = residue.id
                if hetflag != " ":
                    continue
                if resname not in ("CYS", "CYX", "CYD", "CYM", "CSE") and not resname.startswith("CY"):
                    continue
                atomdict = {}
                for an in ("SG", "CB", "CA"):
                    if an in residue:
                        a = residue[an]
                        atomdict[an] = AtomRef(
                            ch, resname, int(resseq), icode.strip() if isinstance(icode, str) else "", an, a.coord.copy()
                        )
                if "SG" in atomdict:
                    sg_list.append(atomdict["SG"])
                    resmap[(ch, int(resseq), atomdict["SG"].icode)] = atomdict
    return sg_list, resmap

# -------------------------- Disulfide detection -------------------------- #

@dataclass
class DisulfideCriteria:
    ss_min: float = 1.95
    ss_max: float = 2.25
    cb_sg_sg_angle_min: float = 85.0
    cb_sg_sg_angle_max: float = 130.0
    sg_sg_cb_angle_min: float = 85.0
    sg_sg_cb_angle_max: float = 130.0
    chi_ss_abs: float = 120.0

@dataclass
class LooseCriteria:
    ss_min: float = 1.8
    ss_max: float = 3.5
    use_angles: bool = False

# ------------------------------ WT Analyzer ------------------------------ #

class WTDisulfideAnalyzer:
    def __init__(
        self,
        strict: DisulfideCriteria,
        loose: LooseCriteria,
        rotamer_scan: bool = False,
        rot_angles: Optional[List[float]] = None,
        rot_threshold: float = 10.0
    ):
        self.strict = strict
        self.loose = loose
        self.rotamer_scan = rotamer_scan
        self.rot_angles = rot_angles or [-60.0, 60.0, 180.0]
        self.rot_threshold = rot_threshold

    def _min_rotamer_distance(self, a_map: Dict[str, AtomRef], b_map: Dict[str, AtomRef]) -> Tuple[float, Tuple[float, float]]:
        """Coarse χ1-like scan to check reachability between two Cys SG atoms."""
        req = {"CA", "CB", "SG"}
        if not (req.issubset(a_map.keys()) and req.issubset(b_map.keys())):
            return float('inf'), (0.0, 0.0)
        CAa, CBa, SGa = a_map['CA'].coord, a_map['CB'].coord, a_map['SG'].coord
        CAb, CBb, SGb = b_map['CA'].coord, b_map['CB'].coord, b_map['SG'].coord
        axisA = CBa - CAa
        axisB = CBb - CAb
        best_d = float('inf')
        best_ang = (0.0, 0.0)
        for angA in self.rot_angles:
            SGa_rot = _rotate_point_around_axis(SGa, CBa, axisA, angA)
            for angB in self.rot_angles:
                SGb_rot = _rotate_point_around_axis(SGb, CBb, axisB, angB)
                d = dist(SGa_rot, SGb_rot)
                if d < best_d:
                    best_d = d
                    best_ang = (angA, angB)
        return best_d, best_ang

    def analyze(self, pdb_path: str, native_mode: str = 'STRICT', dump_pairs: Optional[str] = None, debug: bool = False) -> pd.DataFrame:
        native_mode = native_mode.upper()
        rows = []
        _, resmap = read_cysteine_sg_atoms(pdb_path)
        keys = list(resmap.keys())
        if debug:
            print(f"Found {len(keys)} cysteines with SG; chains: {sorted(set(k[0] for k in keys))}")
        dump_rows = []
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                akey, bkey = keys[i], keys[j]
                a_map, b_map = resmap[akey], resmap[bkey]
                hasSG = ('SG' in a_map and 'SG' in b_map)
                d0 = dist(a_map['SG'].coord, b_map['SG'].coord) if hasSG else np.nan
                min_rot_d, min_rot_ang = (np.nan, (0.0, 0.0))
                rot_reachable = False

                if self.rotamer_scan and hasSG and not np.isnan(d0) and d0 <= self.rot_threshold:
                    min_rot_d, min_rot_ang = self._min_rotamer_distance(a_map, b_map)
                    if self.loose.ss_min <= min_rot_d <= self.loose.ss_max:
                        rot_reachable = True

                if native_mode == 'STRICT':
                    ok, met = self._strict_pair(a_map, b_map)
                else:
                    ok, met = self._loose_pair(a_map, b_map)

                if ok:
                    row = {
                        "chainA": akey[0], "resA": akey[1], "icodeA": akey[2],
                        "chainB": bkey[0], "resB": bkey[1], "icodeB": bkey[2],
                        "ok": True
                    }
                    row.update(met)
                    row.update({"min_rot_d_ss": float(min_rot_d), "rotamer_reachable": bool(rot_reachable)})
                    rows.append(row)

                dump_rows.append({
                    'chainA': akey[0], 'resA': akey[1], 'chainB': bkey[0], 'resB': bkey[1],
                    'd_ss': float(d0) if not np.isnan(d0) else np.nan,
                    'min_rot_d_ss': float(min_rot_d) if not np.isnan(min_rot_d) else np.nan,
                    'min_rot_ang_A': min_rot_ang[0] if not np.isnan(min_rot_d) else np.nan,
                    'min_rot_ang_B': min_rot_ang[1] if not np.isnan(min_rot_d) else np.nan,
                })
        if dump_pairs:
            pd.DataFrame(dump_rows).to_csv(dump_pairs, index=False)
        return pd.DataFrame(rows)

    def _strict_pair(self, a_map: Dict[str, AtomRef], b_map: Dict[str, AtomRef]):
        if "SG" not in a_map or "SG" not in b_map:
            return False, {}
        sgA, sgB = a_map['SG'].coord, b_map['SG'].coord
        d_ss = dist(sgA, sgB)
        met = {"d_ss": d_ss}
        ok = (self.strict.ss_min <= d_ss <= self.strict.ss_max)
        if ok and ("CB" in a_map and "CB" in b_map):
            a1 = angle(a_map['CB'].coord, a_map['SG'].coord, b_map['SG'].coord)
            a2 = angle(a_map['SG'].coord, b_map['SG'].coord, b_map['CB'].coord)
            met.update({"angle_cb_sg_sg": a1, "angle_sg_sg_cb": a2})
            ok = (self.strict.cb_sg_sg_angle_min <= a1 <= self.strict.cb_sg_sg_angle_max) and \
                 (self.strict.sg_sg_cb_angle_min <= a2 <= self.strict.sg_sg_cb_angle_max)
            if ok:
                chi = dihedral(a_map['CB'].coord, a_map['SG'].coord, b_map['SG'].coord, b_map['CB'].coord)
                met['chi_ss'] = chi
                ok = abs(chi) <= self.strict.chi_ss_abs
        return ok, met

    def _loose_pair(self, a_map: Dict[str, AtomRef], b_map: Dict[str, AtomRef]):
        if "SG" not in a_map or "SG" not in b_map:
            return False, {}
        d = dist(a_map['SG'].coord, b_map['SG'].coord)
        ok = self.loose.ss_min <= d <= self.loose.ss_max
        met = {"d_ss": d}
        return ok, met

    @staticmethod
    def export_mutlist_simple(df_native: pd.DataFrame, variants_file: str, out_tsv: str, default_chain: str = 'A', native_sites_file: Optional[str] = None):
        """From native disulfides + optional known sites, list Cys->nonCys as potential disruptions."""
        posset = set()
        if df_native is not None and not df_native.empty:
            for _, r in df_native.iterrows():
                posset.add((str(r.chainA), int(r.resA)))
                posset.add((str(r.chainB), int(r.resB)))
        if native_sites_file:
            try:
                with open(native_sites_file, 'r') as fh:
                    for line in fh:
                        s = line.strip()
                        if not s or s.startswith('#'):
                            continue
                        pos = int(s)
                        posset.add((default_chain, pos))
            except Exception:
                pass
        if not posset:
            pd.DataFrame(columns=["variant_id", "reason"]).to_csv(out_tsv, sep='\t', index=False)
            return

        vdf = pd.read_csv(variants_file, header=None, names=['tok'])
        rows = []
        for _, v in vdf.iterrows():
            s = str(v['tok']).strip()
            m = re.match(r'^([A-Z])(\d+)([A-Z])$', s)
            if not m:
                continue
            wt, pos, mut = m.groups()
            pos = int(pos)
            ch = default_chain
            if wt == 'C' and (ch, pos) in posset and mut != 'C':
                rows.append({
                    'variant_id': f'{wt}{ch}{pos}{mut}',
                    'chain': ch,
                    'pos': pos,
                    'wt': wt,
                    'mut': mut,
                    'reason': 'disulfide_disruption(native_Cys->nonCys)'
                })
        pd.DataFrame(rows).to_csv(out_tsv, sep='\t', index=False)

# ------------------------------ De Novo (worker) ------------------------- #

@dataclass
class PotentialCriteria:
    ss_min: float = 2.25
    ss_max: float = 3.5
    angle_min: float = 70.0
    angle_max: float = 140.0
    chi_abs: float = 150.0

class DeNovoScanner:
    def __init__(
        self,
        strict: DisulfideCriteria,
        loose: PotentialCriteria,
        majority: float = 0.5,
        model_id_regex: str = r"_(?P<mid>\d+)\.pdb$",
        model_dir_regex: str = r"model_(?P<mid>\d+)$",
        skip_non_cys: bool = True
    ):
        self.strict = strict
        self.loose = loose
        self.majority = majority
        self._re_file_id = re.compile(model_id_regex, re.IGNORECASE)
        self._re_dir_id = re.compile(model_dir_regex, re.IGNORECASE)
        self.skip_non_cys = skip_non_cys

    @staticmethod
    def parse_folder_mutation(folder_name: str) -> Optional[Dict[str, str]]:
        """
        Accepts:
         - IA1124  (WT=I, chain=A, pos=1124) → implicit mut='C'
         - VA260_C (WT=V, chain=A, pos=260, mut=C) and similar older styles
        """
        # IA1124 implicit C
        m = re.match(r"^(?P<wt>[A-Z])(?P<chain>[A-Za-z])(?P<pos>\d+)$", folder_name)
        if m:
            d = m.groupdict()
            return {"wt": d["wt"], "chain": d["chain"], "pos": d["pos"], "mut": "C"}

        # Older styles
        pats = [
            re.compile(r"^(?P<wt>[A-Z])(?P<chain>[A-Za-z])(?P<pos>\d+)_?(?P<mut>[A-Z])$"),
            re.compile(r"^(?P<chain>[A-Za-z])(?P<pos>\d+)(?P<wt>[A-Z])->?(?P<mut>[A-Z])$"),
            re.compile(r"^(?P<wt>[A-Z])_(?P<chain>[A-Za-z])_(?P<pos>\d+)_?(?P<mut>[A-Z])$"),
        ]
        for p in pats:
            mm = p.match(folder_name)
            if mm:
                d = mm.groupdict()
                return {"wt": d.get("wt"), "chain": d.get("chain"), "pos": d.get("pos"), "mut": d.get("mut")}
        return None

    def _model_id_from_path(self, path: str) -> str:
        m = self._re_file_id.search(os.path.basename(path))
        if m:
            return m.group('mid')
        parent = os.path.basename(os.path.dirname(path))
        m2 = self._re_dir_id.match(parent)
        if m2:
            return m2.group('mid')
        return os.path.splitext(os.path.basename(path))[0]

    def assess_pair(self, a_map, b_map):
        """Return (formed_strict, metrics_strict, potential_loose, metrics_loose)."""
        formed, met_strict = WTDisulfideAnalyzer._strict_pair(self, a_map, b_map)
        potential = False
        met_loose = {}
        if 'SG' in a_map and 'SG' in b_map:
            d = dist(a_map['SG'].coord, b_map['SG'].coord)
            if self.loose.ss_min <= d <= self.loose.ss_max:
                potential = True
            met_loose['d_ss'] = d
        return formed, met_strict, potential, met_loose

    def scan_variant_folder(self, folder: str) -> Dict:
        meta = self.parse_folder_mutation(os.path.basename(folder)) or {}
        if not meta:
            return {
                "variant_folder": os.path.basename(folder),
                "n_models_total": 0, "n_models": 0, "n_positive": 0,
                "fraction_positive": 0.0, "classification": "invalid_folder", "details": []
            }
        if self.skip_non_cys and meta.get("mut") != "C":
            return {
                "variant_folder": os.path.basename(folder), **meta,
                "n_models_total": 0, "n_models": 0, "n_positive": 0,
                "fraction_positive": 0.0, "classification": "not_applicable", "details": []
            }

        # Collect ONLY mutant models like ..._X_Y.pdb (skip bare WT *.pdb)
        models = []
        include_re = re.compile(r'.*_\d+_\d+\.pdb$', re.IGNORECASE)
        for root, _, files in os.walk(folder):
            for fn in files:
                if fn.lower().endswith('.pdb') and include_re.match(fn):
                    models.append(os.path.join(root, fn))
        models.sort()

        n_total = len(models)  # all mutant PDBs found by naming
        n = 0                  # models where the specified site is actually Cys
        n_pos = 0
        details = []

        for mpath in models:
            _, resmap = read_cysteine_sg_atoms(mpath)
            keys = list(resmap.keys())

            # Find the mutant site key (chain, position)
            mut_key = None
            if meta.get('pos', '').isdigit():
                posi = int(meta['pos'])
                chain = meta.get('chain', 'A')
                for k in keys:
                    if k[0] == chain and k[1] == posi:
                        mut_key = k
                        break

            if mut_key not in keys:
                # mutant site is not Cys in this model → exclude from denominator
                details.append({
                    "model": self._model_id_from_path(mpath),
                    "positive": False,
                    "formed_pairs": 0,
                    "potential_pairs": 0,
                    "mutant_present": False,
                    "d_ss": np.nan,
                })
                continue

            # Only test pairs involving the mutant Cys
            pair_iter = [(mut_key, k) for k in keys if k != mut_key]
            n += 1  # this model is counted

            model_flag = False
            formed_count = 0
            potential_count = 0
            best = {"d_ss": 1e9}
            for akey, bkey in pair_iter:
                a_map, b_map = resmap[akey], resmap[bkey]
                formed, met_strict, potential, met_loose = self.assess_pair(a_map, b_map)
                met = met_strict if formed else met_loose
                if met and met.get("d_ss", 1e9) < best.get("d_ss", 1e9):
                    best = {
                        **met,
                        "pair": f"{akey[0]}:{akey[1]}-{bkey[0]}:{bkey[1]}",
                        "formed": formed, "potential": potential
                    }
                if formed:
                    formed_count += 1
                if potential:
                    potential_count += 1
                if formed or potential:
                    model_flag = True

            if model_flag:
                n_pos += 1

            details.append({
                "model": self._model_id_from_path(mpath),
                "positive": model_flag,
                "formed_pairs": formed_count,
                "potential_pairs": potential_count,
                "mutant_present": True,
                **best,
            })

        frac = (n_pos / n) if n else 0.0
        label = "de_novo_disulfide" if (n > 0 and frac >= self.majority) else ("neutral" if n > 0 else "not_applicable")
        return {
            "variant_folder": os.path.basename(folder), **meta,
            "n_models_total": n_total,
            "n_models": n,
            "n_positive": n_pos,
            "fraction_positive": round(frac, 3),
            "classification": label,
            "details": details,
        }

# ------------------------------ CLI utilities --------------------------- #

def _human_time(sec: float) -> str:
    m, s = divmod(int(sec), 60)
    h, m = divmod(m, 60)
    if h:
        return f"{h}h{m}m{s}s"
    if m:
        return f"{m}m{s}s"
    return f"{s}s"

# ------------------------------ CLIs ------------------------------------- #

def cli_wt():
    ap = argparse.ArgumentParser(description='Detect native disulfide bonds (AF-friendly), optional rotamer scan; export mutlist; optional pair dump.')
    ap.add_argument('pdb', help='WT PDB/AF model file')
    ap.add_argument('--variants-file', default=None, help='Plain list: A260V, C360V ... (no header)')
    ap.add_argument('--chain', default='A', help='Assumed chain for simple list tokens')
    ap.add_argument('--native-sites', default=None, help='Optional file with native cysteine positions (one integer per line); chain assumed from --chain')
    ap.add_argument('--out-csv', default='native_disulfides.csv')
    ap.add_argument('--out-mutlist', default='mutlist_disulfide_disruptions.tsv')
    ap.add_argument('--native-mode', choices=['STRICT', 'LOOSE'], default='LOOSE')
    ap.add_argument('--ss-min', type=float, default=1.95)
    ap.add_argument('--ss-max', type=float, default=2.25)
    ap.add_argument('--ss-loose-min', type=float, default=1.8)
    ap.add_argument('--ss-loose-max', type=float, default=3.5)
    ap.add_argument('--no-angle', action='store_true', help='In LOOSE mode, skip angle constraints')
    # Rotamer scan controls
    ap.add_argument('--rotamer-scan', action='store_true', help='Enable coarse rotamer scan for close pairs')
    ap.add_argument('--rot-angles', type=str, default='-60,60,180', metavar='ANGLES', help='Comma-separated χ1-like angles (deg) — quote the list')
    ap.add_argument('--rot-threshold', type=float, default=10.0, help='Only scan pairs with initial SG–SG ≤ this (Å)')
    ap.add_argument('--dump-pairs', default=None, help='Write all Cys–Cys pairs to this CSV (with distances & rotamer min distance)')
    ap.add_argument('--debug', action='store_true')
    args = ap.parse_args()

    strict = DisulfideCriteria(ss_min=args.ss_min, ss_max=args.ss_max)
    loose = LooseCriteria(ss_min=args.ss_loose_min, ss_max=args.ss_loose_max, use_angles=(not args.no_angle))

    try:
        rot_angles = [float(x) for x in args.rot_angles.split(',')] if args.rotamer_scan else None
    except Exception:
        rot_angles = [-60.0, 60.0, 180.0]

    analyzer = WTDisulfideAnalyzer(strict, loose, rotamer_scan=args.rotamer_scan, rot_angles=rot_angles, rot_threshold=args.rot_threshold)
    df = analyzer.analyze(args.pdb, native_mode=args.native_mode, dump_pairs=args.dump_pairs, debug=args.debug)
    df.to_csv(args.out_csv, index=False)

    if args.debug:
        print(df if not df.empty else 'No native disulfides detected with current settings.')

    if args.variants_file is not None:
        WTDisulfideAnalyzer.export_mutlist_simple(
            df, args.variants_file, args.out_mutlist, default_chain=args.chain, native_sites_file=args.native_sites
        )

def _worker_scan_folder(payload: Dict) -> Dict:
    """Subprocess worker wrapper to scan one variant folder."""
    scanner = DeNovoScanner(
        strict=DisulfideCriteria(ss_min=payload['ss_min'], ss_max=payload['ss_max']),
        loose=PotentialCriteria(ss_min=payload['pot_min'], ss_max=payload['pot_max']),
        majority=payload['majority'],
        model_id_regex=payload['model_id_regex'],
        model_dir_regex=payload['model_dir_regex'],
        skip_non_cys=payload['skip_non_cys'],
    )
    return scanner.scan_variant_folder(payload['folder'])

def cli_denovo():
    ap = argparse.ArgumentParser(
        description=('Scan per-variant model folders for de novo disulfides (majority vote). '
                     'Handles IA1124 (implicit Cys) and older styles; expects model_N/*.pdb or filenames like *_X_Y.pdb.')
    )
    ap.add_argument('root', help='Root with per-variant folders (e.g., IA1124, KA792, VA260_C)')
    ap.add_argument('--out', default='denovo_disulfide_summary.csv')
    ap.add_argument('--details', default='denovo_disulfide_details.csv')
    ap.add_argument('--majority', type=float, default=0.5, help='Fraction of positive models to call de_novo_disulfide')
    ap.add_argument('--ss-min', type=float, default=1.95, help='Strict SG–SG min (Å) for formed disulfides')
    ap.add_argument('--ss-max', type=float, default=2.25, help='Strict SG–SG max (Å) for formed disulfides')
    ap.add_argument('--pot-min', type=float, default=2.25, help='Loose SG–SG min (Å) for potential disulfides')
    ap.add_argument('--pot-max', type=float, default=3.5, help='Loose SG–SG max (Å) for potential disulfides')
    ap.add_argument('--model-id-regex', type=str, default=r'_(?P<mid>\d+)\.pdb$', help='Regex to extract model id from filename; must define group "mid"')
    ap.add_argument('--model-dir-regex', type=str, default=r'model_(?P<mid>\d+)$', help='Regex to extract model id from parent dir; must define group "mid"')
    ap.add_argument('--jobs', type=int, default=1, help='Parallel workers for variants (use 1 for serial)')
    ap.add_argument('--verbose', action='store_true', help='Print per-variant progress messages')
    ap.add_argument('--quiet', action='store_true', help='Minimal output')
    ap.add_argument('--no-skip-non-cys', dest='skip_non_cys', action='store_false', help='Also scan variants that do NOT introduce a cysteine')
    args = ap.parse_args()

    folders = [os.path.join(args.root, d) for d in os.listdir(args.root) if os.path.isdir(os.path.join(args.root, d))]
    folders.sort()

    if not args.quiet:
        print(f"🔍 Scanning {len(folders)} variant folders under {args.root} ...")
    t0 = time.time()

    payload_common = {
        'ss_min': args.ss_min,
        'ss_max': args.ss_max,
        'pot_min': args.pot_min,
        'pot_max': args.pot_max,
        'majority': args.majority,
        'model_id_regex': args.model_id_regex,
        'model_dir_regex': args.model_dir_regex,
        'skip_non_cys': args.skip_non_cys,
    }

    summaries, details = [], []

    if args.jobs and args.jobs > 1:
        if args.verbose and not args.quiet:
            print(f"⚙️  Parallel mode with {args.jobs} workers")
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            futures = {}
            for i, fd in enumerate(folders, start=1):
                payload = dict(payload_common)
                payload['folder'] = fd
                futures[ex.submit(_worker_scan_folder, payload)] = (i, fd)
            for fut in as_completed(futures):
                i, fd = futures[fut]
                res = fut.result()
                summaries.append({k: v for k, v in res.items() if k != 'details'})
                for dr in res['details']:
                    details.append({'variant_folder': res.get('variant_folder', ''), **dr})
                if args.verbose and not args.quiet:
                    print(f"▶ ({i}/{len(folders)}) {os.path.basename(fd)} → models: {res['n_models']}, positives: {res['n_positive']} "
                          f"({res['fraction_positive']*100:.1f}%), class: {res['classification']}")
    else:
        for i, fd in enumerate(folders, start=1):
            payload = dict(payload_common)
            payload['folder'] = fd
            res = _worker_scan_folder(payload)
            summaries.append({k: v for k, v in res.items() if k != 'details'})
            for dr in res['details']:
                details.append({'variant_folder': res.get('variant_folder', ''), **dr})
            if args.verbose and not args.quiet:
                print(f"▶ ({i}/{len(folders)}) {os.path.basename(fd)} → models: {res['n_models']}, positives: {res['n_positive']} "
                      f"({res['fraction_positive']*100:.1f}%), class: {res['classification']}")

    pd.DataFrame(summaries).to_csv(args.out, index=False)
    pd.DataFrame(details).to_csv(args.details, index=False)

    if not args.quiet:
        dt = time.time() - t0
        print(f"\n✅ Completed in {_human_time(dt)} — results saved to:\n  {args.out}\n  {args.details}")

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1 and sys.argv[1] in {'wt', 'denovo'}:
        sub = sys.argv.pop(1)
        if sub == 'wt':
            cli_wt()
        else:
            cli_denovo()
    else:
        print('Usage: python disulfide_module.py wt <WT.pdb> [--rotamer-scan --rot-angles "-60,60,180"] [--native-mode STRICT|LOOSE]')
        print('   or: python disulfide_module.py denovo <ROOT> [--jobs N] [--verbose]')

