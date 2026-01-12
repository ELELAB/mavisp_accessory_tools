#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MAVISp Disulfide Module — WT (static/FoldX) + De Novo (parallel)
================================================================

WT mode:
- Static native detection on a single WT structure (STRICT/LOOSE).
- Optional FoldX WT self-mutation scan: recursively read PDBs under --foldx-wt and
  evaluate Cys–Cys geometry across those models. Choose which source(s) define
  native disulfides via --native-source (static|foldx|both).

De novo mode:
- Variant folders like IA1124 / YA739 (implicit Cys at chain A, position) supported.
- Skip folders starting with 'C' (WT cysteine positions).
- Only mutant PDBs named like *_X_Y.pdb are scanned (WT PDBs ignored).
- Majority vote over models where the specified site is actually Cys.
- Parallel with --jobs, verbose progress with --verbose.
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

# ---------------------------- PDB utilities ------------------------------ #

def read_cysteine_sg_atoms(pdb_path: str):
    """Return (list_of_SG_atoms, residue_map) for cysteines (CYS/CYX/CYM/CYD/CSE...)."""
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
                for an in ("SG","CB","CA"):
                    if an in residue:
                        a = residue[an]
                        atomdict[an] = AtomRef(
                            ch, resname, int(resseq),
                            icode.strip() if isinstance(icode, str) else "",
                            an, a.coord.copy()
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
    use_angles: bool = False  # often distance-only in loose mode

# ------------------------------ WT Analyzer ------------------------------ #

class WTDisulfideAnalyzer:
    def __init__(self,
                 strict: DisulfideCriteria,
                 loose: LooseCriteria,
                 rotamer_scan: bool = False,
                 rot_angles: Optional[List[float]] = None,
                 rot_threshold: float = 10.0):
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
            # rotate SG around CB-CA
            SGa_rot = _rotate_point_around_axis(SGa, CBa, axisA, angA)
            for angB in self.rot_angles:
                SGb_rot = _rotate_point_around_axis(SGb, CBb, axisB, angB)
                d = dist(SGa_rot, SGb_rot)
                if d < best_d:
                    best_d = d
                    best_ang = (angA, angB)
        return best_d, best_ang

    def analyze_static(self, pdb_path: str, native_mode: str = 'STRICT',
                       dump_pairs: Optional[str] = None, debug: bool = False) -> pd.DataFrame:
        """Analyze native disulfides on a single WT structure (static)."""
        native_mode = native_mode.upper()
        rows = []
        _, resmap = read_cysteine_sg_atoms(pdb_path)
        keys = list(resmap.keys())
        if debug:
            print(f"[Static] Found {len(keys)} cysteines with SG; chains: {sorted(set(k[0] for k in keys))}")
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
                    ok, met = self._loose_pair(a_map, b_map, use_angles=False)

                if ok:
                    row = {
                        "chainA": akey[0], "resA": akey[1], "icodeA": akey[2],
                        "chainB": bkey[0], "resB": bkey[1], "icodeB": bkey[2],
                        "ok": True, **met, "source": "static",
                        "min_rot_d_ss": float(min_rot_d), "rotamer_reachable": bool(rot_reachable)
                    }
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

    def _loose_pair(self, a_map: Dict[str, AtomRef], b_map: Dict[str, AtomRef], use_angles: bool = False):
        if "SG" not in a_map or "SG" not in b_map:
            return False, {}
        d = dist(a_map['SG'].coord, b_map['SG'].coord)
        ok = self.loose.ss_min <= d <= self.loose.ss_max
        met = {"d_ss": d}
        if ok and use_angles and ("CB" in a_map and "CB" in b_map):
            a1 = angle(a_map['CB'].coord, a_map['SG'].coord, b_map['SG'].coord)
            a2 = angle(a_map['SG'].coord, b_map['SG'].coord, b_map['CB'].coord)
            met.update({"angle_cb_sg_sg": a1, "angle_sg_sg_cb": a2})
            ok = (self.strict.cb_sg_sg_angle_min <= a1 <= self.strict.cb_sg_sg_angle_max) and \
                 (self.strict.sg_sg_cb_angle_min <= a2 <= self.strict.sg_sg_cb_angle_max)
        return ok, met

# --- Rotation helper used by WTDisulfideAnalyzer._min_rotamer_distance --- #

def _rotate_point_around_axis(point: np.ndarray, axis_point: np.ndarray, axis_dir: np.ndarray, angle_deg: float) -> np.ndarray:
    p = point - axis_point
    u = axis_dir / (np.linalg.norm(axis_dir) + 1e-12)
    theta = math.radians(angle_deg)
    p_rot = p * math.cos(theta) + np.cross(u, p) * math.sin(theta) + u * np.dot(u, p) * (1 - math.cos(theta))
    return axis_point + p_rot

# ----------------------- FoldX WT self-mutation scan --------------------- #

def _list_pdbs_recursive(root: str) -> List[str]:
    """Recursively collect .pdb/.PDB under root (follow symlinks)."""
    root_abs = os.path.abspath(os.path.expanduser(root))
    pdbs = []
    for r, dnames, fnames in os.walk(root_abs, followlinks=True):
        for fn in fnames:
            if fn.lower().endswith('.pdb'):
                pdbs.append(os.path.join(r, fn))
    pdbs.sort()
    return pdbs

def analyze_foldx_wt(foldx_root: str, mode: str,
                     strict: DisulfideCriteria, loose: LooseCriteria,
                     debug: bool = False) -> pd.DataFrame:
    """
    Scan all .pdb files recursively under foldx_root as alternative WT models
    built by FoldX (self-mutations). Return pairs that satisfy criteria in
    at least one model. mode: 'STRICT' or 'LOOSE'.
    """
    mode = mode.upper()
    rows = []
    root_abs = os.path.abspath(os.path.expanduser(foldx_root))

    if debug:
        print(f"[FoldX] CWD: {os.getcwd()}")
        print(f"[FoldX] Root (as passed): {foldx_root}")
        print(f"[FoldX] Root (resolved):  {root_abs}")
        if os.path.isdir(root_abs):
            preview = [d for d in os.listdir(root_abs) if os.path.isdir(os.path.join(root_abs, d))]
            print(f"[FoldX] First-level subfolders (up to 10): {preview[:10]}")
        else:
            print(f"[FoldX] WARNING: resolved path does not exist or is not a directory")

    pdbs = _list_pdbs_recursive(root_abs)
    if debug:
        print(f"[FoldX] Found {len(pdbs)} PDB(s) under {root_abs}")
        for p in pdbs[:8]:
            print(f"[FoldX] example: {p}")

    for mpath in pdbs:
        try:
            _, resmap = read_cysteine_sg_atoms(mpath)
        except Exception:
            continue
        keys = list(resmap.keys())
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                akey, bkey = keys[i], keys[j]
                a_map, b_map = resmap[akey], resmap[bkey]
                if mode == 'STRICT':
                    ok, met = WTDisulfideAnalyzer(strict, loose)._strict_pair(a_map, b_map)
                else:
                    ok, met = WTDisulfideAnalyzer(strict, loose)._loose_pair(a_map, b_map, use_angles=False)
                if ok:
                    rows.append({
                        "chainA": akey[0], "resA": akey[1], "icodeA": akey[2],
                        "chainB": bkey[0], "resB": bkey[1], "icodeB": bkey[2],
                        "ok": True, **met, "source": "foldx", "model": os.path.basename(mpath)
                    })
    if debug:
        print(f"[FoldX] Valid pairs found: {len(rows)}")
    return pd.DataFrame(rows)

# -------------------------- Mutlist export helper ------------------------ #

def export_mutlist_simple(df_native: pd.DataFrame, variants_file: str, out_tsv: str, default_chain: str = 'A'):
    """
    From a table of native Cys–Cys pairs (any source), export Cys->nonCys disruptions.
    Input variants_file is a simple list like: A260V (no header).
    Assumes chain given by default_chain (e.g., 'A').
    """
    posset = set()
    if df_native is not None and not df_native.empty:
        for _, r in df_native.iterrows():
            chainA = str(r.get('chainA', 'A'))
            chainB = str(r.get('chainB', 'A'))
            resA = int(r.get('resA', -1))
            resB = int(r.get('resB', -1))
            if resA > 0: posset.add((chainA, resA))
            if resB > 0: posset.add((chainB, resB))

    if not posset:
        pd.DataFrame(columns=["variant_id","reason"]).to_csv(out_tsv, sep='\t', index=False)
        return

    vdf = pd.read_csv(variants_file, header=None, names=['tok'])
    rows = []
    for _, v in vdf.iterrows():
        s = str(v['tok']).strip()
        m = re.match(r'^([A-Z])(\d+)([A-Z])$', s)   # e.g. A260V
        if not m:
            continue
        wt, pos, mut = m.groups()
        pos = int(pos)
        ch = default_chain
        if wt == 'C' and (ch, pos) in posset and mut != 'C':
            rows.append({
                'variant_id': f'{wt}{ch}{pos}{mut}',
                'chain': ch, 'pos': pos, 'wt': wt, 'mut': mut,
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
    def __init__(self,
                 strict: DisulfideCriteria,
                 loose: PotentialCriteria,
                 majority: float = 0.5,
                 model_id_regex: str = r"_(?P<mid>\d+)\.pdb$",
                 model_dir_regex: str = r"model_(?P<mid>\d+)$",
                 skip_non_cys: bool = True):
        self.strict = strict
        self.loose = loose
        self.majority = majority
        self._re_file_id = re.compile(model_id_regex, re.IGNORECASE)
        self._re_dir_id  = re.compile(model_dir_regex, re.IGNORECASE)
        self.skip_non_cys = skip_non_cys

    @staticmethod
    def parse_folder_mutation(folder_name: str) -> Optional[Dict[str, str]]:
        """
        Accepts:
          - IA1124 / YA739 / WA381 ...  (WT=<letter>, chain=<letter>, pos=<int>) → implicit mut='C'
          - Also supports VA260_C / AA273_A / A_A_273_C styles.
        """
        m = re.match(r"^(?P<wt>[A-Z])(?P<chain>[A-Za-z])(?P<pos>\d+)$", folder_name)
        if m:
            d = m.groupdict()
            return {"wt": d["wt"], "chain": d["chain"], "pos": d["pos"], "mut": "C"}
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

        n_total = len(models)   # all mutant PDBs found by naming
        n = 0                   # models where the specified site is actually Cys
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
    if h: return f"{h}h{m}m{s}s"
    if m: return f"{m}m{s}s"
    return f"{s}s"

# ------------------------------ CLIs ------------------------------------- #

def cli_wt():
    ap = argparse.ArgumentParser(
        description='WT disulfide analysis: static (STRICT/LOOSE) + optional FoldX WT models; export disruption mutlist.'
    )
    ap.add_argument('pdb', help='WT PDB/AF model file')
    ap.add_argument('--out-csv', default='native_disulfides.csv', help='Output CSV of native SS pairs')
    # Static options
    ap.add_argument('--native-mode', choices=['STRICT','LOOSE'], default='LOOSE', help='Static WT evaluation criteria')
    ap.add_argument('--ss-min', type=float, default=1.95, help='Strict SG–SG min (Å)')
    ap.add_argument('--ss-max', type=float, default=2.25, help='Strict SG–SG max (Å)')
    ap.add_argument('--ss-loose-min', type=float, default=1.8, help='LOOSE SG–SG min (Å)')
    ap.add_argument('--ss-loose-max', type=float, default=3.5, help='LOOSE SG–SG max (Å)')
    ap.add_argument('--no-angle', action='store_true', help='In LOOSE mode, do not require angles')
    # Rotamer scan (optional; independent of FoldX)
    ap.add_argument('--rotamer-scan', action='store_true', help='Enable coarse rotamer scan for close pairs')
    ap.add_argument('--rot-angles', type=str, default='-60,60,180', metavar='ANGLES', help='Comma-separated χ1-like angles (deg), quote the list')
    ap.add_argument('--rot-threshold', type=float, default=10.0, help='Only scan pairs with initial SG–SG ≤ this (Å)')
    ap.add_argument('--dump-pairs', default=None, help='Write all Cys–Cys pairs to this CSV')
    # FoldX WT support
    ap.add_argument('--foldx-wt', default=None, help='Root folder with FoldX WT self-mutation models (recursively scanned)')
    ap.add_argument('--foldx-mode', choices=['STRICT','LOOSE'], default='LOOSE', help='Criteria for evaluating FoldX WT models')
    ap.add_argument('--native-source', choices=['static','foldx','both'], default='both', help='Which sources define native SS sites')
    # Mutlist export
    ap.add_argument('--variants-file', default=None, help='Plain list: A260V, C360V ... (no header)')
    ap.add_argument('--chain', default='A', help='Assumed chain for simple list tokens')
    ap.add_argument('--out-mutlist', default='mutlist_disulfide_disruptions.tsv')
    ap.add_argument('--debug', action='store_true')
    args = ap.parse_args()

    strict = DisulfideCriteria(ss_min=args.ss_min, ss_max=args.ss_max)
    loose = LooseCriteria(ss_min=args.ss_loose_min, ss_max=args.ss_loose_max, use_angles=(not args.no_angle))

    # Static WT evaluation
    static_df = pd.DataFrame()
    if args.native_source in ('static', 'both'):
        try:
            rot_angles = [float(x) for x in args.rot_angles.split(',')] if args.rotamer_scan else None
        except Exception:
            rot_angles = [-60.0, 60.0, 180.0]
        analyzer = WTDisulfideAnalyzer(strict, loose, rotamer_scan=args.rotamer_scan,
                                       rot_angles=rot_angles, rot_threshold=args.rot_threshold)
        static_df = analyzer.analyze_static(args.pdb, native_mode=args.native_mode,
                                            dump_pairs=args.dump_pairs, debug=args.debug)

    # FoldX WT evaluation
    foldx_df = pd.DataFrame()
    if args.foldx_wt and args.native_source in ('foldx', 'both'):
        foldx_df = analyze_foldx_wt(args.foldx_wt, args.foldx_mode, strict, loose, debug=args.debug)

    # Combine according to native-source
    if args.native_source == 'static':
        out_df = static_df
    elif args.native_source == 'foldx':
        out_df = foldx_df
    else:
        out_df = pd.concat([static_df, foldx_df], ignore_index=True) if (not static_df.empty or not foldx_df.empty) else pd.DataFrame(columns=[
            "chainA","resA","icodeA","chainB","resB","icodeB","ok","d_ss","source"
        ])

    # De-duplicate by pair + source
    if not out_df.empty:
        out_df = out_df.drop_duplicates(subset=["chainA","resA","chainB","resB","source"])

    out_df.to_csv(args.out_csv, index=False)
    if args.debug:
        print(out_df if not out_df.empty else "No native disulfides detected with current settings/sources.")

    if args.variants_file:
        export_mutlist_simple(out_df, args.variants_file, args.out_mutlist, default_chain=args.chain)

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
                     'Handles IA1124/YA739 implicit Cys; skips folders starting with C; expects *_X_Y.pdb mutant files.')
    )
    ap.add_argument('root', help='Root with per-variant folders (e.g., IA1124, YA739, VA260_C)')
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

    # Discover variant folders and skip those starting with 'C' (WT Cys positions)
    raw_folders = [
        os.path.join(args.root, d)
        for d in os.listdir(args.root)
        if os.path.isdir(os.path.join(args.root, d))
    ]
    folders = [f for f in raw_folders if not os.path.basename(f).startswith('C')]
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
    if len(sys.argv) > 1 and sys.argv[1] in {'wt','denovo'}:
        sub = sys.argv.pop(1)
        if sub == 'wt':
            cli_wt()
        else:
            cli_denovo()
    else:
        print('Usage:')
        print('  WT:     python disulfide_module.py wt <WT.pdb> [--foldx-wt DIR] [--native-source static|foldx|both] [--debug]')
        print('  Denovo: python disulfide_module.py denovo <ROOT> [--jobs N] [--verbose]')

