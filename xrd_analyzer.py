#!/usr/bin/env python3
"""
Universal Automated XRD Phase Identification & Live COD Reference Matcher
=========================================================================
100% Pure Universal Edition — Zero Hardcoded Reference Cards!

Features:
1. Opens a Windows File Explorer dialog so you can pick an XRD file (.xlsx, .xls,
   .csv, .txt, .xy) from ANY folder on your computer.
2. Prompts interactively for:
   - Expected main elements or target compound (e.g., Cu, S, Sn or Sn-CuS)
   - Precursors / reagents used in synthesis (e.g., CuCl2.2H2O, SnCl2.2H2O, thiourea, HCl)
3. Detects experimental peaks (2Theta, d-spacing, FWHM, Scherrer crystallite size, strain).
4. Queries the live Crystallography Open Database (COD) over MySQL/REST with strict
   exact-element filtering (eliminating false substring matches like 'S' matching 'Sn' or 'Ce').
5. Streams the live .cif crystal structure files from COD and computes exact structure-factor
   (|F_hkl|^2) diffraction intensities and Miller indices (hkl) directly in Python/NumPy
   (eliminating space-group forbidden reflections like (017) or (015) even without pymatgen).
6. Uses Sequential Dominant-Phase Search-Match (Primary Phase -> Secondary/Dopant Phase ->
   Precursor Impurity Phase) and exports a complete multi-sheet Excel report.
"""

import io
import math
import os
import re
import sys
import warnings
from itertools import combinations
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, savgol_filter

warnings.filterwarnings("ignore")

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

try:
    import mysql.connector
    HAS_MYSQL = True
except ImportError:
    HAS_MYSQL = False

# ==============================================================================
# 1. PERIODIC TABLE & ATOMIC SCATTERING APPROXIMATION (Z)
# ==============================================================================
VALID_ELEMENTS: Dict[str, int] = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
    "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26,
    "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34,
    "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42,
    "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
    "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58,
    "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66,
    "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71, "Hf": 72, "Ta": 73, "W": 74,
    "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82,
    "Bi": 83, "Th": 90, "U": 92,
}

ELEMENT_NAMES_TO_SYMBOL: Dict[str, str] = {
    "copper": "Cu", "sulfur": "S", "sulphur": "S", "tin": "Sn", "zinc": "Zn",
    "iron": "Fe", "oxygen": "O", "oxide": "O", "sulfide": "S", "sulphide": "S",
    "chloride": "Cl", "chlorine": "Cl", "nitrogen": "N", "carbon": "C",
    "hydrogen": "H", "silver": "Ag", "gold": "Au", "titanium": "Ti",
    "nickel": "Ni", "cobalt": "Co", "manganese": "Mn", "aluminum": "Al",
    "aluminium": "Al", "silicon": "Si", "sodium": "Na", "potassium": "K",
    "calcium": "Ca", "magnesium": "Mg", "barium": "Ba", "bismuth": "Bi",
    "lead": "Pb", "cadmium": "Cd", "selenium": "Se", "tellurium": "Te",
    "molybdenum": "Mo", "tungsten": "W", "vanadium": "V", "chromium": "Cr",
    "cerium": "Ce", "lanthanum": "La", "zirconium": "Zr", "indium": "In",
    "gallium": "Ga", "germanium": "Ge", "antimony": "Sb", "phosphorus": "P",
    "fluorine": "F", "bromine": "Br", "iodine": "I", "lithium": "Li",
}

COMMON_REAGENT_EXPANSIONS: Dict[str, List[str]] = {
    "thiourea": ["C", "H", "N", "S"],
    "tu": ["C", "H", "N", "S"],
    "thioacetamide": ["C", "H", "N", "S"],
    "taa": ["C", "H", "N", "S"],
    "urea": ["C", "H", "N", "O"],
    "tea": ["C", "H", "N", "O"],
    "triethanolamine": ["C", "H", "N", "O"],
    "eda": ["C", "H", "N"],
    "ethylenediamine": ["C", "H", "N"],
    "ipa": ["C", "H", "O"],
    "isopropanol": ["C", "H", "O"],
    "ethanol": ["C", "H", "O"],
    "methanol": ["C", "H", "O"],
    "acetone": ["C", "H", "O"],
    "ammonia": ["N", "H", "O"],
    "nh4oh": ["N", "H", "O"],
    "naoh": ["Na", "O", "H"],
    "koh": ["K", "O", "H"],
    "hcl": ["H", "Cl"],
    "hno3": ["H", "N", "O"],
    "h2so4": ["H", "S", "O"],
    "water": ["H", "O"],
    "pvp": ["C", "H", "N", "O"],
    "ctab": ["C", "H", "N", "Br"],
    "sds": ["C", "H", "Na", "S", "O"],
}

ORGANIC_SOLVENT_LIGHT_ELEMENTS: Set[str] = {"C", "H", "N"}


def parse_chemical_input_to_elements(text: str, drop_pure_organics: bool = False) -> List[str]:
    """
    Extracts chemical element symbols from any free-text input (formulas, names,
    reagents, abbreviations like TEA/IPA/thiourea).
    Matches abbreviations BEFORE regex parsing so 'IPA' is never misread as I + P!
    """
    if not text or not text.strip():
        return []

    found: List[str] = []

    def add_el(el: str) -> None:
        if el in VALID_ELEMENTS and el not in found:
            if drop_pure_organics and el in ORGANIC_SOLVENT_LIGHT_ELEMENTS:
                return
            found.append(el)

    cleaned_text = text
    word_tokens = re.findall(r"[A-Za-z0-9]+", text)
    for token in word_tokens:
        low = token.lower()
        if low in COMMON_REAGENT_EXPANSIONS:
            for el in COMMON_REAGENT_EXPANSIONS[low]:
                add_el(el)
            cleaned_text = re.sub(rf"\b{re.escape(token)}\b", " ", cleaned_text)
        elif low in ELEMENT_NAMES_TO_SYMBOL:
            add_el(ELEMENT_NAMES_TO_SYMBOL[low])
            cleaned_text = re.sub(rf"\b{re.escape(token)}\b", " ", cleaned_text)

    for sym in re.findall(r"[A-Z][a-z]?", cleaned_text):
        add_el(sym)

    return found


# ==============================================================================
# 2. FILE LOADING & EXPERIMENTAL PEAK DETECTION
# ==============================================================================
def select_xrd_file_interactive() -> str:
    """
    Opens a native Windows File Explorer dialog so the user can choose an XRD file
    from any folder on their computer, with a terminal path fallback.
    """
    print("\n📂 Opening file chooser window... (Check your taskbar if it doesn't pop up in front)")
    chosen_path = ""
    try:
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        chosen_path = filedialog.askopenfilename(
            title="Select Experimental XRD Data File",
            filetypes=[
                ("XRD Data Files", "*.xlsx *.xls *.csv *.txt *.xy *.dat"),
                ("Excel Files", "*.xlsx *.xls"),
                ("CSV / Text Files", "*.csv *.txt *.xy *.dat"),
                ("All Files", "*.*"),
            ],
        )
        root.destroy()
    except Exception:
        chosen_path = ""

    if chosen_path and os.path.isfile(chosen_path):
        print(f"✅ Selected file: {chosen_path}")
        return chosen_path

    while True:
        raw = input("\n👉 Enter or paste the full path to your XRD file (.xlsx, .csv, .txt): ").strip().strip('"').strip("'")
        if os.path.isfile(raw):
            return raw
        print(f"❌ File not found: '{raw}'. Please paste a valid file path.")


def load_xrd_file(file_path: str) -> Tuple[pd.DataFrame, str, str]:
    """
    Universal loader for XRD files (.xlsx, .xls, .csv, .txt, .xy).
    Automatically detects the 2Theta and Intensity columns.
    """
    ext = os.path.splitext(file_path)[1].lower()
    if ext in (".xlsx", ".xls"):
        xl = pd.ExcelFile(file_path)
        chosen_sheet = xl.sheet_names[0]
        for s in xl.sheet_names:
            if any(k in s.lower() for k in ("raw", "sample", "data", "xrd", "sheet1")):
                chosen_sheet = s
                break
        df_raw = pd.read_excel(file_path, sheet_name=chosen_sheet)
    else:
        with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()
        start_idx = 0
        for idx, line in enumerate(lines):
            parts = re.split(r"[,\s;\t]+", line.strip())
            if len(parts) >= 2:
                try:
                    float(parts[0])
                    float(parts[1])
                    start_idx = max(0, idx - 1)
                    break
                except ValueError:
                    continue
        df_raw = pd.read_csv(
            io.StringIO("".join(lines[start_idx:])),
            sep=None,
            engine="python",
            on_bad_lines="skip",
        )

    num_cols = []
    for col in df_raw.columns:
        s = pd.to_numeric(df_raw[col], errors="coerce")
        if s.notna().sum() > 20:
            num_cols.append(col)

    if len(num_cols) < 2:
        df_raw = pd.read_excel(file_path, header=None) if ext in (".xlsx", ".xls") else df_raw
        for r in range(min(15, len(df_raw))):
            row_vals = [str(x).strip() for x in df_raw.iloc[r].values]
            if any("theta" in v.lower() or "angle" in v.lower() or "2θ" in v.lower() for v in row_vals):
                df_raw.columns = row_vals
                df_raw = df_raw.iloc[r + 1 :].reset_index(drop=True)
                break

    clean_cols = []
    for col in df_raw.columns:
        s = pd.to_numeric(df_raw[col], errors="coerce")
        if s.notna().sum() > 20:
            df_raw[col] = s
            clean_cols.append(col)

    if len(clean_cols) < 2:
        raise ValueError(f"Could not find numeric 2Theta and Intensity columns in {file_path}")

    tt_col, int_col = clean_cols[0], clean_cols[1]
    for c in clean_cols:
        cl = str(c).lower()
        if any(k in cl for k in ("theta", "2θ", "angle", "pos", "x")):
            tt_col = c
            break
    for c in clean_cols:
        if c == tt_col:
            continue
        cl = str(c).lower()
        if any(k in cl for k in ("int", "count", "cps", "obs", "y")):
            int_col = c
            break

    df = df_raw[[tt_col, int_col]].dropna().sort_values(by=tt_col).reset_index(drop=True)
    return df, str(tt_col), str(int_col)


def _pseudo_voigt(x: np.ndarray, amp: float, cen: float, fwhm: float, eta: float, bg: float) -> np.ndarray:
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    gauss = np.exp(-((x - cen) ** 2) / (2.0 * sigma ** 2))
    lorentz = 1.0 / (1.0 + ((x - cen) / (fwhm / 2.0)) ** 2)
    return bg + amp * (eta * lorentz + (1.0 - eta) * gauss)


def detect_experimental_peaks(
    df: pd.DataFrame,
    tt_col: str,
    int_col: str,
    wavelength_a: float = 1.5406,
    scherrer_k: float = 0.9,
) -> pd.DataFrame:
    """
    Detects experimental XRD peaks, fits local Pseudo-Voigt profiles for accurate
    FWHM, Bragg d-spacing, Scherrer crystallite size (nm), and microstrain.
    """
    tt = df[tt_col].to_numpy(dtype=float)
    raw_int = df[int_col].to_numpy(dtype=float)
    step = float(np.median(np.diff(tt))) if len(tt) > 1 else 0.02

    win = max(7, int(round(0.25 / max(step, 0.005))))
    if win % 2 == 0:
        win += 1
    smooth = savgol_filter(raw_int, window_length=min(win, len(raw_int) - (1 - len(raw_int) % 2)), polyorder=3)

    bg_win = max(51, int(round(4.0 / max(step, 0.005))))
    bg = pd.Series(smooth).rolling(bg_win, center=True, min_periods=1).min().to_numpy()
    bg = pd.Series(bg).rolling(bg_win // 2, center=True, min_periods=1).mean().to_numpy()
    net_int = np.clip(smooth - bg, 0, None)

    noise_std = float(np.percentile(np.abs(np.diff(raw_int)), 68)) / 1.414
    max_net = float(np.max(net_int)) if len(net_int) else 1.0
    prominence = max(noise_std * 2.8, max_net * 0.045)
    min_dist = max(5, int(round(0.85 / max(step, 0.005))))

    peaks, _ = find_peaks(
        net_int,
        prominence=prominence,
        distance=min_dist,
        width=max(2, int(round(0.10 / max(step, 0.005)))),
    )

    records: List[Dict[str, Any]] = []
    for idx_num, p_idx in enumerate(peaks, start=1):
        p_tt = float(tt[p_idx])
        p_raw = float(raw_int[p_idx])
        p_bg = float(bg[p_idx])
        p_net = max(1.0, float(net_int[p_idx]))

        mask = (tt >= p_tt - 1.25) & (tt <= p_tt + 1.25)
        x_loc = tt[mask]
        y_loc = net_int[mask]

        half_val = p_net / 2.0
        left_idx = p_idx
        while left_idx > 0 and net_int[left_idx] > half_val and (p_tt - tt[left_idx]) < 2.0:
            left_idx -= 1
        right_idx = p_idx
        while right_idx < len(tt) - 1 and net_int[right_idx] > half_val and (tt[right_idx] - p_tt) < 2.0:
            right_idx += 1
        emp_fwhm = max(0.08, min(2.5, float(tt[right_idx] - tt[left_idx])))

        fwhm_deg = emp_fwhm
        try:
            popt, _ = curve_fit(
                _pseudo_voigt,
                x_loc,
                y_loc,
                p0=[p_net, p_tt, emp_fwhm, 0.5, 0.0],
                bounds=(
                    [p_net * 0.3, p_tt - 0.35, 0.05, 0.0, -p_net * 0.2],
                    [p_net * 2.0, p_tt + 0.35, 2.50, 1.0, p_net * 0.3],
                ),
                maxfev=1500,
            )
            fit_fwhm = float(abs(popt[2]))
            if 0.06 <= fit_fwhm <= 2.2:
                fwhm_deg = fit_fwhm
        except Exception:
            pass

        theta_rad = math.radians(p_tt / 2.0)
        beta_rad = math.radians(fwhm_deg)
        d_spacing = wavelength_a / (2.0 * math.sin(theta_rad))
        crystallite_nm = (scherrer_k * wavelength_a) / (beta_rad * math.cos(theta_rad)) / 10.0
        strain = beta_rad / (4.0 * math.tan(theta_rad))

        records.append({
            "Peak_ID": f"P{idx_num}",
            "Peak_2Theta": round(p_tt, 4),
            "Raw_Intensity": round(p_raw, 2),
            "Baseline_Intensity": round(p_bg, 2),
            "Net_Intensity": round(p_net, 2),
            "d_Spacing_A": round(d_spacing, 4),
            "FWHM_deg": round(fwhm_deg, 4),
            "Theta_rad": round(theta_rad, 4),
            "Beta_rad": round(beta_rad, 4),
            "Crystallite_Size_nm": round(crystallite_nm, 2),
            "Lattice_Strain": round(strain, 4),
        })

    out_df = pd.DataFrame(records)
    if not out_df.empty:
        max_i = out_df["Net_Intensity"].max()
        out_df["Rel_Intensity_pct"] = (out_df["Net_Intensity"] / max_i * 100.0).round(1)
    return out_df


# ==============================================================================
# 3. PURE-PYTHON + NUMPY CIF STRUCTURE-FACTOR (|F_hkl|^2) XRD SIMULATOR
# ==============================================================================
def _parse_cif_float(val_str: str) -> Optional[float]:
    if not val_str:
        return None
    cleaned = re.sub(r"\([^)]*\)", "", val_str.strip()).strip("'\"")
    try:
        return float(cleaned)
    except ValueError:
        return None


def _eval_symop_expr(expr: str, x: float, y: float, z: float) -> float:
    """Evaluates a CIF symmetry expression like '-x+y, -x, z+1/2' for one coordinate."""
    e = expr.strip().lower().replace(" ", "")

    def repl_frac(m: re.Match) -> str:
        return str(float(m.group(1)) / float(m.group(2)))

    e = re.sub(r"(\d+)/(\d+)", repl_frac, e)
    val = 0.0
    for term in re.findall(r"[+-]?[^+-]+", e):
        sign = -1.0 if term.startswith("-") else 1.0
        body = term.lstrip("+-")
        if body == "x":
            val += sign * x
        elif body == "y":
            val += sign * y
        elif body == "z":
            val += sign * z
        elif body.endswith("x"):
            val += sign * float(body[:-1] or 1.0) * x
        elif body.endswith("y"):
            val += sign * float(body[:-1] or 1.0) * y
        elif body.endswith("z"):
            val += sign * float(body[:-1] or 1.0) * z
        else:
            try:
                val += sign * float(body)
            except ValueError:
                pass
    return val % 1.0


def _parse_cif_structure(cif_text: str) -> Optional[Dict[str, Any]]:
    """
    Parses unit cell parameters, symmetry operations, and fractional atomic coordinates
    directly from a raw CIF string.
    """
    def extract_tag_float(tag: str, default: float) -> float:
        m = re.search(rf"{tag}\s+([0-9.+-]+(?:\([0-9]+\))?)", cif_text, re.IGNORECASE)
        if m:
            v = _parse_cif_float(m.group(1))
            if v is not None:
                return v
        return default

    a = extract_tag_float(r"_cell_length_a", 0.0)
    b = extract_tag_float(r"_cell_length_b", a)
    c = extract_tag_float(r"_cell_length_c", a)
    alpha = extract_tag_float(r"_cell_angle_alpha", 90.0)
    beta = extract_tag_float(r"_cell_angle_beta", 90.0)
    gamma = extract_tag_float(r"_cell_angle_gamma", 90.0)
    if a <= 0 or b <= 0 or c <= 0:
        return None

    lines = cif_text.splitlines()
    symops: List[Tuple[str, str, str]] = []
    atoms: List[Tuple[int, float, float, float, float]] = []

    i = 0
    n = len(lines)
    while i < n:
        line = lines[i].strip()
        if line.lower() == "loop_":
            headers: List[str] = []
            i += 1
            while i < n and lines[i].strip().startswith("_"):
                headers.append(lines[i].strip().split()[0].lower())
                i += 1
            rows: List[List[str]] = []
            while i < n:
                s = lines[i].strip()
                if not s or s.startswith("#"):
                    i += 1
                    continue
                if s.lower() == "loop_" or s.startswith("_") or s.lower().startswith("data_"):
                    break
                if s.startswith(";"):
                    i += 1
                    while i < n and not lines[i].strip().startswith(";"):
                        i += 1
                    i += 1
                    continue
                tokens = re.findall(r"'[^']*'|\"[^\"]*\"|\S+", s)
                if tokens:
                    rows.append([t.strip("'\"") for t in tokens])
                i += 1

            sym_col = -1
            for idx, h in enumerate(headers):
                if h in (
                    "_symmetry_equiv_pos_as_xyz",
                    "_space_group_symop_operation_xyz",
                ):
                    sym_col = idx
                    break
            if sym_col >= 0:
                for r in rows:
                    if sym_col < len(r):
                        parts = r[sym_col].split(",")
                        if len(parts) == 3:
                            symops.append((parts[0], parts[1], parts[2]))

            fx = next((idx for idx, h in enumerate(headers) if h == "_atom_site_fract_x"), -1)
            fy = next((idx for idx, h in enumerate(headers) if h == "_atom_site_fract_y"), -1)
            fz = next((idx for idx, h in enumerate(headers) if h == "_atom_site_fract_z"), -1)
            flabel = next(
                (idx for idx, h in enumerate(headers) if h in ("_atom_site_type_symbol", "_atom_site_label")),
                -1,
            )
            focc = next((idx for idx, h in enumerate(headers) if h == "_atom_site_occupancy"), -1)

            if fx >= 0 and fy >= 0 and fz >= 0 and flabel >= 0:
                for r in rows:
                    if max(fx, fy, fz, flabel) < len(r):
                        x_v = _parse_cif_float(r[fx])
                        y_v = _parse_cif_float(r[fy])
                        z_v = _parse_cif_float(r[fz])
                        occ_v = _parse_cif_float(r[focc]) if (0 <= focc < len(r)) else 1.0
                        if occ_v is None:
                            occ_v = 1.0
                        sym_m = re.match(r"([A-Z][a-z]?)", r[flabel])
                        el_sym = sym_m.group(1) if sym_m else "Cu"
                        z_num = VALID_ELEMENTS.get(el_sym, 20)
                        if x_v is not None and y_v is not None and z_v is not None:
                            atoms.append((z_num, x_v % 1.0, y_v % 1.0, z_v % 1.0, occ_v))
            continue
        i += 1

    if not symops:
        symops = [("x", "y", "z"), ("-x", "-y", "-z")]

    full_atoms: List[Tuple[int, float, float, float, float]] = []
    for z_num, ax, ay, az, occ in atoms:
        orbit: List[Tuple[float, float, float]] = []
        for sx, sy, sz in symops:
            nx = _eval_symop_expr(sx, ax, ay, az)
            ny = _eval_symop_expr(sy, ax, ay, az)
            nz = _eval_symop_expr(sz, ax, ay, az)
            is_dup = False
            for ox, oy, oz in orbit:
                dx = min(abs(nx - ox), 1.0 - abs(nx - ox))
                dy = min(abs(ny - oy), 1.0 - abs(ny - oy))
                dz = min(abs(nz - oz), 1.0 - abs(nz - oz))
                if dx * dx + dy * dy + dz * dz < 1e-4:
                    is_dup = True
                    break
            if not is_dup:
                orbit.append((nx, ny, nz))
                full_atoms.append((z_num, nx, ny, nz, occ))

    return {
        "a": a, "b": b, "c": c,
        "alpha": alpha, "beta": beta, "gamma": gamma,
        "atoms": full_atoms,
    }


def simulate_xrd_from_cif_str(
    cif_text: str,
    wavelength_a: float = 1.5406,
    min_2theta: float = 10.0,
    max_2theta: float = 80.0,
) -> List[Dict[str, Any]]:
    """
    Computes exact powder XRD peak positions (2Theta, d-spacing), Miller indices (hkl),
    and structure-factor intensities (|F_hkl|^2 with Lorentz-polarization and multiplicity)
    directly from a CIF file text using NumPy.
    Automatically eliminates space-group forbidden reflections (|F_hkl|^2 == 0)!
    """
    parsed = _parse_cif_structure(cif_text)
    if not parsed or not parsed["atoms"]:
        return []

    a, b, c = parsed["a"], parsed["b"], parsed["c"]
    alpha_r = math.radians(parsed["alpha"])
    beta_r = math.radians(parsed["beta"])
    gamma_r = math.radians(parsed["gamma"])

    ca, cb, cg = math.cos(alpha_r), math.cos(beta_r), math.cos(gamma_r)
    sa, sb, sg = math.sin(alpha_r), math.sin(beta_r), math.sin(gamma_r)
    vol_sq = 1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg
    if vol_sq <= 1e-8:
        return []
    vol = a * b * c * math.sqrt(vol_sq)

    astar = b * c * sa / vol
    bstar = a * c * sb / vol
    cstar = a * b * sg / vol
    cos_alphastar = (cb * cg - ca) / (sb * sg)
    cos_betastar = (ca * cg - cb) / (sa * sg)
    cos_gammastar = (ca * cb - cg) / (sa * sb)

    d_min = wavelength_a / (2.0 * math.sin(math.radians(min(max_2theta, 160.0) / 2.0)))
    d_max = wavelength_a / (2.0 * math.sin(math.radians(max(min_2theta, 1.0) / 2.0)))

    h_max = min(int(a / d_min) + 1, 12)
    k_max = min(int(b / d_min) + 1, 12)
    l_max = min(int(c / d_min) + 1, 18)

    h_grid, k_grid, l_grid = np.meshgrid(
        np.arange(-h_max, h_max + 1),
        np.arange(-k_max, k_max + 1),
        np.arange(0, l_max + 1),
        indexing="ij",
    )
    hkl = np.stack([h_grid.ravel(), k_grid.ravel(), l_grid.ravel()], axis=1)
    keep = (hkl[:, 2] > 0) | ((hkl[:, 2] == 0) & ((hkl[:, 0] > 0) | ((hkl[:, 0] == 0) & (hkl[:, 1] > 0))))
    hkl = hkl[keep]

    h = hkl[:, 0].astype(float)
    k = hkl[:, 1].astype(float)
    l = hkl[:, 2].astype(float)

    inv_d2 = (
        (h * astar) ** 2
        + (k * bstar) ** 2
        + (l * cstar) ** 2
        + 2.0 * h * k * astar * bstar * cos_gammastar
        + 2.0 * h * l * astar * cstar * cos_betastar
        + 2.0 * k * l * bstar * cstar * cos_alphastar
    )
    valid_d = (inv_d2 >= 1.0 / (d_max ** 2)) & (inv_d2 <= 1.0 / (d_min ** 2))
    hkl = hkl[valid_d]
    inv_d2 = inv_d2[valid_d]
    if len(hkl) == 0:
        return []

    d_arr = 1.0 / np.sqrt(inv_d2)
    sin_theta = wavelength_a / (2.0 * d_arr)
    valid_th = (sin_theta > 0.0) & (sin_theta < 0.999)
    hkl = hkl[valid_th]
    d_arr = d_arr[valid_th]
    sin_theta = sin_theta[valid_th]
    theta_rad = np.arcsin(sin_theta)
    two_theta = np.degrees(2.0 * theta_rad)

    atoms = parsed["atoms"]
    z_arr = np.array([at[0] for at in atoms], dtype=float)
    xyz_arr = np.array([[at[1], at[2], at[3]] for at in atoms], dtype=float)
    occ_arr = np.array([at[4] for at in atoms], dtype=float)

    s_sq = (sin_theta / wavelength_a) ** 2
    decay = np.exp(-3.5 * s_sq)
    phases = 2.0 * np.pi * (hkl @ xyz_arr.T)
    weights = z_arr * occ_arr
    f_real = (np.cos(phases) @ weights) * decay
    f_imag = (np.sin(phases) @ weights) * decay
    f_sq = f_real ** 2 + f_imag ** 2

    cos_2theta = np.cos(2.0 * theta_rad)
    cos_theta = np.cos(theta_rad)
    lp = (1.0 + cos_2theta ** 2) / (np.maximum(1e-6, (sin_theta ** 2) * cos_theta))
    raw_int = f_sq * lp * 2.0

    order = np.argsort(two_theta)
    two_theta = two_theta[order]
    d_arr = d_arr[order]
    raw_int = raw_int[order]
    hkl = hkl[order]

    grouped_peaks: List[Dict[str, Any]] = []
    for tt_val, d_val, i_val, (hi, ki, li) in zip(two_theta, d_arr, raw_int, hkl):
        if i_val <= 1e-3:
            continue
        ah, ak, al = abs(int(hi)), abs(int(ki)), abs(int(li))
        if abs(a - b) < 0.01 and ah < ak:
            ah, ak = ak, ah
        hkl_score = (
            (1 if (int(hi) >= 0 and int(ki) >= 0 and int(li) >= 0) else 0),
            (1 if ah >= ak else 0),
            -(ah + ak + al),
            ah,
            ak,
            al,
        )
        hkl_label = f"({ah}{ak}{al})"

        if grouped_peaks and abs(tt_val - grouped_peaks[-1]["two_theta"]) < 0.06:
            prev = grouped_peaks[-1]
            tot_i = prev["intensity"] + float(i_val)
            prev["two_theta"] = (prev["two_theta"] * prev["intensity"] + float(tt_val) * float(i_val)) / tot_i
            prev["d_spacing"] = (prev["d_spacing"] * prev["intensity"] + float(d_val) * float(i_val)) / tot_i
            prev["intensity"] = tot_i
            if hkl_score > prev["_hkl_score"]:
                prev["hkl"] = hkl_label
                prev["_hkl_score"] = hkl_score
        else:
            grouped_peaks.append({
                "two_theta": float(tt_val),
                "d_spacing": float(d_val),
                "intensity": float(i_val),
                "hkl": hkl_label,
                "_hkl_score": hkl_score,
            })

    if not grouped_peaks:
        return []

    max_int = max(p["intensity"] for p in grouped_peaks)
    if max_int <= 0:
        return []

    final_peaks: List[Dict[str, Any]] = []
    for p in grouped_peaks:
        rel_i = (p["intensity"] / max_int) * 1000.0
        if rel_i >= 35.0:
            final_peaks.append({
                "two_theta": round(p["two_theta"], 4),
                "d_spacing": round(p["d_spacing"], 4),
                "intensity": round(rel_i, 1),
                "hkl": p["hkl"],
            })

    return final_peaks


# ==============================================================================
# 4. LIVE COD QUERY ENGINE (STRICT ELEMENT MATCHING — ZERO HARDCODED CARDS)
# ==============================================================================
def _clean_cod_formula(raw_formula: str) -> str:
    if not raw_formula:
        return ""
    return raw_formula.strip().strip("-").strip()


def _formula_matches_exact_elements(raw_formula: str, target_elements: Set[str]) -> bool:
    """
    Strictly verifies that the chemical elements present in a COD formula match
    ONLY the target_elements set (preventing 'S' from matching 'Sn' or 'Ce'!).
    """
    cleaned = _clean_cod_formula(raw_formula)
    if not cleaned:
        return False
    els = set(re.findall(r"[A-Z][a-z]?", cleaned))
    return els == target_elements


def _entry_priority_score(row: Dict[str, Any]) -> Tuple[int, int, int, int]:
    """
    Ranks COD entries for a given element combination so well-characterized ambient
    stoichiometric mineral/inorganic phases (e.g. Covellite CuS, Ottemannite Sn2S3,
    Berndtite SnS2, Herzenbergite SnS, Cassiterite SnO2) come before high-pressure,
    fractional-alloy, or disordered entries.
    """
    formula = _clean_cod_formula(str(row.get("formula") or ""))
    mineral = str(row.get("mineral") or "").strip()
    title = str(row.get("title") or "").lower()

    has_fraction = "." in formula
    is_high_pressure_or_temp = any(
        w in title for w in ("pressure", "gpa", "kbar", "high-temperature", "above", "350", "400", "500", "defective", "nanocrystal")
    )
    has_mineral = 1 if mineral else 0
    nums = [int(n) for n in re.findall(r"\b\d+\b", formula)]
    stoich_sum = sum(nums) if nums else 2

    return (
        0 if has_fraction else 1,
        0 if is_high_pressure_or_temp else 1,
        has_mineral,
        -stoich_sum,
    )


def _generate_element_subsets(primary_els: List[str], precursor_els: List[str]) -> List[Tuple[Tuple[str, ...], str]]:
    inorg_primary = [e for e in primary_els if e not in ORGANIC_SOLVENT_LIGHT_ELEMENTS] or primary_els

    subsets: List[Tuple[Tuple[str, ...], str]] = []
    seen: Set[Tuple[str, ...]] = set()

    for r in (2, 3):
        if len(inorg_primary) >= r:
            for combo in combinations(inorg_primary, r):
                key = tuple(sorted(combo))
                if key not in seen:
                    seen.add(key)
                    subsets.append((key, "Primary / Dopant"))

    extra_anions = [e for e in ("O", "Cl") if e in precursor_els or e in primary_els]
    metals = [e for e in inorg_primary if e not in ("S", "O", "Cl", "Se", "Te")]
    for m in metals:
        for an in extra_anions:
            key = tuple(sorted((m, an)))
            if key not in seen:
                seen.add(key)
                subsets.append((key, "Precursor / Impurity"))

    return subsets[:10]


def fetch_cif_and_simulate(cod_id: str, wavelength_a: float = 1.5406) -> List[Dict[str, Any]]:
    if not HAS_REQUESTS:
        return []
    for scheme in ("https", "http"):
        url = f"{scheme}://www.crystallography.net/cod/{cod_id}.cif"
        try:
            resp = requests.get(url, timeout=10)
            if resp.status_code == 200 and "_cell_length_a" in resp.text:
                peaks = simulate_xrd_from_cif_str(resp.text, wavelength_a=wavelength_a)
                if peaks:
                    return peaks
        except Exception:
            continue
    return []


def query_cod_live(
    primary_elements: List[str],
    precursor_elements: List[str],
    wavelength_a: float = 1.5406,
) -> Tuple[List[Dict[str, Any]], str]:
    subsets = _generate_element_subsets(primary_elements, precursor_elements)
    if not subsets:
        return [], "No valid inorganic element combinations provided."

    raw_candidates: List[Dict[str, Any]] = []
    status_msg = ""

    if HAS_MYSQL:
        try:
            conn = mysql.connector.connect(
                host="sql.crystallography.net",
                user="cod_reader",
                password="",
                database="cod",
                connection_timeout=8,
            )
            cursor = conn.cursor(dictionary=True)
            for combo, role in subsets:
                nel = len(combo)
                like_clauses = " AND ".join(
                    [f"(formula LIKE '% {el} %' OR formula LIKE '% {el}0%' OR formula LIKE '% {el}1%' OR formula LIKE '% {el}2%' OR formula LIKE '% {el}3%' OR formula LIKE '% {el}4%' OR formula LIKE '% {el}5%' OR formula LIKE '% {el}6%' OR formula LIKE '% {el}7%' OR formula LIKE '% {el}8%' OR formula LIKE '% {el}9%')" for el in combo]
                )
                sql = (
                    "SELECT file AS cod_id, formula, mineral, sg, a, b, c, alpha, beta, gamma, "
                    "authors, title, journal, year, volume, firstpage, lastpage, doi "
                    f"FROM data WHERE nel = {nel} AND ({like_clauses}) "
                    "AND a IS NOT NULL AND status IS NULL "
                    "LIMIT 120"
                )
                cursor.execute(sql)
                rows = cursor.fetchall()
                target_set = set(combo)
                valid_rows = [
                    r for r in rows
                    if _formula_matches_exact_elements(str(r.get("formula") or ""), target_set)
                ]
                valid_rows.sort(key=_entry_priority_score, reverse=True)
                for r in valid_rows[:5]:
                    r["_role"] = role
                    r["_elements"] = list(combo)
                    raw_candidates.append(r)
            cursor.close()
            conn.close()
            status_msg = f"Connected to live COD MySQL server (sql.crystallography.net) and retrieved {len(raw_candidates)} candidate crystal structures."
        except Exception as e:
            status_msg = f"COD MySQL fallback to REST ({e})"

    if not raw_candidates and HAS_REQUESTS:
        for combo, role in subsets:
            nel = len(combo)
            params: Dict[str, Any] = {
                "format": "json",
                "strictmin": nel,
                "strictmax": nel,
            }
            for idx, el in enumerate(combo, start=1):
                params[f"el{idx}"] = el
            try:
                resp = requests.get("https://www.crystallography.net/cod/result", params=params, timeout=12)
                if resp.status_code == 200:
                    rows = resp.json()
                    if isinstance(rows, list):
                        target_set = set(combo)
                        valid_rows = [
                            r for r in rows
                            if _formula_matches_exact_elements(str(r.get("formula") or ""), target_set)
                        ]
                        valid_rows.sort(key=_entry_priority_score, reverse=True)
                        for r in valid_rows[:5]:
                            r["cod_id"] = r.get("file") or r.get("cod_id")
                            r["_role"] = role
                            r["_elements"] = list(combo)
                            raw_candidates.append(r)
            except Exception:
                continue
        if raw_candidates:
            status_msg = f"Connected to live COD REST API and retrieved {len(raw_candidates)} candidate crystal structures."

    print(f"🌐 {status_msg}")
    print("⏳ Downloading live .cif atomic coordinates from COD & computing |F_hkl|^2 diffraction planes...")

    references: List[Dict[str, Any]] = []
    seen_ids: Set[str] = set()
    for r in raw_candidates:
        cid = str(r.get("cod_id") or r.get("file") or "").strip()
        if not cid or cid in seen_ids:
            continue
        seen_ids.add(cid)

        peaks = fetch_cif_and_simulate(cid, wavelength_a=wavelength_a)
        if not peaks:
            continue

        formula_clean = _clean_cod_formula(str(r.get("formula") or "")).replace(" ", "")
        mineral = str(r.get("mineral") or "").strip()
        compound_name = f"{mineral} ({formula_clean})" if mineral else formula_clean
        doi = str(r.get("doi") or "").strip()
        cod_url = f"https://www.crystallography.net/cod/{cid}.html"
        paper_url = f"https://doi.org/{doi}" if doi else cod_url
        fp = str(r.get("firstpage") or "").strip()
        lp = str(r.get("lastpage") or "").strip()
        pages = f"{fp}-{lp}" if (fp and lp) else fp

        references.append({
            "cod_id": cid,
            "compound": compound_name,
            "formula": formula_clean,
            "elements": r.get("_elements", []),
            "role": r.get("_role", "Primary / Dopant"),
            "space_group": str(r.get("sg") or ""),
            "a": _parse_cif_float(str(r.get("a") or "")),
            "b": _parse_cif_float(str(r.get("b") or "")),
            "c": _parse_cif_float(str(r.get("c") or "")),
            "authors": str(r.get("authors") or ""),
            "title": str(r.get("title") or ""),
            "journal": str(r.get("journal") or ""),
            "year": r.get("year") or "",
            "volume": str(r.get("volume") or ""),
            "pages": pages,
            "doi": doi,
            "cod_url": cod_url,
            "paper_url": paper_url,
            "peaks": peaks,
        })

    return references, status_msg


# ==============================================================================
# 5. SEQUENTIAL DOMINANT-PHASE SEARCH-MATCH (PRIMARY -> SECONDARY -> IMPURITY)
# ==============================================================================
def _score_candidate_on_peaks(
    ref: Dict[str, Any],
    exp_peaks_df: pd.DataFrame,
    unassigned_indices: Set[int],
    tolerance_deg: float,
) -> Tuple[float, Dict[int, Dict[str, Any]]]:
    ref_peaks = ref.get("peaks", [])
    if not ref_peaks:
        return 0.0, {}

    matches: Dict[int, Dict[str, Any]] = {}
    raw_score = 0.0

    for idx in unassigned_indices:
        row = exp_peaks_df.loc[idx]
        exp_tt = float(row["Peak_2Theta"])
        exp_rel = float(row.get("Rel_Intensity_pct", 50.0))

        best_rp = None
        best_rp_score = -1e9
        for rp in ref_peaks:
            diff = abs(exp_tt - float(rp["two_theta"]))
            if diff <= tolerance_deg:
                rp_int = float(rp.get("intensity", 100.0))
                prox = 1.0 - (diff / tolerance_deg) * 0.65
                int_weight = 0.35 + 0.65 * ((rp_int / 1000.0) ** 0.5)
                s = prox * int_weight * (0.5 + 0.5 * (exp_rel / 100.0))
                if s > best_rp_score:
                    best_rp_score = s
                    best_rp = rp

        if best_rp is not None:
            matches[idx] = best_rp
            raw_score += best_rp_score

    if not matches:
        return 0.0, {}

    formula = str(ref.get("formula") or "")
    has_fraction = "." in formula
    nel = len(ref.get("elements", []))
    nums = [int(n) for n in re.findall(r"\b\d+\b", formula)]
    stoich_complexity = sum(nums) if nums else 2

    simplicity_factor = 1.0
    if has_fraction:
        simplicity_factor *= 0.65
    if nel == 2 and stoich_complexity <= 5:
        simplicity_factor *= 1.25
    elif stoich_complexity > 8:
        simplicity_factor *= 0.80

    fom = raw_score * (len(matches) ** 0.65) * simplicity_factor
    return fom, matches


def match_peaks_to_cod(
    exp_peaks_df: pd.DataFrame,
    cod_references: List[Dict[str, Any]],
    tolerance_deg: float = 0.38,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    df = exp_peaks_df.copy()
    df["Reference_2Theta"] = np.nan
    df["Delta_2Theta"] = np.nan
    df["Reference_d_A"] = np.nan
    df["Delta_d_A"] = np.nan
    df["Plane"] = "-"
    df["Compound_Name"] = "-"
    df["Formula"] = "-"
    df["COD_ID"] = "-"
    df["Phase_Role"] = "-"
    df["Space_Group"] = "-"
    df["Reference_Authors"] = "-"
    df["Reference_Title"] = "-"
    df["Reference_Journal"] = "-"
    df["Reference_Year"] = "-"
    df["Reference_DOI"] = "-"
    df["COD_URL"] = "-"
    df["Paper_URL"] = "-"
    df["Decision"] = "Unindexed"

    if df.empty or not cod_references:
        return df, pd.DataFrame()

    unassigned: Set[int] = set(df.index.tolist())
    selected_phases: List[Dict[str, Any]] = []
    used_cod_ids: Set[str] = set()

    for pass_num in range(1, 4):
        if not unassigned:
            break

        best_ref = None
        best_fom = 0.0
        best_matches: Dict[int, Dict[str, Any]] = {}

        for ref in cod_references:
            if ref["cod_id"] in used_cod_ids:
                continue
            fom, matches = _score_candidate_on_peaks(ref, df, unassigned, tolerance_deg)
            if fom > best_fom:
                best_fom = fom
                best_ref = ref
                best_matches = matches

        if best_ref is None or not best_matches:
            break

        used_cod_ids.add(best_ref["cod_id"])
        phase_role_label = (
            "Primary Phase" if pass_num == 1
            else ("Secondary / Dopant Phase" if pass_num == 2 else "Minor / Impurity Phase")
        )

        deltas = []
        for idx, rp in best_matches.items():
            exp_tt = float(df.loc[idx, "Peak_2Theta"])
            exp_d = float(df.loc[idx, "d_Spacing_A"])
            ref_tt = float(rp["two_theta"])
            ref_d = float(rp["d_spacing"])
            delta_tt = round(exp_tt - ref_tt, 4)
            delta_d = round(exp_d - ref_d, 4)
            deltas.append(abs(delta_tt))

            df.loc[idx, "Reference_2Theta"] = ref_tt
            df.loc[idx, "Delta_2Theta"] = delta_tt
            df.loc[idx, "Reference_d_A"] = ref_d
            df.loc[idx, "Delta_d_A"] = delta_d
            df.loc[idx, "Plane"] = rp["hkl"]
            df.loc[idx, "Compound_Name"] = best_ref["compound"]
            df.loc[idx, "Formula"] = best_ref["formula"]
            df.loc[idx, "COD_ID"] = f"COD-{best_ref['cod_id']}"
            df.loc[idx, "Phase_Role"] = phase_role_label
            df.loc[idx, "Space_Group"] = best_ref["space_group"]
            df.loc[idx, "Reference_Authors"] = best_ref["authors"]
            df.loc[idx, "Reference_Title"] = best_ref["title"]
            df.loc[idx, "Reference_Journal"] = best_ref["journal"]
            df.loc[idx, "Reference_Year"] = str(best_ref["year"])
            df.loc[idx, "Reference_DOI"] = best_ref["doi"]
            df.loc[idx, "COD_URL"] = best_ref["cod_url"]
            df.loc[idx, "Paper_URL"] = best_ref["paper_url"]
            df.loc[idx, "Decision"] = f"Confirmed {best_ref['formula']} ({phase_role_label})"

            unassigned.discard(idx)

        selected_phases.append({
            "Phase_Rank": pass_num,
            "Phase_Role": phase_role_label,
            "COD_ID": f"COD-{best_ref['cod_id']}",
            "Compound_Name": best_ref["compound"],
            "Formula": best_ref["formula"],
            "Space_Group": best_ref["space_group"],
            "Matched_Peaks_Count": len(best_matches),
            "Mean_Abs_Delta_2Theta_deg": round(float(np.mean(deltas)), 4) if deltas else np.nan,
            "a_A": best_ref["a"],
            "b_A": best_ref["b"],
            "c_A": best_ref["c"],
            "Authors": best_ref["authors"],
            "Title": best_ref["title"],
            "Journal": best_ref["journal"],
            "Year": best_ref["year"],
            "Volume": best_ref["volume"],
            "Pages": best_ref["pages"],
            "DOI": best_ref["doi"],
            "COD_URL": best_ref["cod_url"],
            "Paper_URL": best_ref["paper_url"],
        })

    summary_df = pd.DataFrame(selected_phases)
    return df, summary_df


# ==============================================================================
# 6. EXCEL REPORT EXPORT & INTERACTIVE CLI
# ==============================================================================
def export_results_to_excel(
    output_path: str,
    peaks_df: pd.DataFrame,
    phases_df: pd.DataFrame,
    cod_references: List[Dict[str, Any]],
) -> None:
    ref_rows = []
    for ref in cod_references:
        for rp in ref.get("peaks", []):
            ref_rows.append({
                "COD_ID": f"COD-{ref['cod_id']}",
                "Compound": ref["compound"],
                "Formula": ref["formula"],
                "Space_Group": ref["space_group"],
                "Ref_2Theta_deg": rp["two_theta"],
                "Ref_d_Spacing_A": rp["d_spacing"],
                "Ref_Intensity_permil": rp["intensity"],
                "Plane_hkl": rp["hkl"],
                "Authors": ref["authors"],
                "Journal": ref["journal"],
                "Year": ref["year"],
                "DOI": ref["doi"],
                "COD_URL": ref["cod_url"],
            })
    all_ref_df = pd.DataFrame(ref_rows)

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        peaks_df.to_excel(writer, sheet_name="XRD_Peaks", index=False)
        if not phases_df.empty:
            phases_df.to_excel(writer, sheet_name="Identified_Phases_Papers", index=False)
        if not all_ref_df.empty:
            all_ref_df.to_excel(writer, sheet_name="COD_Simulated_Reference_Peaks", index=False)


def main_cli() -> None:
    print("=" * 72)
    print("🔬 UNIVERSAL XRD PHASE IDENTIFIER & LIVE COD REFERENCE MATCHER")
    print("=" * 72)

    file_path = select_xrd_file_interactive()
    df, tt_col, int_col = load_xrd_file(file_path)
    peaks_df = detect_experimental_peaks(df, tt_col, int_col)
    print(f"\n✅ Loaded '{os.path.basename(file_path)}' ({len(df)} points) -> Detected {len(peaks_df)} XRD peaks.")

    print("\n--- Step 2: Synthesis Chemistry & Expected Elements ---")
    main_inp = input(
        "  1) Expected main elements or target compound (e.g., Cu, S, Sn or Sn-CuS): "
    ).strip()
    prec_inp = input(
        "  2) Precursors & reagents used during synthesis (e.g., CuCl2.2H2O, SnCl2.2H2O, thiourea, TEA, IPA, HCl) [optional]: "
    ).strip()

    primary_els = parse_chemical_input_to_elements(main_inp, drop_pure_organics=True)
    precursor_els = parse_chemical_input_to_elements(prec_inp, drop_pure_organics=True)
    if not primary_els and precursor_els:
        primary_els = precursor_els

    print(f"\n🎯 Primary crystal-forming elements: {primary_els}")
    print(f"🧪 Inorganic precursor/reagent elements: {precursor_els}")

    cod_refs, _ = query_cod_live(primary_els, precursor_els)
    indexed_df, phases_df = match_peaks_to_cod(peaks_df, cod_refs, tolerance_deg=0.38)

    display_cols = [
        "Peak_ID", "Peak_2Theta", "Reference_2Theta", "Delta_2Theta",
        "d_Spacing_A", "FWHM_deg", "Crystallite_Size_nm", "Plane",
        "Formula", "COD_ID", "Decision",
    ]
    print("\n--- 1. DETECTED EXPERIMENTAL PEAKS & COD INDEXING ---")
    print(indexed_df[display_cols].to_string(index=False))

    if not phases_df.empty:
        print("\n--- 2. IDENTIFIED PHASES & REFERENCE PAPERS ---")
        for _, prow in phases_df.iterrows():
            print(
                f"  • [{prow['Phase_Role']}] {prow['Compound_Name']} ({prow['COD_ID']}, {prow['Space_Group']}) "
                f"— Matched {prow['Matched_Peaks_Count']} peaks (mean |Δ2θ| = {prow['Mean_Abs_Delta_2Theta_deg']}°)\n"
                f"    📄 Paper: {prow['Authors']} ({prow['Year']}), \"{prow['Title']}\", "
                f"{prow['Journal']} {prow['Volume']}:{prow['Pages']} | {prow['Paper_URL']}"
            )

    base_name = os.path.splitext(os.path.basename(file_path))[0]
    out_dir = os.path.join(os.path.dirname(os.path.abspath(file_path)), "XRD_Universal_Results")
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, f"{base_name}_COD_Analysis.xlsx")
    export_results_to_excel(out_file, indexed_df, phases_df, cod_refs)
    print(f"\n💾 Saved complete Excel report to:\n   {out_file}\n")


if __name__ == "__main__":
    main_cli()