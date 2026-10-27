import re
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.signal import find_peaks, savgol_filter
from scipy.ndimage import uniform_filter1d

# Optional: Direct COD .cif structure factor calculation if pymatgen is installed
try:
    from pymatgen.core import Structure
    from pymatgen.analysis.diffraction.xrd import XRDCalculator
    HAS_PYMATGEN = True
except ImportError:
    HAS_PYMATGEN = False

# ============================================================
# 1. FILE LOCATIONS & SETTINGS
# ============================================================

raw_data_folder = Path(
    r"C:\Users\user\OneDrive\Desktop\Research_Analytics_Project"
    r"\01_Raw_Data\XRD"
)

reference_folder = Path(
    r"C:\Users\user\OneDrive\Desktop\Research_Analytics_Project"
    r"\03_Reference_Data"
)

output_folder = Path(
    r"C:\Users\user\OneDrive\Desktop\Research_Analytics_Project"
    r"\06_Python\XRD_Universal_Results"
)

# Analysis Parameters
min_peak_distance_deg = 0.15
min_snr = 3.0
tolerance = 0.3
minimum_reference_relative_intensity = 5.0


# ============================================================
# 2. AUTOMATIC PHASE & METADATA EXTRACTION FUNCTIONS
# ============================================================

# Regex matching valid chemical formulas (e.g., CuS, Sn2S3, Cu39S28, Cu1.9375S, SnCl2)
FORMULA_TOKEN_REGEX = re.compile(
    r"^(?:[A-Z][a-z]?\d*(?:\.\d+)?)+$"
)

# Common words that look like element symbols if capitalized
NON_FORMULA_WORDS = {
    "PDF", "COD", "ICDD", "ICSD", "CARD", "PHASE", "RAW", "XRD",
    "CSV", "DATA", "FILE", "REF", "INT", "HKL", "SYN", "ALPHA", "BETA"
}


def normalise_card_number(value):
    return re.sub(r"[^0-9]", "", str(value))


def is_chemical_formula(token):
    """Returns True if a string token looks like a valid chemical formula."""
    clean = token.strip(" ()[]{}\"',;:")
    if not clean or clean.upper() in NON_FORMULA_WORDS:
        return False
    # Must match element + stoichiometry pattern and contain at least one uppercase letter
    if not FORMULA_TOKEN_REGEX.match(clean):
        return False
    # Reject normal English words (e.g., "Covellite") that have 2+ lowercase letters in a row
    if re.search(r"[a-z]{2,}", clean):
        return False
    return True


def check_contains_sn(formula_str, phase_str=""):
    """Detects whether the phase contains Tin (Sn) without false-matching Sm, Se, Si, etc."""
    combined = f"{formula_str} {phase_str}"
    return bool(re.search(r"Sn(?![a-z])", combined)) or ("tin" in combined.lower())


def assign_dynamic_priority(formula_str):
    """
    Automatically assigns matching priority based on crystal chemistry:
      1 -> Primary host phase (CuS)
      2 -> Stoichiometric binary sulfides (e.g., Sn2S3, SnS2, SnS, Cu2S)
      3 -> Non-stoichiometric / supercell sulfides (e.g., Cu39S28, Cu1.9375S)
      4 -> Halides / secondary salts / others (e.g., CuCl, SnCl2)
    """
    clean_f = formula_str.replace(" ", "")
    if clean_f == "CuS":
        return 1
    if "Cl" in clean_f or "O" in clean_f:
        return 4
    if "." in clean_f or any(int(n) > 10 for n in re.findall(r"\d+", clean_f)):
        return 3
    return 2


def read_text_file(file_path):
    encodings = ["utf-8-sig", "utf-8", "cp1252", "latin1"]
    for encoding in encodings:
        try:
            with open(file_path, "r", encoding=encoding) as file:
                return file.readlines()
        except UnicodeDecodeError:
            continue
    raise ValueError(f"Could not read reference file: {file_path}")


def extract_metadata_from_reference(file_path, header_lines):
    """
    Automatically extracts Card Number, Formula, Phase, and Mineral Name
    from CSV header lines and/or filename without any hardcoded registry.
    """
    card_number = normalise_card_number(file_path.stem) or file_path.stem
    formula = None
    mineral_name = None
    phase_name = None

    # 1. Scan inside the file header lines for labeled or unlabeled metadata
    for line in header_lines:
        clean = line.strip().strip('",;')
        if not clean:
            continue
        lower = clean.lower()

        # Check key-value pairs separated by ':' or ',' or '=' or '\t'
        if any(sep in clean for sep in (":", ",", "=", "\t")):
            parts = re.split(r"[:=\t,]+", clean, maxsplit=1)
            if len(parts) == 2:
                key = parts[0].strip(' "#').lower()
                val = parts[1].strip(' ",;\t')

                if val:
                    # Card / Entry ID
                    if any(k in key for k in ("pdf", "card", "cod", "icsd", "entry", "ref_code")):
                        extracted_num = normalise_card_number(val)
                        if extracted_num:
                            card_number = extracted_num

                    # Chemical Formula
                    elif any(k in key for k in ("formula", "chem", "empirical")):
                        formula = re.sub(r"\s+", "", val)

                    # Mineral / Compound / Phase Name
                    elif any(k in key for k in ("mineral", "compound", "name", "phase", "material")):
                        # If the field actually holds a formula, assign to formula
                        if is_chemical_formula(val.replace(" ", "")) and not formula:
                            formula = val.replace(" ", "")
                        else:
                            mineral_name = val

        # Also scan free-text comment/header lines for tokens that match chemical formulas
        tokens = re.split(r"[\s,;()\[\]\"']+", clean)
        for tok in tokens:
            if not formula and is_chemical_formula(tok):
                formula = tok

    # 2. Fallback: Extract from filename (e.g., "000752233_CuS_Covellite.csv" or "COD_9008386_CuS.csv")
    stem_tokens = [
        t for t in re.split(r"[_\-\s]+", file_path.stem)
        if t and not t.isdigit() and t.upper() not in NON_FORMULA_WORDS
    ]
    for tok in stem_tokens:
        if not formula and is_chemical_formula(tok):
            formula = tok
        elif not mineral_name and not is_chemical_formula(tok) and len(tok) > 2:
            mineral_name = tok

    # 3. Consolidate final fields
    if not formula and phase_name and is_chemical_formula(phase_name):
        formula = phase_name
    if not formula:
        formula = mineral_name if mineral_name else file_path.stem

    phase_name = formula
    if not mineral_name:
        mineral_name = phase_name

    contains_sn = check_contains_sn(formula, mineral_name)
    priority = assign_dynamic_priority(formula)

    return {
        "Reference_Card": card_number,
        "Phase": phase_name,
        "Mineral_Name": mineral_name,
        "Formula": formula,
        "Contains_Sn": contains_sn,
        "Priority": priority
    }


def parse_reference_csv(file_path):
    """
    Parses any reference card CSV file (ICDD, COD, HighScore, Jade, etc.)
    and dynamically identifies Phase, Formula, and Mineral Name from the file itself.
    """
    lines = read_text_file(file_path)

    # Locate the start of the numerical diffraction peak list
    diffraction_start = None
    diffraction_pattern = re.compile(
        r"^\s*"
        r"([+-]?\d+(?:\.\d+)?)"          # 2theta
        r"[\s,\t;]+"
        r"([+-]?\d+(?:\.\d+)?)"          # d-spacing
        r"[\s,\t;]+"
        r"([+-]?\d+(?:\.\d+)?)"          # intensity
        r"(?:[\s,\t;]+(.+?))?"           # hkl (optional if 3-column CSV)
        r"\s*$"
    )

    for index, line in enumerate(lines):
        cleaned_line = line.strip().replace('"', "")
        lower_line = cleaned_line.lower()

        if "2theta" in lower_line and ("int" in lower_line or "d" in lower_line):
            diffraction_start = index + 1
            break
        if diffraction_pattern.match(cleaned_line):
            diffraction_start = index
            break

    if diffraction_start is None:
        diffraction_start = 0

    # Extract phase identity from all header lines above the diffraction table
    meta = extract_metadata_from_reference(file_path, lines[:diffraction_start])

    extracted_rows = []
    for line in lines[diffraction_start:]:
        cleaned_line = line.strip().replace('"', "")
        if not cleaned_line:
            continue
        match = diffraction_pattern.match(cleaned_line)
        if match is None:
            continue

        raw_hkl = match.group(4) if match.group(4) else "-"
        hkl = re.sub(r"[\t,;]+", " ", raw_hkl.strip())
        hkl = re.sub(r"\s+", " ", hkl)

        extracted_rows.append({
            "Reference_Card": meta["Reference_Card"],
            "Phase": meta["Phase"],
            "Mineral_Name": meta["Mineral_Name"],
            "Formula": meta["Formula"],
            "Contains_Sn": meta["Contains_Sn"],
            "Phase_Priority": meta["Priority"],
            "Reference_2theta": float(match.group(1)),
            "Reference_d": float(match.group(2)),
            "Reference_Intensity": float(match.group(3)),
            "HKL": hkl,
            "Reference_File": file_path.name
        })

    if not extracted_rows:
        return None

    ref_df = pd.DataFrame(extracted_rows)

    # Merge symmetry-equivalent reflections at identical 2-theta
    ref_df["2theta_bin"] = ref_df["Reference_2theta"].round(2)
    ref_df = (
        ref_df.sort_values("Reference_Intensity", ascending=False)
        .groupby(["Reference_Card", "2theta_bin"], as_index=False)
        .agg({
            "Phase": "first",
            "Mineral_Name": "first",
            "Formula": "first",
            "Contains_Sn": "first",
            "Phase_Priority": "first",
            "Reference_2theta": "first",
            "Reference_d": "first",
            "Reference_Intensity": "sum",
            "HKL": lambda s: " / ".join(dict.fromkeys(s)),
            "Reference_File": "first"
        })
        .drop(columns=["2theta_bin"])
        .sort_values("Reference_2theta")
        .reset_index(drop=True)
    )

    max_int = ref_df["Reference_Intensity"].max()
    ref_df["Reference_Relative_Intensity"] = (
        (ref_df["Reference_Intensity"] / max_int * 100.0) if max_int > 0 else 0.0
    )

    return ref_df


def parse_cod_cif_file(file_path):
    """
    Directly reads a COD .cif file, extracts crystal metadata, and calculates
    the Cu K-alpha XRD reference pattern if pymatgen is installed.
    """
    if not HAS_PYMATGEN:
        print(f"  - Note: Skipping {file_path.name} (install 'pymatgen' to read .cif files directly)")
        return None

    lines = read_text_file(file_path)
    cod_id = normalise_card_number(file_path.stem) or file_path.stem
    formula = None
    mineral_name = None

    for line in lines:
        s = line.strip()
        if s.startswith("_chemical_formula_sum") or s.startswith("_chemical_formula_structural"):
            parts = s.split(maxsplit=1)
            if len(parts) > 1:
                formula = parts[1].strip(" '\"").replace(" ", "")
        elif s.startswith("_chemical_name_mineral") or s.startswith("_chemical_name_common") or s.startswith("_chemical_name_systematic"):
            parts = s.split(maxsplit=1)
            if len(parts) > 1 and parts[1].strip(" ?'\"."):
                mineral_name = parts[1].strip(" '\"")
        elif s.startswith("_cod_database_code") or s.startswith("_database_code_ICSD"):
            parts = s.split(maxsplit=1)
            if len(parts) > 1:
                cod_id = parts[1].strip(" '\"")

    structure = Structure.from_file(str(file_path))

    # If CIF tags omitted the formula, compute reduced formula directly from the unit cell structure
    if not formula:
        formula = structure.composition.reduced_formula
    if not mineral_name:
        mineral_name = formula

    contains_sn = check_contains_sn(formula, mineral_name)
    priority = assign_dynamic_priority(formula)

    calc = XRDCalculator(wavelength="CuKa")
    pattern = calc.get_pattern(structure, scaled=True, two_theta_range=(5, 90))

    rows = []
    for two_theta, intensity, hkls, d_val in zip(pattern.x, pattern.y, pattern.hkls, pattern.d_hkls):
        hkl_list = [" ".join(str(i) for i in h.get("hkl", ())) for h in hkls if h.get("hkl")]
        rows.append({
            "Reference_Card": cod_id,
            "Phase": formula,
            "Mineral_Name": mineral_name,
            "Formula": formula,
            "Contains_Sn": contains_sn,
            "Phase_Priority": priority,
            "Reference_2theta": float(two_theta),
            "Reference_d": float(d_val),
            "Reference_Intensity": float(intensity),
            "Reference_Relative_Intensity": float(intensity),
            "HKL": " / ".join(hkl_list) if hkl_list else "-",
            "Reference_File": file_path.name
        })

    return pd.DataFrame(rows) if rows else None


# ============================================================
# 3. LOAD ALL AVAILABLE REFERENCE CARDS (.CSV AND .CIF)
# ============================================================

print(f"\nScanning reference folder: {reference_folder}")
reference_tables = []

for ref_file in sorted(reference_folder.glob("*.csv")):
    try:
        df_card = parse_reference_csv(ref_file)
        if df_card is not None:
            reference_tables.append(df_card)
            print(
                f"  - Loaded CSV Card {df_card['Reference_Card'].iloc[0]}: "
                f"Phase={df_card['Phase'].iloc[0]} | "
                f"Mineral={df_card['Mineral_Name'].iloc[0]} | "
                f"Sn={df_card['Contains_Sn'].iloc[0]} | "
                f"Priority={df_card['Phase_Priority'].iloc[0]}"
            )
    except Exception as e:
        print(f"  - Warning: Could not parse {ref_file.name}: {e}")

for cif_file in sorted(reference_folder.glob("*.cif")):
    try:
        df_cif = parse_cod_cif_file(cif_file)
        if df_cif is not None:
            reference_tables.append(df_cif)
            print(
                f"  - Loaded COD CIF  {df_cif['Reference_Card'].iloc[0]}: "
                f"Phase={df_cif['Phase'].iloc[0]} | "
                f"Mineral={df_cif['Mineral_Name'].iloc[0]} | "
                f"Sn={df_cif['Contains_Sn'].iloc[0]} | "
                f"Priority={df_cif['Phase_Priority'].iloc[0]}"
            )
    except Exception as e:
        print(f"  - Warning: Could not parse {cif_file.name}: {e}")

if not reference_tables:
    raise ValueError("No valid reference cards found in the reference folder!")

all_reference_peaks = pd.concat(reference_tables, ignore_index=True)


# ============================================================
# 4. BATCH PROCESS ALL RAW XRD FILES & IDENTIFY CRYSTALS
# ============================================================

raw_files = sorted(list(raw_data_folder.glob("*.xlsx")) + list(raw_data_folder.glob("*.xls")))

if not raw_files:
    raise FileNotFoundError(f"No raw XRD Excel files found in: {raw_data_folder}")

print(f"\nFound {len(raw_files)} raw XRD file(s) for universal crystal identification.\n")

for raw_file_path in raw_files:
    print("=" * 75)
    print(f"Processing File: {raw_file_path.name}")
    print("=" * 75)

    raw_df = pd.read_excel(raw_file_path)
    required_columns = ["X", "d (Å)", "Y (counts)"]

    if not all(col in raw_df.columns for col in required_columns):
        print(f"Skipping {raw_file_path.name}: Missing required columns {required_columns}")
        continue

    raw_df = raw_df.dropna(subset=required_columns).copy()
    x = raw_df["X"].to_numpy()
    y = raw_df["Y (counts)"].to_numpy()

    # --- 4A. Chemistry Filter Based on Sample Filename ---
    clean_stem = raw_file_path.stem.lower().replace(" ", "").replace("_", "-")
    if clean_stem.startswith("cus") and "sn" not in clean_stem:
        sample_refs = all_reference_peaks[~all_reference_peaks["Contains_Sn"]].copy()
        sample_mode = "Pure CuS (Sn reference phases excluded)"
    else:
        sample_refs = all_reference_peaks.copy()
        sample_mode = "Sn-CuS / Universal (All reference phases enabled)"

    print(f"Sample Mode: {sample_mode}")

    # --- 4B. Peak Detection Dynamics (IUPAC 5-sigma + Baseline-Aware) ---
    step_size = float(np.median(np.diff(x)))
    dynamic_distance = max(1, int(round(min_peak_distance_deg / step_size)))

    sg_window = int(round(0.25 / step_size))
    if sg_window % 2 == 0:
        sg_window += 1
    sg_window = max(7, sg_window)

    smoothed_y = savgol_filter(y, window_length=sg_window, polyorder=3)
    baseline_pts = max(15, int(round(2.5 / step_size)))
    baseline = uniform_filter1d(smoothed_y, size=baseline_pts, mode="nearest")
    net_intensity = np.maximum(0.0, smoothed_y - baseline)

    noise_residuals = y - smoothed_y
    sigma_noise = float(np.std(noise_residuals))
    prominence_cutoff = 5.0 * sigma_noise
    relative_floor = 0.025 * float(np.max(net_intensity))
    dynamic_prominence = max(prominence_cutoff, relative_floor)

    peak_indices, peak_properties = find_peaks(
        smoothed_y,
        distance=dynamic_distance,
        prominence=dynamic_prominence,
        width=3
    )

    detected_peaks = raw_df.iloc[peak_indices][["X", "d (Å)", "Y (counts)"]].copy()
    detected_peaks = detected_peaks.rename(columns={
        "X": "Detected_2theta",
        "d (Å)": "Detected_d",
        "Y (counts)": "Experimental_Intensity"
    }).reset_index(drop=True)

    detected_peaks.insert(0, "Peak_ID", range(1, len(detected_peaks) + 1))
    detected_peaks["Prominence"] = peak_properties["prominences"].round(2)
    max_peak_intensity = detected_peaks["Experimental_Intensity"].max()
    detected_peaks["Experimental_Relative_Intensity"] = (
        detected_peaks["Experimental_Intensity"] / max_peak_intensity * 100.0
    ).round(2)

    # --- 4C. Filter Reference Peaks to Scan Range ---
    exp_min, exp_max = x.min(), x.max()
    ref_for_matching = sample_refs[
        (sample_refs["Reference_2theta"] >= exp_min) &
        (sample_refs["Reference_2theta"] <= exp_max) &
        (sample_refs["Reference_Relative_Intensity"] >= minimum_reference_relative_intensity)
    ].copy()

    # --- 4D. Match Experimental Peaks to Reference Crystal Phases ---
    all_matches = []
    unmatched_peaks = []

    for _, det_row in detected_peaks.iterrows():
        det_2theta = det_row["Detected_2theta"]
        candidates = ref_for_matching.copy()
        candidates["Difference"] = (candidates["Reference_2theta"] - det_2theta).abs()
        possible = candidates[candidates["Difference"] <= tolerance].copy()

        if possible.empty:
            nearest = candidates.loc[candidates["Difference"].idxmin()]
            unmatched_peaks.append({
                "Peak_ID": int(det_row["Peak_ID"]),
                "Detected_2theta": det_2theta,
                "Detected_d": det_row["Detected_d"],
                "Experimental_Intensity": det_row["Experimental_Intensity"],
                "Experimental_Relative_Intensity": det_row["Experimental_Relative_Intensity"],
                "Nearest_Phase": nearest["Phase"],
                "Nearest_Mineral": nearest["Mineral_Name"],
                "Nearest_Reference_2theta": nearest["Reference_2theta"],
                "Nearest_Difference": round(nearest["Difference"], 4),
                "Status": "NO MATCH"
            })
            continue

        possible = possible.sort_values(
            by=["Phase_Priority", "Difference", "Reference_Relative_Intensity"],
            ascending=[True, True, False]
        )

        for rank, (_, ref_row) in enumerate(possible.iterrows(), start=1):
            all_matches.append({
                "Peak_ID": int(det_row["Peak_ID"]),
                "Detected_2theta": det_2theta,
                "Detected_d": det_row["Detected_d"],
                "Experimental_Intensity": det_row["Experimental_Intensity"],
                "Experimental_Relative_Intensity": det_row["Experimental_Relative_Intensity"],
                "Phase": ref_row["Phase"],
                "Mineral_Name": ref_row["Mineral_Name"],
                "Formula": ref_row["Formula"],
                "HKL": ref_row["HKL"],
                "Reference_2theta": ref_row["Reference_2theta"],
                "Reference_d": ref_row["Reference_d"],
                "Difference": round(ref_row["Difference"], 4),
                "Reference_Intensity": ref_row["Reference_Intensity"],
                "Reference_Relative_Intensity": round(ref_row["Reference_Relative_Intensity"], 2),
                "Reference_Card": ref_row["Reference_Card"],
                "Phase_Priority": ref_row["Phase_Priority"],
                "Match_Rank_For_Peak": rank,
                "Status": "MATCH"
            })

    all_matches_df = pd.DataFrame(all_matches)
    unmatched_peaks_df = pd.DataFrame(unmatched_peaks)

    # Pick the single best crystal match for each peak
    if not all_matches_df.empty:
        best_matches_df = (
            all_matches_df
            .sort_values(
                by=["Peak_ID", "Phase_Priority", "Difference", "Reference_Relative_Intensity"],
                ascending=[True, True, True, False]
            )
            .drop_duplicates(subset=["Peak_ID"], keep="first")
            .reset_index(drop=True)
        )
    else:
        best_matches_df = pd.DataFrame()

    # --- 4E. Build the Unified Master Table (Peaks + Crystal Identification) ---
    if not best_matches_df.empty:
        crystal_id_df = detected_peaks.merge(
            best_matches_df[[
                "Peak_ID",
                "Phase",
                "Mineral_Name",
                "Formula",
                "HKL",
                "Reference_2theta",
                "Reference_d",
                "Difference",
                "Reference_Relative_Intensity",
                "Reference_Card",
                "Status"
            ]],
            on="Peak_ID",
            how="left"
        )
    else:
        crystal_id_df = detected_peaks.copy()
        for col in ["Phase", "Mineral_Name", "Formula", "HKL", "Reference_Card"]:
            crystal_id_df[col] = "Unidentified"
        for col in ["Reference_2theta", "Reference_d", "Difference", "Reference_Relative_Intensity"]:
            crystal_id_df[col] = np.nan
        crystal_id_df["Status"] = "NO MATCH"

    crystal_id_df["Phase"] = crystal_id_df["Phase"].fillna("Unidentified")
    crystal_id_df["Mineral_Name"] = crystal_id_df["Mineral_Name"].fillna("-")
    crystal_id_df["Formula"] = crystal_id_df["Formula"].fillna("-")
    crystal_id_df["HKL"] = crystal_id_df["HKL"].fillna("-")
    crystal_id_df["Status"] = crystal_id_df["Status"].fillna("NO MATCH")

    # Print the unified Crystal Identification Table directly in the terminal
    print("\nCRYSTAL IDENTIFICATION RESULT:")
    display_cols = [
        "Peak_ID",
        "Detected_2theta",
        "Detected_d",
        "Experimental_Relative_Intensity",
        "Phase",
        "Mineral_Name",
        "HKL",
        "Reference_2theta",
        "Difference",
        "Reference_Card",
        "Status"
    ]
    print(crystal_id_df[display_cols].to_string(index=False))

    # --- 4F. Save Results ---
    out_file = output_folder / f"{raw_file_path.stem}_Universal_Analysis.xlsx"
    output_folder.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(out_file, engine="openpyxl") as writer:
        crystal_id_df.to_excel(writer, sheet_name="Crystal_Identification", index=False)
        if not all_matches_df.empty:
            all_matches_df.to_excel(writer, sheet_name="All_Possible_Matches", index=False)
        if not unmatched_peaks_df.empty:
            unmatched_peaks_df.to_excel(writer, sheet_name="Unmatched_Peaks", index=False)
        detected_peaks.to_excel(writer, sheet_name="Raw_Detected_Peaks", index=False)

    print(f"\n -> Saved to Excel: {out_file.name}\n")

print("All files processed successfully!")