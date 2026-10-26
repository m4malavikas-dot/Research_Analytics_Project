import re
from pathlib import Path

import pandas as pd
from scipy.signal import find_peaks


# ============================================================
# 1. FILE LOCATIONS
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
    r"\06_Python\XRD_Results"
)

# ============================================================
# 2. REFERENCE CARD INFORMATION
# ============================================================

reference_phases = {
    "000752233": {
        "Phase": "CuS",
        "Mineral_Name": "Covellite",
        "Formula": "CuS",
        "Contains_Sn": False
    },

    "000360380": {
        "Phase": "Cu39S28",
        "Mineral_Name": "Spionkopite",
        "Formula": "Cu39S28",
        "Contains_Sn": False
    },

    "000711383": {
        "Phase": "Cu1.9375S",
        "Mineral_Name": "Djurleite",
        "Formula": "Cu1.9375S",
        "Contains_Sn": False
    },

    "000822117": {
        "Phase": "CuCl",
        "Mineral_Name": None,
        "Formula": "CuCl",
        "Contains_Sn": False
    },

    "000752183": {
        "Phase": "Sn2S3",
        "Mineral_Name": None,
        "Formula": "Sn2S3",
        "Contains_Sn": True
    },

    "000890474": {
        "Phase": "SnCl2",
        "Mineral_Name": None,
        "Formula": "SnCl2",
        "Contains_Sn": True
    }
}

# ============================================================
# 3. SETTINGS
# ============================================================

# Physical minimum separation between distinct peaks in degrees 2-theta
# (typical nanocrystalline XRD peak widths are ~0.15° to 0.25°)
min_peak_distance_deg = 0.15

# Minimum signal-to-noise ratio: peak prominence must be >= 3 * sigma_noise
min_snr = 3.0

# Minimum relative prominence as a percentage of the highest peak (1% to 2%)
min_relative_prominence = 1.5

# Maximum allowed difference between experimental and reference 2θ
tolerance = 0.3

# Ignore extremely weak reference peaks below this percentage
minimum_reference_relative_intensity = 5.0

# ============================================================
# 4. HELPER FUNCTIONS
# ============================================================

def normalise_card_number(value):
    """
    Removes spaces and hyphens from a reference-card number.

    Example:
    00-071-1383 becomes 000711383
    """

    return re.sub(r"[^0-9]", "", str(value))


def read_text_file(file_path):
    """
    Reads the reference-card file.

    Several encodings are attempted because exported reference files
    may not always use UTF-8 encoding.
    """

    encodings = ["utf-8-sig", "utf-8", "cp1252", "latin1"]

    for encoding in encodings:
        try:
            with open(file_path, "r", encoding=encoding) as file:
                return file.readlines()
        except UnicodeDecodeError:
            continue

    raise ValueError(f"Could not read reference file: {file_path}")


def find_reference_file(folder, card_number):
    """
    Searches the reference-data folder for a file whose name matches
    the reference-card number.

    It works whether the filename contains hyphens or not.

    Examples:
    00-071-1383.csv
    000711383.csv
    """

    required_number = normalise_card_number(card_number)

    for file_path in folder.glob("*.csv"):
        file_number = normalise_card_number(file_path.stem)

        if file_number == required_number:
            return file_path

    return None


def extract_reference_peaks(file_path, card_number, phase_information):
    """
    Extracts the diffraction data from a reference-card file.

    The function searches for the table beginning with a heading such as:

    2theta    d[A]    Int.    hkl

    It ignores the descriptive material information above the table.
    """

    lines = read_text_file(file_path)

    diffraction_start = None

    # Find the diffraction-data column heading
    for index, line in enumerate(lines):
        cleaned_line = line.lower().replace('"', "")

        if (
            "2theta" in cleaned_line
            and "int" in cleaned_line
            and "hkl" in cleaned_line
        ):
            diffraction_start = index + 1
            break

    if diffraction_start is None:
        raise ValueError(
            f"Could not locate the diffraction-data table in:\n{file_path}"
        )

    extracted_rows = []

    # Pattern used to read:
    # 2theta, d-spacing, intensity and remaining hkl text
    diffraction_pattern = re.compile(
        r"^\s*"
        r"([+-]?\d+(?:\.\d+)?)"          # 2theta
        r"[\s,\t;]+"
        r"([+-]?\d+(?:\.\d+)?)"          # d-spacing
        r"[\s,\t;]+"
        r"([+-]?\d+(?:\.\d+)?)"          # intensity
        r"[\s,\t;]+"
        r"(.+?)"                          # hkl
        r"\s*$"
    )

    for line in lines[diffraction_start:]:
        cleaned_line = line.strip().replace('"', "")

        if not cleaned_line:
            continue

        match = diffraction_pattern.match(cleaned_line)

        if match is None:
            continue

        reference_2theta = float(match.group(1))
        reference_d = float(match.group(2))
        reference_intensity = float(match.group(3))

        hkl = match.group(4).strip()

        # Standardise hkl formatting
        hkl = re.sub(r"[\t,;]+", " ", hkl)
        hkl = re.sub(r"\s+", " ", hkl)

        extracted_rows.append({
            "Reference_Card": card_number,
            "Phase": phase_information["Phase"],
            "Mineral_Name": phase_information["Mineral_Name"],
            "Formula": phase_information["Formula"],
            "Reference_2theta": reference_2theta,
            "Reference_d": reference_d,
            "Reference_Intensity": reference_intensity,
            "HKL": hkl,
            "Reference_File": file_path.name
        })

    if not extracted_rows:
        raise ValueError(
            f"No diffraction peaks were extracted from:\n{file_path}"
        )

    reference_df = pd.DataFrame(extracted_rows)

    # Normalise reference intensity separately for each reference card
    maximum_intensity = reference_df["Reference_Intensity"].max()

    if maximum_intensity > 0:
        reference_df["Reference_Relative_Intensity"] = (
            reference_df["Reference_Intensity"]
            / maximum_intensity
            * 100
        )
    else:
        reference_df["Reference_Relative_Intensity"] = 0

    return reference_df

def choose_raw_xrd_file(folder):
    """
    Shows all Excel files in the raw XRD folder
    and asks the user to choose one.
    """

    raw_files = sorted(
        list(folder.glob("*.xlsx")) +
        list(folder.glob("*.xls"))
    )

    if not raw_files:
        raise FileNotFoundError(
            f"No Excel XRD files were found in:\n{folder}"
        )

    print("\nAVAILABLE RAW XRD FILES")

    for number, file_path in enumerate(raw_files, start=1):
        print(f"{number}. {file_path.name}")

    while True:
        try:
            choice = int(
                input(
                    "\nEnter the number of the raw XRD file "
                    "you want to analyse: "
                )
            )

            if 1 <= choice <= len(raw_files):
                selected_file = raw_files[choice - 1]

                print(
                    f"\nSelected raw file: {selected_file.name}"
                )

                return selected_file

            print(
                f"Please enter a number between "
                f"1 and {len(raw_files)}."
            )

        except ValueError:
            print("Please enter a valid number.")


def detect_sample_type(raw_file_path):
    """
    Automatically determines the sample type from the filename.

    Recognised types:
    - CuS files
    - Sn-CuS files

    Unknown materials are stopped to prevent wrong comparison.
    """

    file_name = raw_file_path.stem.lower()
    clean_file_name = file_name.replace(" ", "").replace("_", "-")

    if (
        clean_file_name.startswith("sn-cus")
        or clean_file_name.startswith("sncus")
    ):
        return "Sn-CuS"

    if clean_file_name.startswith("cus"):
        return "CuS"

    raise ValueError(
        f"\nCould not identify the sample type from the filename: "
        f"{raw_file_path.name}\n"
        f"This program currently recognises only CuS and Sn-CuS files.\n"
        f"Add a new sample-type rule before analysing another material."
    )

def select_reference_phases(reference_dictionary, sample_type):
    """
    Selects reference cards based on sample chemistry.

    CuS:
        Uses only references without Sn.

    Sn-CuS:
        Uses all references.
    """

    if sample_type == "CuS":

        selected_references = {
            card_number: phase_information
            for card_number, phase_information
            in reference_dictionary.items()
            if phase_information["Contains_Sn"] is False
        }

    elif sample_type == "Sn-CuS":

        selected_references = reference_dictionary.copy()

    else:
        raise ValueError(
            f"Unknown sample type: {sample_type}"
        )

    return selected_references

# ============================================================
# 5. READ EXPERIMENTAL XRD DATA
# ============================================================

# Ask the user which raw XRD file should be analysed
raw_file_path = choose_raw_xrd_file(
    raw_data_folder
)

# Automatically detect the sample type from the selected filename
sample_type = detect_sample_type(
    raw_file_path
)

# Select only the chemically relevant reference phases
selected_reference_phases = select_reference_phases(
    reference_phases,
    sample_type
)

print(f"\nSample type detected automatically: {sample_type}")

print("\nREFERENCE PHASES THAT WILL BE USED")
for card_number, phase_information in selected_reference_phases.items():
    print(
        f"{card_number}: "
        f"{phase_information['Phase']} "
        f"({phase_information['Formula']})"
    )

# Create the output filename automatically
safe_sample_type = sample_type.replace("-", "_")
output_file_name = (
    f"{raw_file_path.stem}_"
    f"{safe_sample_type}_XRD_Analysis.xlsx"
)
output_file_path = output_folder / output_file_name

if not raw_file_path.exists():
    raise FileNotFoundError(
        f"Experimental XRD file was not found:\n{raw_file_path}"
    )

raw_df = pd.read_excel(raw_file_path)

required_columns = ["X", "d (Å)", "Y (counts)"]
missing_columns = [
    col for col in required_columns
    if col not in raw_df.columns
]
if missing_columns:
    raise ValueError(
        f"Missing columns in the experimental file: {missing_columns}"
    )

# Remove rows containing missing values in essential columns
raw_df = raw_df.dropna(subset=required_columns).copy()

x = raw_df["X"].to_numpy()
y = raw_df["Y (counts)"].to_numpy()


# ============================================================
# 6. DETECT EXPERIMENTAL PEAKS (IUPAC 5σ LOQ + BASELINE-AWARE)
# ============================================================

import numpy as np
from scipy.signal import find_peaks, savgol_filter
from scipy.ndimage import uniform_filter1d

# 1. Determine instrument step size (Δ2θ)
step_size = float(np.median(np.diff(x)))

# 2. Derive point distance from instrumental resolution floor (0.15°)
# Ensures the distance setting works on scans with any step size
dynamic_distance = max(1, int(round(0.15 / step_size)))

# 3. Savitzky-Golay filter to suppress high-frequency detector noise
# Uses ~0.25° window (~13 points for 0.02° step), polyorder 3
sg_window = int(round(0.25 / step_size))
if sg_window % 2 == 0:
    sg_window += 1
sg_window = max(7, sg_window)

smoothed_y = savgol_filter(y, window_length=sg_window, polyorder=3)

# 4. Background baseline estimation (rolling window ~2.5° wide)
# Subtracting this baseline removes glass/amorphous substrate humps <25°
baseline_pts = max(15, int(round(2.5 / step_size)))
baseline = uniform_filter1d(smoothed_y, size=baseline_pts, mode="nearest")
net_intensity = np.maximum(0.0, smoothed_y - baseline)

# 5. Measure instrument noise floor (σ_noise) from high-frequency residuals
noise_residuals = y - smoothed_y
sigma_noise = float(np.std(noise_residuals))

# 6. Prominence threshold: IUPAC 5σ Limit of Quantitation (LOQ)
# Evaluates against net prominence so baseline curvature doesn't trigger peaks
prominence_cutoff = 5.0 * sigma_noise

# Optional floor: at least 2.5% of the strongest net reflection
max_net_signal = float(np.max(net_intensity))
relative_floor = 0.025 * max_net_signal
dynamic_prominence = max(prominence_cutoff, relative_floor)

print(f"\n[Peak Detection Dynamics]")
print(f"  - Step size: {step_size:.4f}° 2θ")
print(f"  - Calculated point distance: {dynamic_distance} points (0.15°)")
print(f"  - Noise floor (σ): {sigma_noise:.2f} counts")
print(f"  - IUPAC 5σ LOQ threshold: {prominence_cutoff:.2f} counts")
print(f"  - Effective prominence cutoff: {dynamic_prominence:.2f} counts")

# 7. Identify peaks on the smoothed curve using net prominence
peak_indices, peak_properties = find_peaks(
    smoothed_y,
    distance=dynamic_distance,
    prominence=dynamic_prominence,
    width=3
)

detected_peaks = raw_df.iloc[peak_indices][
    ["X", "d (Å)", "Y (counts)"]
].copy()

detected_peaks = detected_peaks.rename(columns={
    "X": "Detected_2theta",
    "d (Å)": "Detected_d",
    "Y (counts)": "Experimental_Intensity"
})

detected_peaks = detected_peaks.reset_index(drop=True)
detected_peaks.insert(0, "Peak_ID", range(1, len(detected_peaks) + 1))
detected_peaks["Prominence"] = peak_properties["prominences"]

# Normalize relative intensity (Base 100 on strongest detected reflection)
max_peak_intensity = detected_peaks["Experimental_Intensity"].max()
detected_peaks["Experimental_Relative_Intensity"] = (
    detected_peaks["Experimental_Intensity"]
    / max_peak_intensity
    * 100.0
)

print(f"  - Total true peaks detected: {len(detected_peaks)}")


# ============================================================
# 7. READ ALL REFERENCE CARDS
# ============================================================

reference_tables = []
missing_reference_files = []

for card_number, phase_information in selected_reference_phases.items():

    reference_file = find_reference_file(
        reference_folder,
        card_number
    )

    if reference_file is None:
        missing_reference_files.append(card_number)

        print(
            f"WARNING: Reference file not found for card "
            f"{card_number}"
        )

        continue

    print(
        f"Reading {phase_information['Phase']} "
        f"from {reference_file.name}"
    )

    card_df = extract_reference_peaks(
        file_path=reference_file,
        card_number=card_number,
        phase_information=phase_information
    )

    reference_tables.append(card_df)


if not reference_tables:
    raise ValueError(
        "No reference-card diffraction data could be loaded."
    )

all_reference_peaks = pd.concat(
    reference_tables,
    ignore_index=True
)


# ============================================================
# 8. FILTER REFERENCE PEAKS
# ============================================================

experimental_minimum = raw_df["X"].min()
experimental_maximum = raw_df["X"].max()

# Keep only reference peaks within the experimental scan range
reference_peaks_in_range = all_reference_peaks[
    (
        all_reference_peaks["Reference_2theta"]
        >= experimental_minimum
    )
    &
    (
        all_reference_peaks["Reference_2theta"]
        <= experimental_maximum
    )
].copy()

# Remove extremely weak reference peaks
reference_peaks_for_matching = reference_peaks_in_range[
    reference_peaks_in_range["Reference_Relative_Intensity"]
    >= minimum_reference_relative_intensity
].copy()


# ============================================================
# 9. COMPARE EVERY EXPERIMENTAL PEAK WITH EVERY REFERENCE PHASE
# ============================================================

all_matches = []
unmatched_peaks = []

for _, detected_row in detected_peaks.iterrows():

    detected_2theta = detected_row["Detected_2theta"]

    candidate_references = reference_peaks_for_matching.copy()

    candidate_references["Difference"] = abs(
        candidate_references["Reference_2theta"]
        - detected_2theta
    )

    possible_matches = candidate_references[
        candidate_references["Difference"] <= tolerance
    ].copy()

    if possible_matches.empty:

        # Find the nearest reference peak for information purposes
        nearest_index = candidate_references["Difference"].idxmin()

        nearest_peak = candidate_references.loc[nearest_index]

        unmatched_peaks.append({
            "Peak_ID": detected_row["Peak_ID"],
            "Detected_2theta": detected_2theta,
            "Detected_d": detected_row["Detected_d"],
            "Experimental_Intensity":
                detected_row["Experimental_Intensity"],
            "Experimental_Relative_Intensity":
                detected_row["Experimental_Relative_Intensity"],
            "Nearest_Reference_2theta":
                nearest_peak["Reference_2theta"],
            "Nearest_Phase": nearest_peak["Phase"],
            "Nearest_Reference_Card":
                nearest_peak["Reference_Card"],
            "Nearest_Difference":
                nearest_peak["Difference"],
            "Status": "NO MATCH"
        })

        continue

    # Sort by smallest 2theta difference first
    possible_matches = possible_matches.sort_values(
        by=[
            "Difference",
            "Reference_Relative_Intensity"
        ],
        ascending=[True, False]
    )

    for match_rank, (_, reference_row) in enumerate(
        possible_matches.iterrows(),
        start=1
    ):

        all_matches.append({
            "Peak_ID": detected_row["Peak_ID"],
            "Detected_2theta": detected_2theta,
            "Detected_d": detected_row["Detected_d"],
            "Experimental_Intensity":
                detected_row["Experimental_Intensity"],
            "Experimental_Relative_Intensity":
                detected_row["Experimental_Relative_Intensity"],
            "Reference_2theta":
                reference_row["Reference_2theta"],
            "Reference_d":
                reference_row["Reference_d"],
            "Reference_Intensity":
                reference_row["Reference_Intensity"],
            "Reference_Relative_Intensity":
                reference_row["Reference_Relative_Intensity"],
            "Difference":
                reference_row["Difference"],
            "HKL":
                reference_row["HKL"],
            "Phase":
                reference_row["Phase"],
            "Mineral_Name":
                reference_row["Mineral_Name"],
            "Formula":
                reference_row["Formula"],
            "Reference_Card":
                reference_row["Reference_Card"],
            "Match_Rank_For_Peak": match_rank,
            "Status": "MATCH"
        })


all_matches_df = pd.DataFrame(all_matches)
unmatched_peaks_df = pd.DataFrame(unmatched_peaks)


# ============================================================
# 10. HIERARCHICAL PHASE-FIRST MATCHING
# ============================================================

if not all_matches_df.empty:

    # 1. Define Phase Priority Hierarchy
    # Host matrix (Covellite CuS) takes precedence over sub-stoichiometric
    # or dopant phases when both match within tolerance.
    phase_priority = {
        "000752233": 1,  # CuS (Primary Host Matrix)
        "000752183": 2,  # Sn2S3 (Primary Secondary/Dopant Phase)
        "000360380": 3,  # Cu39S28 (Spionkopite sub-phase)
        "000711383": 4,  # Cu1.9375S (Djurleite sub-phase)
        "000822117": 5,  # CuCl (Precursor residue)
        "000890474": 6,  # SnCl2 (Precursor residue)
    }

    # Map priority rank (default to 99 if a card isn't explicitly listed)
    all_matches_df["Phase_Priority"] = (
        all_matches_df["Reference_Card"]
        .map(phase_priority)
        .fillna(99)
        .astype(int)
    )

    # 2. Sort by:
    #    a) Peak_ID
    #    b) Phase_Priority (Host CuS tested first)
    #    c) Difference (Smallest delta within the same phase)
    #    d) Reference_Relative_Intensity (Strongest reference peak)
    best_matches_df = (
        all_matches_df
        .sort_values(
            by=[
                "Peak_ID",
                "Phase_Priority",
                "Difference",
                "Reference_Relative_Intensity"
            ],
            ascending=[True, True, True, False]
        )
        .drop_duplicates(
            subset=["Peak_ID"],
            keep="first"
        )
        .reset_index(drop=True)
    )

else:
    best_matches_df = pd.DataFrame()

# ============================================================
# 11. CREATE PHASE SUMMARY
# ============================================================

if not all_matches_df.empty:

    # A single experimental peak may match more than one reference
    # peak from the same phase. Keep only the closest phase match
    # for each experimental peak.
    unique_phase_matches = (
        all_matches_df
        .sort_values("Difference")
        .drop_duplicates(
            subset=["Peak_ID", "Reference_Card"],
            keep="first"
        )
    )

    phase_summary = (
        unique_phase_matches
        .groupby(
            [
                "Reference_Card",
                "Phase",
                "Mineral_Name",
                "Formula"
            ],
            dropna=False
        )
        .agg(
            Matched_Experimental_Peaks=(
                "Peak_ID",
                "nunique"
            ),
            Mean_2theta_Difference=(
                "Difference",
                "mean"
            ),
            Maximum_2theta_Difference=(
                "Difference",
                "max"
            ),
            Experimental_Intensity_Support=(
                "Experimental_Relative_Intensity",
                "sum"
            )
        )
        .reset_index()
    )

    total_detected_peaks = detected_peaks["Peak_ID"].nunique()

    phase_summary["Experimental_Peak_Match_Percentage"] = (
        phase_summary["Matched_Experimental_Peaks"]
        / total_detected_peaks
        * 100
    )

    # Count reference peaks available within scan range
    reference_peak_counts = (
        reference_peaks_for_matching
        .groupby("Reference_Card")
        .size()
        .rename("Reference_Peaks_In_Range")
        .reset_index()
    )

    phase_summary = phase_summary.merge(
        reference_peak_counts,
        on="Reference_Card",
        how="left"
    )

    phase_summary["Reference_Peak_Coverage_Percentage"] = (
        phase_summary["Matched_Experimental_Peaks"]
        / phase_summary["Reference_Peaks_In_Range"]
        * 100
    )

    phase_summary = phase_summary.sort_values(
        by=[
            "Matched_Experimental_Peaks",
            "Experimental_Intensity_Support",
            "Mean_2theta_Difference"
        ],
        ascending=[False, False, True]
    ).reset_index(drop=True)

else:
    phase_summary = pd.DataFrame()


# ============================================================
# 12. CREATE FILE-STATUS TABLE
# ============================================================

reference_file_status = []

for card_number, phase_information in selected_reference_phases.items():

    file_path = find_reference_file(
        reference_folder,
        card_number
    )

    reference_file_status.append({
        "Reference_Card": card_number,
        "Phase": phase_information["Phase"],
        "Mineral_Name": phase_information["Mineral_Name"],
        "Formula": phase_information["Formula"],
        "File_Found": file_path is not None,
        "File_Name": file_path.name if file_path else None
    })

reference_file_status_df = pd.DataFrame(reference_file_status)


# ============================================================
# 13. PRINT IMPORTANT RESULTS
# ============================================================

print("\nDETECTED EXPERIMENTAL PEAKS")
print(detected_peaks.to_string(index=False))

print("\nBEST MATCH FOR EACH EXPERIMENTAL PEAK")
if not best_matches_df.empty:
    # Select key columns so you can instantly see which phase it matched with
    display_cols = [
        "Peak_ID", "Detected_2theta", "Phase", "Formula", 
        "Reference_2theta", "Difference", "HKL"
    ]
    print(best_matches_df[display_cols].to_string(index=False))
else:
    print("No matches recorded.")

if not unmatched_peaks_df.empty:
    print("\nUNMATCHED EXPERIMENTAL PEAKS")
    print(unmatched_peaks_df.to_string(index=False))

# ============================================================
# 14. SAVE RESULTS TO EXCEL
# ============================================================

output_file_path.parent.mkdir(
    parents=True,
    exist_ok=True
)

with pd.ExcelWriter(
    output_file_path,
    engine="openpyxl"
) as writer:

    detected_peaks.to_excel(
        writer,
        sheet_name="Detected_Peaks",
        index=False
    )

    all_reference_peaks.to_excel(
        writer,
        sheet_name="All_Reference_Peaks",
        index=False
    )

    reference_peaks_for_matching.to_excel(
        writer,
        sheet_name="Reference_Peaks_Used",
        index=False
    )

    all_matches_df.to_excel(
        writer,
        sheet_name="All_Possible_Matches",
        index=False
    )

    best_matches_df.to_excel(
        writer,
        sheet_name="Closest_Matches",
        index=False
    )

    unmatched_peaks_df.to_excel(
        writer,
        sheet_name="Unmatched_Peaks",
        index=False
    )

    phase_summary.to_excel(
        writer,
        sheet_name="Phase_Summary",
        index=False
    )

    reference_file_status_df.to_excel(
        writer,
        sheet_name="Reference_File_Status",
        index=False
    )


print("\nAnalysis completed successfully.")
print(f"Results saved to:\n{output_file_path}")