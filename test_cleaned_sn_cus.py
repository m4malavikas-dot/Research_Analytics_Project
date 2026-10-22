import re
from pathlib import Path

import pandas as pd
import numpy as np
from scipy.signal import find_peaks, savgol_filter, peak_prominences


# ============================================================
# 1. FOLDER LOCATIONS
# ============================================================

cleaned_data_folder = Path(
    r"C:\Users\user\OneDrive\Desktop\Research_Analytics_Project"
    r"\02_Cleaned_Data"
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
# 2. SETTINGS
# ============================================================

expected_sn2s3_angle = 21.466

# Reference matching tolerance
tolerance = 0.3

# Main peak detector
main_prominence = 0.01

# Sensitive detector used only near 21.466 degrees
sn_region_prominence = 0.001

# Region examined around the expected Sn2S3 peak
sn_region_minimum = 20.0
sn_region_maximum = 23.0

# Ignore reference lines below this relative intensity
minimum_reference_relative_intensity = 0.0


# ============================================================
# 3. REFERENCE PHASES
# ============================================================

reference_phases = {
    "000752233": {
        "Phase": "CuS",
        "Mineral_Name": "Covellite",
        "Formula": "CuS"
    },

    "000360380": {
        "Phase": "Cu39S28",
        "Mineral_Name": "Spionkopite",
        "Formula": "Cu39S28"
    },

    "000711383": {
        "Phase": "Cu1.9375S",
        "Mineral_Name": "Djurleite",
        "Formula": "Cu1.9375S"
    },

    "000822117": {
        "Phase": "CuCl",
        "Mineral_Name": None,
        "Formula": "CuCl"
    },

    "000752183": {
        "Phase": "Sn2S3",
        "Mineral_Name": None,
        "Formula": "Sn2S3"
    },

    "000890474": {
        "Phase": "SnCl2",
        "Mineral_Name": None,
        "Formula": "SnCl2"
    }
}


# ============================================================
# 4. HELPER FUNCTIONS
# ============================================================

def choose_cleaned_file(folder):
    """
    Shows cleaned CSV and Excel files and asks the user
    which one should be analysed.
    """

    cleaned_files = sorted(
        list(folder.glob("*.xlsx"))
        + list(folder.glob("*.xls"))
        + list(folder.glob("*.csv"))
    )

    if not cleaned_files:
        raise FileNotFoundError(
            f"No CSV or Excel files were found in:\n{folder}"
        )

    print("\nAVAILABLE CLEANED XRD FILES")

    for number, file_path in enumerate(cleaned_files, start=1):
        print(f"{number}. {file_path.name}")

    while True:
        try:
            choice = int(
                input(
                    "\nEnter the number of the cleaned file "
                    "you want to analyse: "
                )
            )

            if 1 <= choice <= len(cleaned_files):
                selected_file = cleaned_files[choice - 1]

                print(
                    f"\nSelected cleaned file: "
                    f"{selected_file.name}"
                )

                return selected_file

            print(
                f"Please enter a number between "
                f"1 and {len(cleaned_files)}."
            )

        except ValueError:
            print("Please enter a valid number.")


def read_cleaned_data(file_path):
    """
    Reads the cleaned XRD file.

    If the file is Excel, all sheet names are displayed
    so the user can select the sheet containing 2theta
    and intensity.
    """

    print(f"\nPython is opening this file:\n{file_path}")

    if file_path.suffix.lower() == ".csv":
        return pd.read_csv(file_path)

    excel_file = pd.ExcelFile(file_path)

    print("\nAVAILABLE SHEETS")

    for number, sheet_name in enumerate(
        excel_file.sheet_names,
        start=1
    ):
        print(f"{number}. {sheet_name}")

    while True:
        try:
            choice = int(
                input(
                    "\nEnter the number of the sheet "
                    "containing 2theta and intensity: "
                )
            )

            if 1 <= choice <= len(excel_file.sheet_names):

                selected_sheet = excel_file.sheet_names[
                    choice - 1
                ]

                print(
                    f"\nSelected sheet: {selected_sheet}"
                )

                dataframe = pd.read_excel(
                    file_path,
                    sheet_name=selected_sheet
                )

                print(
                    f"\nColumns found in selected sheet: "
                    f"{list(dataframe.columns)}"
                )

                return dataframe

            print(
                f"Please enter a number between 1 and "
                f"{len(excel_file.sheet_names)}."
            )

        except ValueError:
            print("Please enter a valid number.")

def identify_cleaned_columns(dataframe):
    """
    Finds the 2theta and intensity columns even if the exact
    spelling or capitalisation is slightly different.
    """

    theta_column = None
    intensity_column = None

    for column in dataframe.columns:

        cleaned_name = (
            str(column)
            .strip()
            .lower()
            .replace(" ", "")
            .replace("_", "")
            .replace("-", "")
            .replace("(", "")
            .replace(")", "")
        )

        if (
            cleaned_name in ["2theta", "twotheta", "x"]
            or "2theta" in cleaned_name
        ):
            theta_column = column

        if (
            cleaned_name in [
                "intensity",
                "y",
                "counts",
                "ycounts"
            ]
            or "intensity" in cleaned_name
            or "count" in cleaned_name
        ):
            intensity_column = column

    if theta_column is None or intensity_column is None:
        raise ValueError(
            "\nCould not identify the required columns.\n"
            f"Columns found: {list(dataframe.columns)}\n"
            "The file must contain one 2theta column and "
            "one intensity column."
        )

    return theta_column, intensity_column


def normalise_card_number(value):
    """
    Removes spaces and hyphens from a reference-card number.
    """

    return re.sub(r"[^0-9]", "", str(value))


def find_reference_file(folder, card_number):
    """
    Finds a reference CSV file by its reference-card number.
    """

    required_number = normalise_card_number(card_number)

    for file_path in folder.glob("*.csv"):

        file_number = normalise_card_number(
            file_path.stem
        )

        if file_number == required_number:
            return file_path

    return None


def read_text_file(file_path):
    """
    Reads the reference-card file using several possible encodings.
    """

    encodings = [
        "utf-8-sig",
        "utf-8",
        "cp1252",
        "latin1"
    ]

    for encoding in encodings:
        try:
            with open(
                file_path,
                "r",
                encoding=encoding
            ) as file:
                return file.readlines()

        except UnicodeDecodeError:
            continue

    raise ValueError(
        f"Could not read reference file:\n{file_path}"
    )


def extract_reference_peaks(
    file_path,
    card_number,
    phase_information
):
    """
    Extracts diffraction peaks from the reference-card CSV.
    """

    lines = read_text_file(file_path)

    diffraction_start = None

    for index, line in enumerate(lines):

        cleaned_line = (
            line.lower()
            .replace('"', "")
        )

        if (
            "2theta" in cleaned_line
            and "int" in cleaned_line
            and "hkl" in cleaned_line
        ):
            diffraction_start = index + 1
            break

    if diffraction_start is None:
        raise ValueError(
            f"Could not locate diffraction data in:\n"
            f"{file_path}"
        )

    diffraction_pattern = re.compile(
        r"^\s*"
        r"([+-]?\d+(?:\.\d+)?)"
        r"[\s,\t;]+"
        r"([+-]?\d+(?:\.\d+)?)"
        r"[\s,\t;]+"
        r"([+-]?\d+(?:\.\d+)?)"
        r"[\s,\t;]+"
        r"(.+?)"
        r"\s*$"
    )

    extracted_rows = []

    for line in lines[diffraction_start:]:

        cleaned_line = (
            line.strip()
            .replace('"', "")
        )

        if not cleaned_line:
            continue

        match = diffraction_pattern.match(
            cleaned_line
        )

        if match is None:
            continue

        hkl = match.group(4).strip()

        hkl = re.sub(
            r"[\t,;]+",
            " ",
            hkl
        )

        hkl = re.sub(
            r"\s+",
            " ",
            hkl
        )

        extracted_rows.append({
            "Reference_Card": card_number,
            "Phase": phase_information["Phase"],
            "Mineral_Name":
                phase_information["Mineral_Name"],
            "Formula": phase_information["Formula"],
            "Reference_2theta":
                float(match.group(1)),
            "Reference_d":
                float(match.group(2)),
            "Reference_Intensity":
                float(match.group(3)),
            "HKL": hkl,
            "Reference_File": file_path.name
        })

    if not extracted_rows:
        raise ValueError(
            f"No diffraction peaks were extracted from:\n"
            f"{file_path}"
        )

    reference_df = pd.DataFrame(
        extracted_rows
    )

    maximum_intensity = (
        reference_df["Reference_Intensity"].max()
    )

    if maximum_intensity > 0:

        reference_df[
            "Reference_Relative_Intensity"
        ] = (
            reference_df["Reference_Intensity"]
            / maximum_intensity
            * 100
        )

    else:

        reference_df[
            "Reference_Relative_Intensity"
        ] = 0

    return reference_df


# ============================================================
# 5. SELECT AND READ CLEANED DATA
# ============================================================

cleaned_file_path = choose_cleaned_file(
    cleaned_data_folder
)

cleaned_df = read_cleaned_data(
    cleaned_file_path
)

theta_column, intensity_column = identify_cleaned_columns(
    cleaned_df
)

print(
    f"\n2theta column identified as: "
    f"{theta_column}"
)

print(
    f"Intensity column identified as: "
    f"{intensity_column}"
)

cleaned_df = cleaned_df[
    [theta_column, intensity_column]
].copy()

cleaned_df = cleaned_df.rename(columns={
    theta_column: "2theta",
    intensity_column: "Intensity"
})

cleaned_df["2theta"] = pd.to_numeric(
    cleaned_df["2theta"],
    errors="coerce"
)

cleaned_df["Intensity"] = pd.to_numeric(
    cleaned_df["Intensity"],
    errors="coerce"
)

cleaned_df = cleaned_df.dropna(
    subset=["2theta", "Intensity"]
).copy()

cleaned_df = cleaned_df.sort_values(
    "2theta"
).reset_index(drop=True)

x = cleaned_df["2theta"].to_numpy()
y = cleaned_df["Intensity"].to_numpy()


# ============================================================
# 6. CHECK THE EXPECTED Sn2S3 REGION DIRECTLY
# ============================================================

sn_region_mask = (
    (x >= sn_region_minimum)
    & (x <= sn_region_maximum)
)

sn_region_df = cleaned_df.loc[
    sn_region_mask
].copy()

print(
    f"\nCHECKING CLEANED DATA BETWEEN "
    f"{sn_region_minimum}° AND "
    f"{sn_region_maximum}°"
)

if sn_region_df.empty:

    print(
        "No cleaned-data points exist in this region."
    )

else:

    highest_sn_point = sn_region_df.loc[
        sn_region_df["Intensity"].idxmax()
    ]

    difference_from_expected = abs(
        highest_sn_point["2theta"]
        - expected_sn2s3_angle
    )

    print("\nHighest cleaned-data point in this region:")

    print(
        f"2theta = "
        f"{highest_sn_point['2theta']}"
    )

    print(
        f"Intensity = "
        f"{highest_sn_point['Intensity']}"
    )

    print(
        f"Difference from expected "
        f"{expected_sn2s3_angle}° = "
        f"{difference_from_expected:.4f}°"
    )


# ============================================================
# 7. SMOOTH CLEANED DATA
# ============================================================

number_of_points = len(y)

if number_of_points < 7:
    raise ValueError(
        "The cleaned-data file contains too few points."
    )

window_length = 11

if window_length >= number_of_points:
    window_length = number_of_points - 1

if window_length % 2 == 0:
    window_length -= 1

smoothed_y = savgol_filter(
    y,
    window_length=window_length,
    polyorder=3
)

cleaned_df["Smoothed_Intensity"] = smoothed_y

# Calculate a suitable prominence automatically
# based on the intensity scale of the cleaned data.
lower_signal_level = np.percentile(
    smoothed_y,
    5
)

upper_signal_level = np.percentile(
    smoothed_y,
    95
)

signal_range = (
    upper_signal_level
    - lower_signal_level
)

adaptive_main_prominence = max(
    0.02,
    signal_range * 0.20
)

print(
    f"\nAutomatically selected main prominence: "
    f"{adaptive_main_prominence:.6f}"
)

# ============================================================
# 8. DETECT MAIN PEAKS
# ============================================================

main_peak_indices, main_properties = find_peaks(
    smoothed_y,
    distance=30,
    prominence=adaptive_main_prominence,
    width=3
)

main_detected_peaks = cleaned_df.iloc[
    main_peak_indices
][
    [
        "2theta",
        "Intensity",
        "Smoothed_Intensity"
    ]
].copy()

main_detected_peaks[
    "Prominence"
] = main_properties["prominences"]

main_detected_peaks[
    "Detection_Type"
] = "MAIN PEAK"


# ============================================================
# 9. SENSITIVE SEARCH NEAR 21.466 DEGREES
# ============================================================

sn_region_positions = np.where(
    sn_region_mask
)[0]

if len(sn_region_positions) > 0:

    sn_region_smoothed_y = smoothed_y[
        sn_region_positions
    ]

    local_sn_indices, local_sn_properties = find_peaks(
        sn_region_smoothed_y,
        distance=5,
        prominence=sn_region_prominence,
        width=1
    )

    global_sn_indices = sn_region_positions[
        local_sn_indices
    ]

    sensitive_sn_peaks = cleaned_df.iloc[
        global_sn_indices
    ][
        [
            "2theta",
            "Intensity",
            "Smoothed_Intensity"
        ]
    ].copy()

    sensitive_sn_peaks[
        "Prominence"
    ] = local_sn_properties["prominences"]

    sensitive_sn_peaks[
        "Difference_From_21.466"
    ] = abs(
        sensitive_sn_peaks["2theta"]
        - expected_sn2s3_angle
    )

    sensitive_sn_peaks[
        "Detection_Type"
    ] = "WEAK Sn REGION PEAK"

    sensitive_sn_peaks = sensitive_sn_peaks.sort_values(
        "Difference_From_21.466"
    )

else:

    sensitive_sn_peaks = pd.DataFrame()
# If find_peaks did not detect a local peak, use the highest
# cleaned-data point in the expected Sn2S3 region as a
# tentative feature for diagnostic comparison.
if sensitive_sn_peaks.empty and not sn_region_df.empty:

    highest_sn_point = sn_region_df.loc[
        sn_region_df["Intensity"].idxmax()
    ]

    tentative_difference = abs(
        highest_sn_point["2theta"]
        - expected_sn2s3_angle
    )

    if tentative_difference <= 0.3:

        sensitive_sn_peaks = pd.DataFrame([{
            "2theta": highest_sn_point["2theta"],
            "Intensity": highest_sn_point["Intensity"],
            "Smoothed_Intensity": cleaned_df.loc[
                highest_sn_point.name,
                "Smoothed_Intensity"
            ],
            "Prominence": None,
            "Difference_From_21.466":
                tentative_difference,
            "Detection_Type":
                "TENTATIVE Sn2S3 FEATURE"
        }])

        print(
            "\nNo local maximum was detected by find_peaks, "
            "but the highest point near 21.466° was retained "
            "as a tentative feature."
        )


print("\nMAIN PEAKS DETECTED FROM CLEANED DATA")

if main_detected_peaks.empty:
    print("No main peaks were detected.")
else:
    print(
        main_detected_peaks.to_string(
            index=False
        )
    )

print(
    "\nSENSITIVE PEAK SEARCH BETWEEN "
    "20° AND 23°"
)

if sensitive_sn_peaks.empty:

    print(
        "No local peak was detected in the Sn2S3 region, "
        f"even with prominence = "
        f"{sn_region_prominence}."
    )

else:

    print(sensitive_sn_peaks)


# ============================================================
# 10. COMBINE DETECTED PEAKS
# ============================================================

all_detected_peaks = pd.concat(
    [
        main_detected_peaks,
        sensitive_sn_peaks
    ],
    ignore_index=True,
    sort=False
)

all_detected_peaks = (
    all_detected_peaks
    .sort_values("2theta")
    .drop_duplicates(
        subset=["2theta"],
        keep="first"
    )
    .reset_index(drop=True)
)

all_detected_peaks.insert(
    0,
    "Peak_ID",
    range(1, len(all_detected_peaks) + 1)
)


# ============================================================
# 11. READ ALL REFERENCE CARDS
# ============================================================

reference_tables = []

for card_number, phase_information in reference_phases.items():

    reference_file = find_reference_file(
        reference_folder,
        card_number
    )

    if reference_file is None:

        print(
            f"WARNING: Reference file not found for "
            f"{card_number}"
        )

        continue

    print(
        f"Reading {phase_information['Phase']} "
        f"from {reference_file.name}"
    )

    card_df = extract_reference_peaks(
        reference_file,
        card_number,
        phase_information
    )

    reference_tables.append(card_df)

if not reference_tables:
    raise ValueError(
        "No reference cards could be loaded."
    )

all_reference_peaks = pd.concat(
    reference_tables,
    ignore_index=True
)

all_reference_peaks = all_reference_peaks[
    all_reference_peaks[
        "Reference_Relative_Intensity"
    ] >= minimum_reference_relative_intensity
].copy()


# ============================================================
# 12. COMPARE CLEANED-DATA PEAKS WITH REFERENCES
# ============================================================

comparison_results = []

for _, detected_row in all_detected_peaks.iterrows():

    detected_angle = detected_row["2theta"]

    candidate_references = all_reference_peaks.copy()

    candidate_references["Difference"] = abs(
        candidate_references["Reference_2theta"]
        - detected_angle
    )

    possible_matches = candidate_references[
        candidate_references["Difference"]
        <= tolerance
    ].copy()

    if possible_matches.empty:

        comparison_results.append({
            "Peak_ID": detected_row["Peak_ID"],
            "Detected_2theta": detected_angle,
            "Detected_Intensity":
                detected_row["Intensity"],
            "Detected_Prominence":
                detected_row["Prominence"],
            "Detection_Type":
                detected_row["Detection_Type"],
            "Reference_2theta": None,
            "Difference": None,
            "Phase": None,
            "Formula": None,
            "HKL": None,
            "Reference_Card": None,
            "Status": "NO MATCH"
        })

        continue

    possible_matches = possible_matches.sort_values(
        [
            "Difference",
            "Reference_Relative_Intensity"
        ],
        ascending=[True, False]
    )

    for _, reference_row in possible_matches.iterrows():

        comparison_results.append({
            "Peak_ID": detected_row["Peak_ID"],
            "Detected_2theta": detected_angle,
            "Detected_Intensity":
                detected_row["Intensity"],
            "Detected_Prominence":
                detected_row["Prominence"],
            "Detection_Type":
                detected_row["Detection_Type"],
            "Reference_2theta":
                reference_row["Reference_2theta"],
            "Difference":
                reference_row["Difference"],
            "Phase":
                reference_row["Phase"],
            "Formula":
                reference_row["Formula"],
            "HKL":
                reference_row["HKL"],
            "Reference_Card":
                reference_row["Reference_Card"],
            "Reference_Relative_Intensity":
                reference_row[
                    "Reference_Relative_Intensity"
                ],
            "Status": "MATCH"
        })


comparison_columns = [
    "Peak_ID",
    "Detected_2theta",
    "Detected_Intensity",
    "Detected_Prominence",
    "Detection_Type",
    "Reference_2theta",
    "Difference",
    "Phase",
    "Formula",
    "HKL",
    "Reference_Card",
    "Reference_Relative_Intensity",
    "Status"
]

comparison_df = pd.DataFrame(
    comparison_results,
    columns=comparison_columns
)

# ============================================================
# 13. SPECIFIC CHECK FOR THE 21.466 DEGREE FEATURE
# ============================================================

expected_region_matches = comparison_df[
    (
        comparison_df["Detected_2theta"]
        >= expected_sn2s3_angle - 0.5
    )
    &
    (
        comparison_df["Detected_2theta"]
        <= expected_sn2s3_angle + 0.5
    )
].copy()

print("\nRESULT NEAR THE EXPECTED 21.466° FEATURE")

if expected_region_matches.empty:

    print(
        "No experimental local peak was detected within "
        "±0.5° of 21.466°."
    )

else:

    print(expected_region_matches)


# ============================================================
# 14. SAVE THE DIAGNOSTIC RESULT
# ============================================================

output_folder.mkdir(
    parents=True,
    exist_ok=True
)

output_file_path = (
    output_folder
    / f"{cleaned_file_path.stem}_Cleaned_Data_Test.xlsx"
)

with pd.ExcelWriter(
    output_file_path,
    engine="openpyxl"
) as writer:

    cleaned_df.to_excel(
        writer,
        sheet_name="Cleaned_Data",
        index=False
    )

    sn_region_df.to_excel(
        writer,
        sheet_name="Sn_Region_Raw_Points",
        index=False
    )

    main_detected_peaks.to_excel(
        writer,
        sheet_name="Main_Detected_Peaks",
        index=False
    )

    sensitive_sn_peaks.to_excel(
        writer,
        sheet_name="Sensitive_Sn_Region",
        index=False
    )

    all_detected_peaks.to_excel(
        writer,
        sheet_name="All_Detected_Peaks",
        index=False
    )

    comparison_df.to_excel(
        writer,
        sheet_name="Reference_Comparison",
        index=False
    )

    expected_region_matches.to_excel(
        writer,
        sheet_name="Expected_21_466_Check",
        index=False
    )

    all_reference_peaks.to_excel(
        writer,
        sheet_name="Reference_Peaks",
        index=False
    )


print("\nCleaned-data diagnostic completed.")

print(
    f"Results saved to:\n"
    f"{output_file_path}"
)