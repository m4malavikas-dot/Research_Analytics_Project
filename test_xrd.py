import pandas as pd
from scipy.signal import find_peaks

file_path = r"C:\Users\user\OneDrive\Desktop\Research_Analytics_Project\01_Raw_Data\XRD\CuS B1.xlsx"

df = pd.read_excel(file_path)

x = df["X"]
y = df["Y (counts)"]

peaks, properties = find_peaks(
    y,
    distance=30,
    prominence=50
)

detected_peaks = df.iloc[peaks][["X", "d (Å)", "Y (counts)"]].copy()

reference_data = [
    {"Reference_2theta": 27.307, "HKL": "100", "Phase": "CuS"},
    {"Reference_2theta": 27.897, "HKL": "101", "Phase": "CuS"},
    {"Reference_2theta": 29.511, "HKL": "102", "Phase": "CuS"},
    {"Reference_2theta": 32.315, "HKL": "006", "Phase": "CuS"},
    {"Reference_2theta": 48.216, "HKL": "110", "Phase": "CuS"},
    {"Reference_2theta": 52.850, "HKL": "108", "Phase": "CuS"},
    {"Reference_2theta": 59.331, "HKL": "116", "Phase": "CuS"},
]

reference_df = pd.DataFrame(reference_data)

tolerance = 0.3

matched_results = []

for _, detected_row in detected_peaks.iterrows():
    detected_2theta = detected_row["X"]

    reference_df["Difference"] = abs(reference_df["Reference_2theta"] - detected_2theta)

    closest_match = reference_df.loc[reference_df["Difference"].idxmin()]

    if closest_match["Difference"] <= tolerance:
        status = "MATCH"
        phase = closest_match["Phase"]
        hkl = closest_match["HKL"]
        reference_2theta = closest_match["Reference_2theta"]
    else:
        status = "NO MATCH"
        phase = None
        hkl = None
        reference_2theta = closest_match["Reference_2theta"]

    matched_results.append({
        "Detected_2theta": detected_2theta,
        "Detected_d": detected_row["d (Å)"],
        "Intensity": detected_row["Y (counts)"],
        "Reference_2theta": reference_2theta,
        "Difference": closest_match["Difference"],
        "HKL": hkl,
        "Phase": phase,
        "Status": status
    })

matched_df = pd.DataFrame(matched_results)

print(matched_df)