import pandas as pd
from scipy.signal import savgol_filter, find_peaks
from pathlib import Path

# --------------------------------------------------
# FOLDER CONTAINING YOUR RAW XRD FILES
# --------------------------------------------------

raw_folder = Path(
    r"C:\Users\user\OneDrive\Desktop\Research_Analytics_Project"
    r"\01_Raw_Data\XRD"
)

samples = [
    "CuS A1.xlsx",
    "CuS A2.xlsx",
    "CuS B1.xlsx",
    "CuS B2.xlsx",
    "Sn-CuS.xlsx"
]

prominences = [100, 75, 50]


# --------------------------------------------------
# TEST EACH SAMPLE
# --------------------------------------------------

for sample in samples:

    print("\n\n" + "=" * 60)
    print(f"SAMPLE: {sample}")
    print("=" * 60)

    file_path = raw_folder / sample

    df = pd.read_excel(file_path)

    x = df["X"]
    y = df["Y (counts)"]

    # Same smoothing as original code
    y_smooth = savgol_filter(
        y,
        window_length=11,
        polyorder=3
    )

    # Test different prominence values
    for prominence in prominences:

        peaks, properties = find_peaks(
            y_smooth,
            distance=30,
            prominence=prominence,
            width=3
        )

        # Only show peaks between 20° and 70°
        detected = [
            round(float(x.iloc[p]), 3)
            for p in peaks
            if 20 <= x.iloc[p] <= 70
        ]

        print(f"\nProminence = {prominence}")
        print(f"Number of peaks = {len(detected)}")
        print(detected)
