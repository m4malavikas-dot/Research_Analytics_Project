# Automated Powder XRD Phase Identification & Crystallographic Analysis

An automated Python analytics pipeline for processing powder X-ray diffraction (XRD) data, designed for pristine Covellite ($\text{CuS}$) and $\text{Sn}$-doped $\text{CuS}$ semiconductor materials.

---

## Project Overview

Manual XRD peak indexing and phase matching across multi-phase and doped chalcogenide systems is prone to subjective baseline thresholding and reference card fragmentation. This pipeline automates the entire workflow from raw diffractometer scans (`.xlsx`) and standard ICDD/JCPDS reference cards (`.csv`) to a structured 8-sheet analytical report.

---

## Core Scientific & Algorithmic Features

### 1. Instrument-Independent Dynamic Peak Detection
* **Automatic Step-Size Derivation ($\Delta 2\theta$):** Computes angular resolution directly from the experimental scan array.
* **Optical Resolution Distance:** Converts physical instrumental FWHM broadening limits ($0.15^\circ\ 2\theta$) into scan-specific data point distances.
* **Savitzky-Golay Filtering:** Suppresses high-frequency detector noise using polynomial smoothing while preserving peak profile geometry.

### 2. Baseline-Aware IUPAC $5\sigma$ Limit of Quantitation (LOQ)
* **Substrate De-trending:** Applies a rolling baseline filter ($2.5^\circ$ window) to strip low-angle amorphous glass substrate scattering ($<25^\circ\ 2\theta$).
* **Statistical Noise Floor ($\sigma_{\text{noise}}$):** Calculates real-time background noise standard deviation from high-frequency residuals ($y - y_{\text{smooth}}$).
* **Universal Cutoff:** Enforces an IUPAC $5\sigma$ LOQ prominence floor, reliably capturing weak dopant/secondary reflections (e.g., $\text{Sn}_2\text{S}_3$ at $\sim 21.4^\circ$ and $\text{CuS}\ (101)$ at $\sim 27.2^\circ$) while rejecting background fluctuations.

### 3. Chemistry-Filtered & Hierarchical Phase Matching
* **Sample-Aware Reference Selection:** Automatically detects sample chemistry (`CuS` vs. `Sn-CuS`) and filters candidate reference cards accordingly.
* **Host-Matrix Priority:** Ranks candidate matches within a $\pm 0.3^\circ\ 2\theta$ tolerance window using a hierarchical phase-first model:
  1. Primary Host Matrix: **$\text{CuS}$ (Covellite)**
  2. Secondary Dopant Phase: **$\text{Sn}_2\text{S}_3$**
  3. Sub-stoichiometric Copper Sulfides: **Cu<sub>39</sub>S<sub>28</sub> (Spionkopite)**, **$\text{Cu}_{1.9375}\text{S}$ (Djurleite)**
  4. Precursor Residuals: **$\text{CuCl}$**, **$\text{SnCl}_2$**

---

## Automated Excel Output (`XRD_Results/`)

Each run generates a comprehensive `.xlsx` workbook containing 8 structured sheets:
1. **`Detected_Peaks`**: Experimental $2\theta$, $d$-spacing, raw counts, prominence, and normalised relative intensities.
2. **`All_Reference_Peaks`**: Complete theoretical library parsed from reference `.csv` cards.
3. **`Reference_Peaks_Used`**: Reference reflections filtered to the experimental scan range and intensity cutoff ($\ge 5\%$).
4. **`All_Possible_Matches`**: Every candidate phase reflection within $\pm 0.3^\circ\ 2\theta$.
5. **`Closest_Matches`**: Priority-resolved primary assignment for each experimental reflection.
6. **`Unmatched_Peaks`**: Unindexed reflections with angular distances to the nearest theoretical peak for strain/defect analysis.
7. **`Phase_Summary`**: Phase coverage percentages, matched peak counts, intensity support, and mean $\Delta 2\theta$ shifts.
8. **`Reference_File_Status`**: Verification table of loaded reference cards.

---

## Repository Structure

* `01_Raw_Data/XRD/` — Experimental diffractometer scans (`CuS A1`, `CuS A2`, `CuS B1`, `CuS B2`, `Sn-CuS`).
* `03_Reference_Data/` — Standard ICDD/JCPDS reference diffraction cards (`.csv`).
* `06_Python/` — Core analysis script (`multi_reference_xrd(1).py`) and output reports (`XRD_Results/`).
