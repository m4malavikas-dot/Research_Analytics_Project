# XRD Multi-Reference Analysis: Version 1 vs. Version 2 (`multi_reference_xrd`)

## Overview
This document explains the architectural differences between `multi_reference_xrd(1).py` and `multi_reference_xrd(2).py` in this repository. Both scripts batch-process experimental raw XRD `.xlsx`/`.xls` files using IUPAC 5-sigma baseline-aware peak detection and match experimental peaks against reference crystallographic cards (`.csv` and `.cif`). 

The primary difference is **how crystal phase metadata (`Phase`, `Formula`, `Mineral_Name`, `Contains_Sn`, and `Phase_Priority`) is identified**.

---

## Side-by-Side Comparison

| Feature | `multi_reference_xrd(1).py` (Version 1) | `multi_reference_xrd(2).py` (Version 2) |
| :--- | :--- | :--- |
| **Phase Identification** | Uses a hardcoded Python dictionary (`reference_phases_registry`) inside the script. | Automatically extracts `Phase`, `Formula`, and `Mineral_Name` directly from the input `.csv` or `.cif` reference files. |
| **Adding New Reference Cards** | Requires manually editing the Python script to add the card ID, formula, mineral name, `Contains_Sn`, and `Priority`. | **Zero code edits required.** Drop any new `.csv` or `.cif` file into `03_Reference_Data` and run the script. |
| **Tin (`Sn`) Filtering (`Contains_Sn`)** | Manually defined as `True` or `False` per card in the dictionary. | Automatically detected from the extracted chemical formula using regex (`Sn(?![a-z])`). |
| **Phase Priority (`Phase_Priority`)** | Manually assigned (`1` to `6`) in the dictionary. | Dynamically calculated from the chemical formula stoichiometry and composition. |
| **CSV Column Flexibility** | Requires 4 columns (`2theta`, `d`, `Intensity`, `hkl`). | Supports both 4-column and 3-column (`hkl` optional) reference peak tables. |
| **CIF Formula Fallback** | Falls back to the filename if `_chemical_formula_sum` is missing. | Computes the reduced chemical formula directly from the unit-cell atomic structure via `pymatgen` (`structure.composition.reduced_formula`) if CIF header tags are missing. |

---

## Detailed Technical Changes

### 1. Removal of Hardcoded `reference_phases_registry`
* **In `multi_reference_xrd(1).py`:**
  Phase metadata was hardcoded at the top of the script for specific ICDD card numbers (`000752233` for `CuS` Covellite, `000752183` for `Sn2S3` Ottemannite, `000360380` for `Cu39S28` Spionkopite, `000711383` for `Cu1.9375S` Djurleite, `000822117` for `CuCl` Nantokite, and `000890474` for `SnCl2`). Any card not listed in that dictionary relied on a basic underscore-split of the filename, which often resulted in missing mineral names or inaccurate priority rankings.

* **In `multi_reference_xrd(2).py`:**
  The hardcoded registry is removed. A new parser (`extract_metadata_from_reference()`) scans all lines above the numerical diffraction table in `.csv` files (and crystallographic tags in `.cif` files) to dynamically extract:
  * **Reference Card ID:** Parsed from `PDF`, `Card`, `COD`, `ICSD`, `Entry` header fields, or normalized digits from the filename.
  * **Chemical Formula:** Parsed from `Formula`, `Chem`, or `Empirical` header fields, free-text header comments, or filename tokens.
  * **Mineral / Compound Name:** Parsed from `Mineral`, `Compound`, `Phase`, or `Name` header fields, or filename tokens.

### 2. Regex Chemical Formula vs. Mineral Name Classifier (`is_chemical_formula`)
To prevent mineral names and chemical formulas from getting swapped when reading unstructured headers or filenames:
* Version 2 validates candidate strings against `^(?:[A-Z][a-z]?\d*(?:\.\d+)?)+$` and rejects normal English words containing two or more consecutive lowercase letters (`[a-z]{2,}`).
* Tokens like `CuS`, `Sn2S3`, `Cu39S28`, `Cu1.9375S`, `CuCl`, and `SnCl2` are classified as **Chemical Formulas**.
* Tokens like `Covellite`, `Ottemannite`, `Spionkopite`, `Djurleite`, and `Nantokite` are classified as **Mineral Names**.

### 3. Automated Chemistry Filter (`check_contains_sn`)
* During batch processing, samples whose filenames start with `CuS` (and do not contain `Sn`) are run in **Pure CuS Mode**, which excludes Tin-containing reference phases to prevent false-positive peak assignments.
* Instead of relying on a manual `Contains_Sn` boolean in the script, Version 2 evaluates `re.search(r"Sn(?![a-z])", combined)` on the extracted formula and mineral name. This detects `Sn` while avoiding false triggers on elements starting with `S` (such as `Sm`, `Se`, or `Si`).

### 4. Dynamic Phase Priority Ranking (`assign_dynamic_priority`)
When multiple reference phases fall within the `2theta` matching tolerance (`0.3°`), candidate matches are sorted by:
1. `Phase_Priority` (ascending)
2. `Difference` in `2theta` (ascending)
3. `Reference_Relative_Intensity` (descending)

Version 2 assigns `Phase_Priority` automatically from the parsed chemical formula:
* **Priority 1:** Primary host phase (`CuS`)
* **Priority 2:** Stoichiometric binary sulfides (e.g., `Sn2S3`, `SnS2`, `SnS`, `Cu2S`)
* **Priority 3:** Non-stoichiometric or complex supercell sulfides (formulas containing decimals or stoichiometric coefficients > 10, such as `Cu1.9375S` or `Cu39S28`)
* **Priority 4:** Halide salts, oxides, or precursor impurities (formulas containing `Cl` or `O`, such as `CuCl` or `SnCl2`)

---

## Reference File Formatting Requirements for Version 2

For `multi_reference_xrd(2).py` to identify phases automatically from `03_Reference_Data`, reference files should follow **any one** of these three formats:

1. **CSV with Metadata Header Lines (Recommended):**
   Include key-value metadata lines above the `2theta, d, Intensity, HKL` table:
   ```csv
   Card: 00-075-2233
   Formula: CuS
   Mineral: Covellite
   2theta,d,Intensity,hkl
   27.68,3.220,28.0,1 0 1
   29.28,3.048,67.0,1 0 2
   31.78,2.813,100.0,1 0 3