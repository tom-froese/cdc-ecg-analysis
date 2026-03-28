# CDC-ECG Analysis

**MATLAB toolbox for Cardiac Duty Cycle (CDC) analysis of ECG data.**
Full reproducibility package for the *Nature Aging* Brief Communication (in prep.).

[![MATLAB](https://img.shields.io/badge/MATLAB-R2024a+-blue)](https://mathworks.com)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Zenodo Data](https://img.shields.io/badge/Data-Zenodo-blue)](https://doi.org/10.5281/zenodo.19246123)

---

## How to reproduce the paper

1. Clone this repository:

   ```bash
   git clone https://github.com/tom-froese/cdc-ecg-analysis.git
   cd cdc-ecg-analysis
   ```

2. Download the preprocessed data archive from [Zenodo](https://doi.org/10.5281/zenodo.19246123) and unzip it into `data/preprocessed/`.

3. Open MATLAB and run:

   ```matlab
   setup
   reproduce_paper
   ```

All statistical results will appear in `results/` and figures in `results/figures/`.

> **Note:** Supplementary Figure 3 (pipeline validation against manual annotations) requires the raw LUDB and QTDB databases from PhysioNet. See [Pipeline validation](#pipeline-validation) below.

---

## Folder structure

```
cdc-ecg-analysis/
├── code/
│   ├── analysis/          ← statistical analyses (analyze_*.m)
│   ├── preprocessing/     ← export scripts and signal readers
│   ├── visualization/     ← figure generation (plot_*.m)
│   └── utils/             ← shared utilities
├── data/
│   └── preprocessed/      ← *_beats.csv files (download from Zenodo)
├── results/
│   └── figures/           ← generated figures (PDF + PNG)
├── config.m               ← central path configuration
├── setup.m                ← adds all project paths
├── reproduce_paper.m      ← master reproduction script
├── LICENSE
└── README.md
```

---

## Required software

- **MATLAB R2024a** or newer
- **Statistics and Machine Learning Toolbox** (for `fitlm`, `fitglme`)

No third-party ECG processing libraries are required.

---

## Pipeline validation

The automatic R-peak and T-end detectors are validated against expert manual annotations in the LUDB and QTDB databases. To reproduce this validation:

1. Download [LUDB](https://doi.org/10.13026/eegm-h675) and [QTDB](https://doi.org/10.13026/C2HS3T) from PhysioNet into `data/raw/`.
2. Run:

   ```matlab
   analyze_gold_standard_validation();
   plot_SI_Fig3();
   ```

---

## Citation

> Froese, T. et al. (in prep.). The healthy heart converges on a thermodynamic optimum (1/*e*) that predicts survival. *Nature Aging*.

> Froese, T. (2026). CDC Analysis — Preprocessed Beat Data (v1.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.19246123

---

## License

Code: [MIT](LICENSE)
Data: [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/) (Zenodo deposit)

---

## Contact

Tom Froese — Embodied Cognitive Science Unit, OIST
tom.froese@oist.jp
