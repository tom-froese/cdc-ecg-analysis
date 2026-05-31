# CDC-ECG Analysis

**MATLAB toolbox for Cardiac Duty Cycle (CDC) analysis of ECG data.**
Full reproducibility package for the *Nature Aging* Brief Communication (submitted; under review).

Authors: Tom Froese, Vaibhav Bhaskar, Ruben Fossion. An independent Python reimplementation by Vaibhav Bhaskar is available at [VAIBHAV-BHASKAR/cardiac-duty-cycle](https://github.com/VAIBHAV-BHASKAR/cardiac-duty-cycle).

[![MATLAB](https://img.shields.io/badge/MATLAB-R2024a+-blue)](https://mathworks.com)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Zenodo Data](https://img.shields.io/badge/Data-Zenodo-blue)](https://doi.org/10.5281/zenodo.19246122)

---

## How to reproduce the paper

1. Clone this repository:

   ```bash
   git clone https://github.com/tom-froese/cdc-ecg-analysis.git
   cd cdc-ecg-analysis
   ```

2. Download the preprocessed data archive from [Zenodo](https://doi.org/10.5281/zenodo.19246122) and unzip it into `data/preprocessed/`.

3. Open MATLAB and run:

   ```matlab
   setup
   reproduce_paper
   ```

All statistical results will appear in `results/` and figures in `results/figures/`.

> **Note:** Supplementary Figure 3 (pipeline validation against manual annotations) requires the raw LUDB and QTDB databases from PhysioNet. See [Pipeline validation](#pipeline-validation) below.

> **Figure outputs:** Each figure is exported as both PDF (vector or raster-in-PDF for large scatter plots) and PNG. A MATLAB-editable `.fig` is also saved for figures with a manageable number of graphic objects (Fig. 2, SI Fig. 1, SI Fig. 2, SI Fig. 4, SI Fig. 5). Figures with very large scatter plots (Fig. 1, SI Figs. 3, 6, 7, 8) skip the `.fig` export, since serialising every point can consume gigabytes of memory; their PDFs are produced via `exportgraphics` with `ContentType='image'` (300 dpi raster) instead.

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

1. Download [LUDB](https://doi.org/10.13026/eegm-h675) and [QTDB](https://doi.org/10.13026/C24K53) from PhysioNet into `data/raw/`.
2. Run:

   ```matlab
   analyze_gold_standard_validation();   % code/analysis/
   plot_SI_Fig3();                       % code/visualization/
   ```

---

## Citation

> Froese, T., Bhaskar, V., & Fossion, R. (2026). The cardiac duty cycle maintains an optimal ratio (1/*e*) in healthy aging. *Nature Aging* (submitted; under review).

> Froese, T. (2026). CDC Analysis — Preprocessed Beat Data (v1.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.19246122

---

## License

Code: [MIT](LICENSE)
Data: [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/) (Zenodo deposit)

---

## Contact

Tom Froese — Embodied Cognitive Science Unit, OIST
tom.froese@oist.jp
