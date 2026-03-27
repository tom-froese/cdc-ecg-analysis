# CDC-ECG Analysis

**MATLAB toolbox for Cardiac Duty Cycle (CDC) analysis of ECG data.**  
Full reproducibility package for the *Nature Aging* Brief Communication (2026).

[![MATLAB](https://img.shields.io/badge/MATLAB-R2024a+-blue)](https://mathworks.com)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Zenodo Data](https://img.shields.io/badge/Data-Zenodo-blue)](https://zenodo.org) 

**One-click reproduction of all figures and tables.**

---

## How to reproduce the paper in < 5 minutes

1. Clone this repository:
   ```bash
   git clone https://github.com/tom-froese/cdc-ecg-analysis.git
   cd cdc-ecg-analysis

2. Download the preprocessed data archive from Zenodo and unzip it into data/preprocessed/.

3. Open MATLAB and run:
setup
reproduce_paper

All statistical results will appear in results/ and figures will appear in results/figures/ and tables in results/tables/.

Folder structure

cdc-ecg-analysis/
├── code/                  ← All MATLAB functions
│   ├── preprocessing/     ← export_* scripts and readers
│   ├── analysis/          ← quality filters, statistical models
│   ├── visualization/     ← figure generation
│   └── utils/             ← shared utilities
├── data/
│   └── preprocessed/      ← *_beats.csv files (required for figures)
├── results/
│   ├── figures/           ← generated figures (Nature style)
│   └── tables/            ← generated tables
├── setup.m
├── config.m
├── reproduce_paper.m      ← master reproduction script
└── README.md

Required software

MATLAB R2024a (or newer)

Citation
Tom Froese, et al. (2026). The healthy heart converges on a thermodynamic optimum (1/e) that predicts survival. Nature Aging (in prep.).
Code DOI: https://github.com/tom-froese/cdc-ecg-analysis
Data DOI: https://doi.org/10.5281/zenodo.19246122

Questions or issues?
Open an issue on this repository or contact:
Tom Froese (tom.froese@oist.jp | @DrTomFroese)