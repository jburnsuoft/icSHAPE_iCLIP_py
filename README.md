# icSHAPE_iCLIP

This repository contains tools to perform the overlap of experimental RNA-Binding Protein binding coordinates (iCLIP/eCLIP) with in vivo and in vitro icSHAPE RNA structural profiles.

Program and data analysis was incorporated into the following publication:
C2H2-zinc-finger transcription factors bind RNA and function in diverse post-transcriptional regulatory processes,
Molecular Cell,
Volume 84, Issue 19,
2024,
Pages 3810-3825.e10,
ISSN 1097-2765,
https://doi.org/10.1016/j.molcel.2024.08.037.
(https://www.sciencedirect.com/science/article/pii/S1097276524007330)

Program conception, final development and troubleshooting was performed by James D. Burns (Greenblatt and Zhang Labs, UofT), initial development was performed in collaboration with student mentee Ryan Denniston (Zhang Lab, UofT).

## Scripts

### CLIPSHAPE.nf
- **Location:** `icSHAPE_iCLIP/icSHAPE_CLIP_Python/CLIPSHAPE.nf`
- **Description:** Nextflow pipeline that orchestrates the analysis workflow. It manages data input, triggers downstream processing (including MAPSHAPE.py), and handles output aggregation. Ideal for scalable execution.
- **Usage:** Run using Nextflow; ensure Nextflow is installed and properly configured.

### MAPSHAPE.py
- **Location:** `icSHAPE_iCLIP/icSHAPE_CLIP_Python/MAPSHAPE.py`
- **Description:** Python script that processes icSHAPE profiles. It loads in vitro and in vivo experimental data along with RNAplFold and shuffle datasets, combines columns from multiple files, computes median profiles, and generates a comparative plot.
- **Usage:** Run with required command-line arguments:
  ```bash
  python MAPSHAPE.py --pname <ProteinName> --vitroShapeFolder <path> --vivoShapeFolder <path> --plFolder <path> --plshuffle <path>
  ```

### CLIPSHAPE_test.py
- **Location:** `icSHAPE_iCLIP/tests/CLIPSHAPE_test.py`
- **Description:** Unit tests for the CLIPSHAPE module. Tests include module import verification and placeholder for future functionality tests.
- **Usage:** Run via:
  ```bash
  python -m unittest icSHAPE_iCLIP/tests/CLIPSHAPE_test.py
  ```
