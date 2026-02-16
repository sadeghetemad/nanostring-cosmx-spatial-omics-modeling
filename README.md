# Multimodal Analysis of NanoString CosMX Spatial Transcriptomics and Fluxomics

This repository provides a structured research workflow for multimodal
analysis of **NanoString CosMX spatial transcriptomics** data combined
with metabolic task inference, spatial statistics (Moran's I),
colocalization analysis, and graph-based modeling.

This project is organized as a reproducible research workflow composed
of Jupyter notebooks and standalone scripts. It is not designed as a
Python package, but as an analysis pipeline and reference
implementation.

------------------------------------------------------------------------

## Project Structure

    notebooks/

    ├── 0_preliminary/
    │   ├── 01_tma_data_analysis.ipynb
    │   ├── 02_tma_preprocessing.ipynb
    │   ├── 03_apply_sc_cellfie.ipynb
    │   ├── 03_single_cell_cell_communication.ipynb
    │   ├── 03_spatial_cell_cell_communication.ipynb
    │   └── 04_cell_detector.ipynb
    │
    └── 1_latest/
        ├── 1_task_analysis.ipynb
        ├── 2_moran_i_analysis.ipynb
        ├── 3_colocalization_analysis.ipynb
        └── gnn_metabolic_trajectory_analysis.ipynb

    r_scripts/
    ├── preprocess_data.R
    └── sccellfie_results/

    src/
    ├── 0_apply_sccellfie/
    └── 1_compute_task_moran_colocalization.py

    data/
    └── (All raw and processed datasets are stored here)

------------------------------------------------------------------------

## Workflow Overview

### 1. Data Preprocessing (R)

`r_scripts/preprocess_data.R`

-   Loads original NanoString TMA data
-   Reformats data into a structure compatible with Python analysis
-   Saves processed files into the `data/` directory

Run with:

``` bash
Rscript r_scripts/preprocess_data.R
```

------------------------------------------------------------------------

### 2. Preliminary Analysis (notebooks/0_preliminary)

**01_tma_data_analysis.ipynb** - Load and inspect TMA data - Exploratory
data analysis

**02_tma_preprocessing.ipynb** - Cleaning and normalization - Filtering
low-quality genes/cells and QC

**03_apply_sc_cellfie.ipynb** - Apply scCellFie - Infer overlapping
genes, metabolic task activity, and reaction-level scores

**03_single_cell_cell_communication.ipynb** - Single-cell communication
analysis

**03_spatial_cell_cell_communication.ipynb** - Spatially-aware
communication modeling

**04_cell_detector.ipynb** - Binary classification (treated vs
untreated) - Model training and evaluation

------------------------------------------------------------------------

### 3. Core Processing Scripts (src)

#### src/0_apply_sccellfie/

-   Applies scCellFie on each tissue separately
-   Creates separate output folders per subject (folder name = subject
    ID)
-   Outputs overlapping genes, metabolic task scores, and reaction-level
    scores

#### src/1_compute_task_moran_colocalization.py

-   Computes Moran's I for metabolic tasks
-   Computes spatial colocalization
-   Concatenates task data across subjects into a unified CSV
-   Calculates colocalization for the top 5 most significant tasks per
    subject

------------------------------------------------------------------------

### 4. Main Analysis (notebooks/1_latest)

**1_task_analysis.ipynb** - Metabolic task exploration - Cross-subject
comparisons

**2_moran_i\_analysis.ipynb** - Spatial autocorrelation (Moran's I) -
Identification of spatially structured tasks

**3_colocalization_analysis.ipynb** - Task-level colocalization -
Visualization of significant spatial interactions

**gnn_metabolic_trajectory_analysis.ipynb** - Graph-based modeling -
Metabolic trajectory analysis using GNN approaches

------------------------------------------------------------------------

## Installation

Clone the repository:

``` bash
git clone https://github.com/Occhipinti-Lab/nanostring-cosmx-spatial-omics-modeling.git
cd nanostring-cosmx-spatial-omics-modeling
```

Install dependencies:

``` bash
pip install -r requirements.txt
```

Or using conda:

``` bash
conda create -n cosmx python=3.10
conda activate cosmx
pip install -r requirements.txt
```

------------------------------------------------------------------------

## Running the Pipeline

1.  (Optional) Run the R preprocessing script.

2.  Run scCellFie processing in `src/0_apply_sccellfie/`.

3.  Run spatial statistics script:

        python src/1_compute_task_moran_colocalization.py

4.  Open notebooks in `notebooks/1_latest/` for downstream analysis and
    visualization.

Start Jupyter:

``` bash
jupyter lab
```

------------------------------------------------------------------------

## Methods Implemented

-   scCellFie metabolic task inference
-   Reaction-level scoring
-   Moran's I spatial autocorrelation
-   Spatial task colocalization
-   Graph Neural Network (GNN) modeling
-   Cell-state classification

------------------------------------------------------------------------

## Requirements

-   Python 3.9+
-   R (for preprocessing)
-   Jupyter Lab / Notebook
-   See `requirements.txt` for full dependency list

------------------------------------------------------------------------

## License

MIT © 2025 Occhipinti-Lab
