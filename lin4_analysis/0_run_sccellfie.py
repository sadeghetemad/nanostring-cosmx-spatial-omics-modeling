#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import warnings
warnings.filterwarnings("ignore")


import os
import gc
from pathlib import Path
from tqdm import tqdm
import scanpy as sc
import sccellfie

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
os.environ["PYTHONWARNINGS"] = "ignore"

# ---------------- CONFIG ----------------
BASE_DIR = Path().resolve()
INPUT_PATH = BASE_DIR / "CRC_lin4_spatial.h5ad"
OUTPUT_DIR = BASE_DIR / "sccellfie_results_new"
ORGANISM = "human"

os.makedirs(OUTPUT_DIR, exist_ok=True)
LOG_PATH = OUTPUT_DIR / "run_log.txt"

# ---------------- SIMPLE LOGGING ----------------
def log(msg):
    print(msg)
    with open(LOG_PATH, "a", encoding="utf-8") as f:
        f.write(msg + "\n")

# ---------------- FUNCTIONS ----------------
def run_sccellfie(adata, folder, name):
    """Run scCellFie pipeline for one subject."""
    os.makedirs(folder, exist_ok=True)
    log(f"Running scCellFie for {name}...")

    adata_out = sccellfie.run_sccellfie_pipeline(
        adata,
        organism=ORGANISM,
        n_counts_col="total_counts",
        process_by_group=False,
        groupby=None,
        neighbors_key="neighbors_30",  
        n_neighbors=30,     
        batch_key="Subject_ID",
        smooth_cells=True,
        alpha=0.33,
        chunk_size=5000,
        disable_pbar=True,
        save_folder=folder,
        save_filename=name,
    )

    del adata
    gc.collect()
    return adata_out


def main():
    log("=" * 60)
    log("Running scCellFie per subject")
    log("=" * 60)

    # Load dataset
    adata = sc.read_h5ad(INPUT_PATH)

    # Ensure var names are string & unique
    adata.var_names = [str(i) for i in adata.var_names]
    adata.var_names_make_unique()

    # Validate metadata availability
    for col in ["Subject_ID", "Treatment_Status"]:
        if col not in adata.obs:
            raise ValueError(f"Missing column: {col}")

    subjects = adata.obs["Subject_ID"].unique()
    log(f" Found {len(subjects)} subjects")

    # Process each subject separately
    for sid in tqdm(subjects, desc="Processing subjects", ncols=90):
        subset = adata[adata.obs["Subject_ID"] == sid].copy()
        treatment = subset.obs["Treatment_Status"].unique()[0]
        out_dir = OUTPUT_DIR / sid

        try:
            run_sccellfie(subset, out_dir, sid)
            log(f" Done: {sid} (Treatment: {treatment})")
        except Exception as e:
            log(f" Error in {sid}: {e}")

        del subset
        gc.collect()

    log(" All subjects processed successfully.")
    log("=" * 60)


if __name__ == "__main__":
    main()
