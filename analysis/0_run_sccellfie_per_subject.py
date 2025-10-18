#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Run scCellFie per Subject
(Simple, CPU-Optimized, with tqdm progress bar)

Author: Sadegh Etemad
"""

import os
import gc
import warnings
from pathlib import Path
from tqdm import tqdm
import scanpy as sc
import sccellfie

# Silence warnings
warnings.filterwarnings("ignore")
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
os.environ["PYTHONWARNINGS"] = "ignore"

# ---------------- CONFIG ----------------
BASE_DIR = Path().resolve()
INPUT_PATH = BASE_DIR / "data/h5ad/filtered_normalized_data.h5ad"
OUTPUT_DIR = BASE_DIR / "results"
ORGANISM = "human"

os.makedirs(OUTPUT_DIR, exist_ok=True)
LOG_PATH = OUTPUT_DIR / "run_log.txt"

# ---------------- SIMPLE LOGGING ----------------
def log(msg):
    """Print message and append to log file."""
    print(msg)
    with open(LOG_PATH, "a", encoding="utf-8") as f:
        f.write(msg + "\n")

# ---------------- FUNCTIONS ----------------
def run_sccellfie(adata, folder, name):
    """Run scCellFie pipeline for one subject."""
    os.makedirs(folder, exist_ok=True)
    log(f"🧬 Running scCellFie for {name}...")
    adata_out = sccellfie.run_sccellfie_pipeline(
        adata,
        organism=ORGANISM,
        n_counts_col='nCount_Nanostring',
        process_by_group=False,
        groupby=None,
        neighbors_key='neighbors',
        n_neighbors=10,
        batch_key='sample',
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
    log("🚀 Running scCellFie per subject")
    log("=" * 60)

    # Load dataset
    adata = sc.read_h5ad(INPUT_PATH)
    adata.var_names = adata.var["gene"].astype(str)
    adata.var_names_make_unique()

    # Basic validation
    if "Subject_ID" not in adata.obs or "Treatment_Status" not in adata.obs:
        raise ValueError("❌ Missing columns: Subject_ID or Treatment_Status")

    subjects = adata.obs["Subject_ID"].unique()

    # tqdm progress bar
    log(f"📦 Found {len(subjects)} subjects")
    for sid in tqdm(subjects, desc="Processing subjects", ncols=90):
        subset = adata[adata.obs["Subject_ID"] == sid].copy()
        treatment = subset.obs["Treatment_Status"].unique()[0]
        out_dir = OUTPUT_DIR / sid

        try:
            run_sccellfie(subset, out_dir, sid)
            log(f"✅ Done: {sid} ({treatment})")
        except Exception as e:
            log(f"⚠️ Error in {sid}: {e}")

    log("🏁 All subjects processed successfully.")
    log("=" * 60)


if __name__ == "__main__":
    main()