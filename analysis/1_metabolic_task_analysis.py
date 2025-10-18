#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analyze scCellFie Results (Treated vs Untreated)
(Log-tracked, Warning-Free, Per-Subject Moran + Metabolic Task aggregation)

Author: Sadegh Etemad
"""

import os
import gc
import traceback
import warnings
from pathlib import Path
import scanpy as sc
import pandas as pd
import squidpy as sq
from datetime import datetime

# ---------------- GLOBAL CLEANUP ----------------
warnings.filterwarnings("ignore")
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
os.environ["PYTHONWARNINGS"] = "ignore"

# ---------------- CONFIG ----------------
BASE_DIR = Path().resolve()
RESULTS_DIR = BASE_DIR / "results"
OUTPUT_DIR = BASE_DIR / "analysis"
os.makedirs(OUTPUT_DIR, exist_ok=True)
timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
LOG_PATH = OUTPUT_DIR / f"logs/analysis_log_{timestamp}.txt"

# ---------------- LOGGING ----------------
def log(msg):
    print(msg)
    with open(LOG_PATH, "a", encoding="utf-8") as f:
        f.write(msg + "\n")

def log_error(err_msg):
    print(f"❌ {err_msg}")
    with open(LOG_PATH, "a", encoding="utf-8") as f:
        f.write(f"❌ ERROR: {err_msg}\n")
        traceback.print_exc(file=f)

# ---------------- FUNCTIONS ----------------
def compute_moran_and_tasks(file_path):
    """Compute Moran’s I and extract metabolic tasks for a single h5ad file."""
    try:
        adata = sc.read_h5ad(file_path)

        # ---- Compute spatial neighbors and Moran’s I ----
        sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=6)
        tasks = adata.var_names
        sq.gr.spatial_autocorr(adata, mode="moran", genes=tasks)
        moran = pd.DataFrame(adata.uns["moranI"])
        moran["Subject_ID"] = adata.obs["Subject_ID"].unique()[0]
        moran["Treatment_Status"] = adata.obs["Treatment_Status"].unique()[0]
        moran = moran.reset_index().rename(columns={"index": "Metabolic_Task"})

        # ---- Extract metabolic task matrix ----
        task_df = adata.to_df().copy()
        task_df["Subject_ID"] = adata.obs["Subject_ID"].values
        task_df["Treatment_Status"] = adata.obs["Treatment_Status"].values
        task_df["Cell_ID"] = adata.obs_names
        task_df["Cell_type"] = adata.obs["cell_type"].values

        gc.collect()
        return moran, task_df

    except Exception as e:
        log_error(f"Failed Moran/Task extraction for {file_path}: {e}")
        return pd.DataFrame(), pd.DataFrame()

def load_and_compute_all():
    """Compute Moran’s I and concatenate metabolic task matrices for all subjects."""
    log("📂 Computing Moran’s I and collecting metabolic task data for each subject...")
    all_morans = []
    all_tasks = []

    for folder in RESULTS_DIR.iterdir():
        if not folder.is_dir():
            continue
        files = list(folder.glob("*metabolic_tasks*.h5ad"))
        if not files:
            continue

        f = files[0]
        log(f"🧠 Processing {folder.name}...")
        moran_df, task_df = compute_moran_and_tasks(f)

        if not moran_df.empty:
            log(f"📈 Moran’s I shape for {folder.name}: {moran_df.shape}")
            all_morans.append(moran_df)
        else:
            log(f"⚠️ No Moran’s I data produced for {folder.name}")

        if not task_df.empty:
            log(f"📊 Task matrix shape for {folder.name}: {task_df.shape}")
            all_tasks.append(task_df)
        else:
            log(f"⚠️ No Task matrix data produced for {folder.name}")

    if not all_morans:
        log("⚠️ No valid Moran’s I results found.")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    # Combine all subjects’ data
    moran_all = pd.concat(all_morans, ignore_index=True)
    task_all = pd.concat(all_tasks, ignore_index=True)

    # Save Moran results
    moran_all.to_csv(OUTPUT_DIR / "data/All_MoranI_combined.csv", index=False)
    log(f"✅ Saved combined Moran’s I results: {OUTPUT_DIR / 'All_MoranI_combined.csv'}")

    # Split by treatment
    moran_treated = moran_all[moran_all["Treatment_Status"] == "Treated"]
    moran_untreated = moran_all[moran_all["Treatment_Status"] == "Untreated"]

    task_treated = task_all[task_all["Treatment_Status"] == "Treated"]
    task_untreated = task_all[task_all["Treatment_Status"] == "Untreated"]

    # Save task matrices
    task_treated.to_csv(OUTPUT_DIR / "data/Metabolic_Tasks_Treated.csv", index=False)
    task_untreated.to_csv(OUTPUT_DIR / "data/Metabolic_Tasks_Untreated.csv", index=False)
    log(f"💾 Saved metabolic task matrices for Treated and Untreated groups.")

    log(f"📊 Treated Moran shape: {moran_treated.shape}")
    log(f"📊 Untreated Moran shape: {moran_untreated.shape}")
    log(f"📊 Treated Task matrix shape: {task_treated.shape}")
    log(f"📊 Untreated Task matrix shape: {task_untreated.shape}")

    return True

# ---------------- MAIN ----------------
def main():
    if LOG_PATH.exists():
        LOG_PATH.unlink()

    log("=" * 70)
    log("📊 Analyzing scCellFie Results (Moran + Metabolic Tasks per Subject)")
    log("=" * 70)

    try:
        output = load_and_compute_all()
        if output:
            log("🏁 Analysis completed successfully.")
    except Exception as e:
        log_error(f"Critical failure: {e}")
    finally:
        gc.collect()
        log("=" * 70)

if __name__ == "__main__":
    main()