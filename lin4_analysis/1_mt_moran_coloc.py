#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import warnings
warnings.filterwarnings("ignore")

import os
import gc
import traceback
from pathlib import Path
import scanpy as sc
import pandas as pd
import squidpy as sq
import numpy as np
import itertools
from datetime import datetime
from joblib import Parallel, delayed
from scipy.stats import ttest_ind   

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
os.environ["PYTHONWARNINGS"] = "ignore"

# ---------------- DIR SETUP ----------------
BASE_DIR = Path().resolve()
RESULTS_DIR = BASE_DIR / "sccellfie_results_new"
OUTPUT_DIR = BASE_DIR
os.makedirs(OUTPUT_DIR / "data", exist_ok=True)
os.makedirs(OUTPUT_DIR / "logs", exist_ok=True)
timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
LOG_PATH = OUTPUT_DIR / f"logs/analysis_log_{timestamp}.txt"


# ---------------- LOGGING ----------------
def log(msg):
    print(msg)
    with open(LOG_PATH, "a", encoding="utf-8") as f:
        f.write(msg + "\n")


def log_error(err_msg):
    print(f"{err_msg}")
    with open(LOG_PATH, "a", encoding="utf-8") as f:
        f.write(f"ERROR: {err_msg}\n")
        traceback.print_exc(file=f)


def use_fixed_neighbors(adata):
    """Attach spatial neighbors using existing 30-neighbor connectivities."""
    if "neighbors_30" not in adata.uns:
        raise ValueError("neighbors_30 is missing in this dataset.")

    adata.uns["spatial_neighbors"] = adata.uns["neighbors_30"]
    adata.obsp["spatial_connectivities"] = adata.obsp["neighbors_30_connectivities"]
    adata.obsp["spatial_distances"] = adata.obsp["neighbors_30_distances"]
    return adata


# -------------------------------------------------------------------
#  MORAN + METABOLIC TASK EXTRACTION
# -------------------------------------------------------------------
def compute_patient_moran_and_tasks(file_base):

    base_path = file_base / f"{file_base.name}.h5ad"
    task_path = file_base / f"{file_base.name}_metabolic_tasks.h5ad"

    if not base_path.exists() or not task_path.exists():
        log(f"⚠ Missing files for {file_base.name}, skipping.")
        return pd.DataFrame(), pd.DataFrame()

    adata = sc.read_h5ad(base_path)
    adata_tasks = sc.read_h5ad(task_path)

    subject = adata.obs["Subject_ID"].unique()[0]
    treatment = adata.obs["Treatment_Status"].unique()[0]

    log(f"\n▶ Moran’s I: {subject} ({treatment}) using k=30 fixed neighbors")

    # Attach 30-neighbor spatial graph
    adata_tasks.obsm["spatial"] = adata.obsm["spatial"]
    adata_tasks = use_fixed_neighbors(adata_tasks)

    # Moran’s I for all metabolic tasks
    sq.gr.spatial_autocorr(
        adata_tasks,
        mode="moran",
        genes=adata_tasks.var_names
    )

    moran = pd.DataFrame(adata_tasks.uns["moranI"])
    moran = moran.reset_index().rename(columns={"index": "Task"})
    moran["Subject_ID"] = subject
    moran["Treatment_Status"] = treatment
    moran["Optimal_k"] = 30

    # z-normalize
    if "I" in moran.columns:
        moran["I_z"] = (moran["I"] - moran["I"].mean()) / moran["I"].std(ddof=0)

    # Extract task matrix with spatial coordinates
    task_df = adata_tasks.to_df().copy()
    task_df["Subject_ID"] = subject
    task_df["Treatment_Status"] = treatment
    task_df["Cell_ID"] = adata.obs_names

    # UPDATED: correct cell type annotation
    task_df["Cell_type"] = adata.obs["Lineage_level4"].values

    spatial = adata_tasks.obsm["spatial"]
    task_df["x"] = spatial[:, 0]
    task_df["y"] = spatial[:, 1]


    del adata, adata_tasks
    gc.collect()
    return moran, task_df


# -------------------------------------------------------------------
#  COLOCALIZATION USING FIXED k=30 NEIGHBORS
# -------------------------------------------------------------------
def compute_pair_score(task_matrix, neigh_list, subject, treatment, t1, t2, t_idx):
    try:
        x = task_matrix[:, t_idx[t1]]
        y = task_matrix[:, t_idx[t2]]

        vals = []
        for i in range(len(x)):
            neigh = neigh_list[i]
            if len(neigh) == 0:
                continue
            vals.append(np.mean((x[neigh] - x[neigh].mean()) * (y[neigh] - y[neigh].mean())))

        score = float(np.mean(vals)) if len(vals) > 0 else 0.0
        return subject, treatment, t1, t2, score

    except Exception as e:
        return subject, treatment, t1, t2, f"ERROR:{str(e)}"


def compute_colocalization_for_filtered_tasks(file_base, selected_tasks, coloc_out_path):

    
    base_path = file_base / f"{file_base.name}.h5ad"
    task_path = file_base / f"{file_base.name}_metabolic_tasks.h5ad"

    if not base_path.exists() or base_path.stat().st_size == 0:
        log(f"⚠ Skipping {file_base.name}: base file missing or empty.")
        return

    if not task_path.exists() or task_path.stat().st_size == 0:
        log(f"⚠ Skipping {file_base.name}: metabolic task file missing or empty.")
        return

    try:
        adata = sc.read_h5ad(base_path)
        adata_tasks = sc.read_h5ad(task_path)
    except Exception as e:
        log(f"⚠ Skipping {file_base.name}: failed to load h5ad ({e})")
        return

    subject = adata.obs["Subject_ID"].unique()[0]
    treatment = adata.obs["Treatment_Status"].unique()[0]

    log(f"\n▶ Colocalization: {subject} ({treatment}) using 30-neighbor graph")

    # Attach spatial & 30-neighbor graph
    adata_tasks.obsm["spatial"] = adata.obsm["spatial"]
    adata_tasks = use_fixed_neighbors(adata_tasks)

    available_tasks = list(adata_tasks.var_names)
    selected_tasks = [t for t in selected_tasks if t in available_tasks]

    if len(selected_tasks) < 2:
        log(f"⚠ Not enough tasks for {subject}, skipping.")
        return

    # ---------- FIXED VERSION: use sparse graph ----------
    adj = adata_tasks.obsp["spatial_connectivities"].tocsr()
    neigh_list = [adj[i].indices for i in range(adj.shape[0])]

    # Task matrix
    df_tasks = adata_tasks.to_df()[selected_tasks].astype(float)
    task_matrix = df_tasks.to_numpy()

    t_idx = {t: i for i, t in enumerate(selected_tasks)}
    pairs = list(itertools.combinations(selected_tasks, 2))

    results = Parallel(n_jobs=10, backend="loky", verbose=3)(
        delayed(compute_pair_score)(
            task_matrix, neigh_list, subject, treatment, t1, t2, t_idx
        )
        for t1, t2 in pairs
    )

    with open(coloc_out_path, "a", encoding="utf-8") as f:
        for r in results:
            s, tstat, t1, t2, score = r
            f.write(f"{s},{tstat},{t1},{t2},{score}\n")

    del adata, adata_tasks
    gc.collect()

# -------------------------------------------------------------------
#  MAIN PIPELINE
# -------------------------------------------------------------------
def main():
    coloc_out_path = OUTPUT_DIR / "data/All_Colocalization_Scores.csv"
    with open(coloc_out_path, "w", encoding="utf-8") as f:
        f.write("Subject_ID,Treatment_Status,Task_1,Task_2,Colocalization_Score\n")

    moran_file = OUTPUT_DIR / "data/All_MoranI_combined.csv"
    task_file = OUTPUT_DIR / "data/All_Metabolic_Tasks.csv"

    # -------------------------------------------------------------
    #  STEP 1 — CHECK IF MORAN & TASK MATRICES ALREADY EXIST
    # -------------------------------------------------------------
    if moran_file.exists() and task_file.exists():
        log("⚠ Files already exist. Skipping Moran + Task extraction.")
        moran_all = pd.read_csv(moran_file)
        task_all = pd.read_csv(task_file)
    else:
        log("▶ Running Moran + Task extraction... (files not found)")
        all_morans, all_tasks = [], []

        for folder in RESULTS_DIR.iterdir():
            if not folder.is_dir():
                continue
            moran, task_df = compute_patient_moran_and_tasks(folder)
            if not moran.empty:
                all_morans.append(moran)
            if not task_df.empty:
                all_tasks.append(task_df)

        moran_all = pd.concat(all_morans, ignore_index=True)
        task_all = pd.concat(all_tasks, ignore_index=True)

        moran_all.to_csv(moran_file, index=False)
        task_all.to_csv(task_file, index=False)

        log("✅ Saved global Moran and metabolic task matrices.")

    # -------------------------------------------------------------
    #  STEP 2 — T-TEST FOR SIGNIFICANT TASKS
    # -------------------------------------------------------------
    treated = task_all[task_all["Treatment_Status"] == "Treated"]
    untreated = task_all[task_all["Treatment_Status"] == "Untreated"]

    task_cols = [
        c for c in task_all.columns
        if c not in ["Subject_ID", "Treatment_Status", "Cell_ID", "Cell_type", "x", "y"]
    ]

    results = []
    for task in task_cols:
        t_vals = treated[task].values
        u_vals = untreated[task].values
        if len(t_vals) > 5 and len(u_vals) > 5:
            stat, p = ttest_ind(t_vals, u_vals, equal_var=False, nan_policy="omit")
            results.append((task, p))

    if not results:
        log("⚠ No tasks significant; skipping colocalization.")
        sig_tasks = []
    else:
        df = pd.DataFrame(results, columns=["Task", "p_value"]).sort_values("p_value")
        sig_tasks = df.head(5)["Task"].tolist()

    log(f"\n📊 Significant tasks to test colocalization: {sig_tasks}")

    # -------------------------------------------------------------
    #  STEP 3 — COLOCALIZATION ONLY
    # -------------------------------------------------------------
    for folder in RESULTS_DIR.iterdir():
        if not folder.is_dir():
            continue
        compute_colocalization_for_filtered_tasks(folder, sig_tasks, coloc_out_path)

    log("\n🎉 Full pipeline completed successfully.")


if __name__ == "__main__":
    main()
