#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import gc
import traceback
import warnings
from pathlib import Path
import scanpy as sc
import pandas as pd
import squidpy as sq
import numpy as np
import itertools
import sccellfie
from datetime import datetime
from joblib import Parallel, delayed
from scipy.stats import ttest_ind   

# ---------------- SETUP ----------------
warnings.filterwarnings("ignore")
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
os.environ["PYTHONWARNINGS"] = "ignore"

BASE_DIR = Path().resolve().parent
RESULTS_DIR = BASE_DIR / "results"
OUTPUT_DIR = BASE_DIR / "analysis"
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
    print(f" {err_msg}")
    with open(LOG_PATH, "a", encoding="utf-8") as f:
        f.write(f" ERROR: {err_msg}\n")
        traceback.print_exc(file=f)

# ---------------- OPTIMAL NEIGHBOR DETECTION ----------------
def find_optimal_neighbors(adata, k_min=3, k_max=15, test_genes=10):
    try:
        genes = adata.var_names[:test_genes]
        mean_morans = []
        for k in range(k_min, k_max + 1):
            sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=k)
            sq.gr.spatial_autocorr(adata, mode="moran", genes=genes)
            mean_I = pd.DataFrame(adata.uns["moranI"])["I"].mean()
            mean_morans.append((k, mean_I))
        df = pd.DataFrame(mean_morans, columns=["k", "mean_I"])
        df["diff"] = df["mean_I"].diff().abs()

        threshold = 0.005
        stable_idx = df.index[df["diff"].rolling(3).mean() < threshold].min()

        optimal_k = int(df.loc[stable_idx, "k"]) if not np.isnan(stable_idx) else int(df["k"][df["mean_I"].idxmax()])
        return optimal_k

    except Exception as e:
        log_error(f"Optimal neighbor failed: {e}")
        return 6

# ---------------- MORAN + TASK EXTRACTION ----------------
def compute_patient_moran_and_tasks(file_base):
    base_path = file_base / f"{file_base.name}.h5ad"
    task_path = file_base / f"{file_base.name}_metabolic_tasks.h5ad"

    if not base_path.exists() or not task_path.exists():
        log(f" ⚠️ Missing files for {file_base.name}, skipping.")
        return pd.DataFrame(), pd.DataFrame()

    adata = sc.read_h5ad(base_path)
    adata.metabolic_tasks = sc.read_h5ad(task_path)

    subject = adata.obs["Subject_ID"].unique()[0]
    treatment = adata.obs["Treatment_Status"].unique()[0]

    log(f"\n▶ Computing Moran’s I for {subject} ({treatment})")
    optimal_k = find_optimal_neighbors(adata, 3, 15, 15)

    adata.metabolic_tasks.obsm["spatial"] = adata.obsm["spatial"]

    sq.gr.spatial_neighbors(adata.metabolic_tasks, coord_type="generic", n_neighs=optimal_k)
    sq.gr.spatial_autocorr(adata.metabolic_tasks, mode="moran", genes=adata.metabolic_tasks.var_names)

    moran = pd.DataFrame(adata.metabolic_tasks.uns["moranI"])
    moran = moran.reset_index().rename(columns={"index": "Task"})
    moran["Subject_ID"] = subject
    moran["Treatment_Status"] = treatment
    moran["Optimal_k"] = optimal_k

    if "I" in moran.columns:
        moran["I_z"] = (moran["I"] - moran["I"].mean()) / moran["I"].std(ddof=0)

    task_df = adata.metabolic_tasks.to_df().copy()
    task_df["Subject_ID"] = subject
    task_df["Treatment_Status"] = treatment
    task_df["Cell_ID"] = adata.obs_names
    task_df["Cell_type"] = adata.obs.get("cell_type", "NA").values

    spatial = adata.metabolic_tasks.obsm["spatial"]
    task_df["x"] = spatial[:, 0]
    task_df["y"] = spatial[:, 1]

    del adata
    gc.collect()
    return moran, task_df

# ---------------- SAFE PARALLEL COLOCALIZATION ----------------
def compute_pair_score(task_matrix, neigh_list, subject, treatment, t1, t2, t_idx):
    try:
        x = task_matrix[:, t_idx[t1]]
        y = task_matrix[:, t_idx[t2]]

        scores = []
        for i in range(len(x)):
            neigh = neigh_list[i]
            if len(neigh) < 3:
                continue
            score = np.mean((x[neigh] - x[neigh].mean()) * (y[neigh] - y[neigh].mean()))
            scores.append(score)

        mean_score = float(np.mean(scores)) if len(scores) > 0 else 0.0
        return subject, treatment, t1, t2, mean_score

    except Exception as e:
        return subject, treatment, t1, t2, "ERROR:" + str(e)

def compute_colocalization_for_filtered_tasks(file_base, selected_tasks, coloc_out_path):
    base_path = file_base / f"{file_base.name}.h5ad"
    task_path = file_base / f"{file_base.name}_metabolic_tasks.h5ad"

    adata = sc.read_h5ad(base_path)
    adata_tasks = sc.read_h5ad(task_path)

    subject = adata.obs["Subject_ID"].unique()[0]
    treatment = adata.obs["Treatment_Status"].unique()[0]

    log(f"\n▶ Running PARALLEL colocalization for {subject} ({treatment})")

    available_tasks = list(adata_tasks.var_names)
    selected_tasks = [t for t in selected_tasks if t in available_tasks]

    if len(selected_tasks) < 2:
        log(f" ⚠️ Not enough tasks for {subject}, skipping colocalization.")
        return

    # Build KNN graph ONCE
    sq.gr.spatial_neighbors(adata_tasks, coord_type="generic", n_neighs=10)
    adj = adata_tasks.obsp["spatial_connectivities"].toarray()
    adj = (adj > 0).astype(int)

    neigh_list = [np.where(adj[i] == 1)[0] for i in range(adj.shape[0])]

    # Task matrix
    df_tasks = adata_tasks.to_df()[selected_tasks].astype(float)
    task_matrix = df_tasks.to_numpy()

    t_idx = {t: i for i, t in enumerate(selected_tasks)}

    pairs = list(itertools.combinations(selected_tasks, 2))

    results = Parallel(n_jobs=15, backend="loky", verbose=3)(
        delayed(compute_pair_score)(
            task_matrix, neigh_list, subject, treatment, t1, t2, t_idx
        )
        for t1, t2 in pairs
    )

    with open(coloc_out_path, "a", encoding="utf-8") as f:
        for r in results:
            subject, treatment, t1, t2, score = r
            f.write(f"{subject},{treatment},{t1},{t2},{score}\n")

    del adata, adata_tasks, task_matrix, df_tasks
    gc.collect()

# ---------------- MAIN PIPELINE ----------------
def main():
    coloc_out_path = OUTPUT_DIR / "data/All_Colocalization_Scores.csv"
    with open(coloc_out_path, "w", encoding="utf-8") as f:
        f.write("Subject_ID,Treatment_Status,Task_1,Task_2,Colocalization_Score\n")

    all_morans = []
    all_tasks = []

    # Step 1 — compute morans
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

    moran_all.to_csv(OUTPUT_DIR / "data/All_MoranI_combined.csv", index=False)
    task_all.to_csv(OUTPUT_DIR / "data/All_Metabolic_Tasks.csv", index=False)

    log(" ✅ Saved treated and untreated metabolic task matrices.")

    #  t-test task selection 
    treated_all = task_all[task_all["Treatment_Status"] == "Treated"]
    untreated_all = task_all[task_all["Treatment_Status"] == "Untreated"]

    task_cols = [
        c for c in task_all.columns
        if c not in ["Subject_ID", "Treatment_Status", "Cell_ID", "Cell_type"]
    ]

    results = []

    for task in task_cols:
        t_vals = treated_all[task].values
        u_vals = untreated_all[task].values

        if len(t_vals) < 5 or len(u_vals) < 5:
            continue

        # Welch's t-test
        stat, p = ttest_ind(t_vals, u_vals, equal_var=False, nan_policy="omit")

        results.append((task, p))

    if len(results) == 0:
        log("\n⚠️ No tasks passed t-test filtering. Check data integrity.")
        sig_tasks = []
    else:
        res_df = pd.DataFrame(results, columns=["Task", "p_value"])
        res_df = res_df.sort_values("p_value")

        sig_tasks = res_df.head(5)["Task"].tolist()

    log(f"\n📊 Total significant tasks selected via t-test: {len(sig_tasks)}")
    log(f" Selected significant tasks: {sig_tasks}")

    # Step 3 — colocalization
    for folder in RESULTS_DIR.iterdir():
        if not folder.is_dir():
            continue
        compute_colocalization_for_filtered_tasks(folder, sig_tasks, coloc_out_path)

    log("\n✅ Full pipeline completed successfully.")

if __name__ == "__main__":
    main()
