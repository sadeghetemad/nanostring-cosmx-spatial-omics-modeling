#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate Clean Violin Plots (Auto Statistical Test, No Warnings)
Each page compares Treated vs Untreated for one metabolic task.

Author: Sadegh Etemad
"""

import os
import gc
import warnings
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import mannwhitneyu, ttest_ind, shapiro

# ---------------- GLOBAL CLEANUP ----------------
warnings.filterwarnings("ignore")
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
os.environ["PYTHONWARNINGS"] = "ignore"

# ---------------- CONFIG ----------------
OUTPUT_DIR = Path("analysis")
DATA_DIR = OUTPUT_DIR / "data"
PDF_DIR = OUTPUT_DIR / "pdf"
PDF_DIR.mkdir(parents=True, exist_ok=True)

MORAN_FILE = DATA_DIR / "All_MoranI_combined.csv"
TASK_TREATED_FILE = DATA_DIR / "Metabolic_Tasks_Treated.csv"
TASK_UNTREATED_FILE = DATA_DIR / "Metabolic_Tasks_Untreated.csv"

# ---------------- LOAD DATA ----------------
moran_df = pd.read_csv(MORAN_FILE)
task_treated = pd.read_csv(TASK_TREATED_FILE)
task_untreated = pd.read_csv(TASK_UNTREATED_FILE)
moran_df.columns = [c.strip() for c in moran_df.columns]

# ---------------- MORAN SUMMARY ----------------
moran_summary = (
    moran_df.groupby(["Treatment_Status", "Task"], as_index=False)["I"]
    .mean()
    .rename(columns={"I": "Mean_MoranI"})
)

# ---------------- FIND COMMON TASKS ----------------
tasks_treated = set(task_treated.drop(columns=["Subject_ID", "Treatment_Status", "Cell_ID","Cell_type"], errors="ignore").columns)
tasks_untreated = set(task_untreated.drop(columns=["Subject_ID", "Treatment_Status", "Cell_ID", "Cell_type"], errors="ignore").columns)
common_tasks = sorted(tasks_treated.intersection(tasks_untreated))

# ---------------- PREPARE OUTPUT FILES ----------------
pdf_all = PdfPages(PDF_DIR / "Violin_All_Tasks.pdf")
pdf_sig = PdfPages(PDF_DIR / "Violin_Significant_Tasks.pdf")
significant_tasks = []

# ---------------- HELPER FUNCTIONS ----------------
def is_normal(data, alpha=0.05):
    """Check normality with Shapiro–Wilk test."""
    if len(data) < 3:
        return False
    stat, p = shapiro(data)
    return p > alpha  # True if normal

# ---------------- MAIN LOOP ----------------
for task in common_tasks:
    treated_vals = task_treated[task].dropna().values
    untreated_vals = task_untreated[task].dropna().values
    if len(treated_vals) < 5 or len(untreated_vals) < 5:
        continue

    # --- Test normality ---
    normal_treated = is_normal(treated_vals)
    normal_untreated = is_normal(untreated_vals)

    # --- Select appropriate test ---
    if normal_treated and normal_untreated:
        test_name = "t-test"
        stat, p_val = ttest_ind(treated_vals, untreated_vals, equal_var=False)
    else:
        test_name = "Mann–Whitney U"
        stat, p_val = mannwhitneyu(treated_vals, untreated_vals, alternative="two-sided")

    # --- Prepare dataframe for plotting ---
    df_plot = pd.DataFrame({
        "Score": np.concatenate([treated_vals, untreated_vals]),
        "Group": ["Treated"] * len(treated_vals) + ["Untreated"] * len(untreated_vals)
    })

    plt.figure(figsize=(6, 6))
    ax = sns.violinplot(
        data=df_plot,
        x="Group",
        y="Score",
        inner="box",
        cut=0,
        linewidth=1.2,
        palette=["#c23b22", "#2255c2"]
    )

    # --- Styling ---
    ymax = df_plot["Score"].max()
    ymin = df_plot["Score"].min()
    plt.ylim(ymin * 0.95, ymax * 1.25)

    plt.ylabel("Metabolic Task Score", fontsize=11)
    plt.xlabel("")
    plt.xticks(fontsize=10)

    # --- Title (separated cleanly from test info) ---
    ax.set_title(task, fontsize=8, loc="center", pad=12, wrap=True)

    # --- Add test result box on top ---
    sig_text = "p < 0.001" if p_val < 0.001 else f"p = {p_val:.3f}"
    color = "darkred" if p_val < 0.05 else "black"
    test_text = f"{test_name} • {sig_text}"

    ax.text(0.5, ymax * 1.18, test_text,
            ha="center", va="bottom",
            fontsize=7, color=color, fontweight="bold")

    # --- Moran I annotations ---
    for group in ["Treated", "Untreated"]:
        mean_I = moran_summary.loc[
            (moran_summary["Treatment_Status"] == group) &
            (moran_summary["Task"] == task),
            "Mean_MoranI"
        ]
        if not mean_I.empty:
            xpos = 0 if group == "Treated" else 1
            ax.text(xpos, ymax * 1.03, f"I={mean_I.values[0]:.3f}",
                    ha="center", va="bottom", fontsize=8, color="darkred")

    plt.tight_layout()
    pdf_all.savefig()
    if p_val < 0.05:
        pdf_sig.savefig()
        significant_tasks.append((task, test_name, p_val, normal_treated, normal_untreated))
    plt.close()
    gc.collect()

pdf_all.close()
pdf_sig.close()

# ---------------- SAVE SUMMARY CSV ----------------
summary_df = pd.DataFrame(significant_tasks,
                          columns=["Task", "Test_Type", "p_value", "Normal_Treated", "Normal_Untreated"])
summary_df.to_csv(DATA_DIR / "DE_Data.csv", index=False)

print("✅ PDFs saved successfully.")
print(f"📘 All Tasks → {PDF_DIR / 'Violin_All_Tasks.pdf'}")
print(f"📕 Significant Tasks → {PDF_DIR / 'Violin_Significant_Tasks.pdf'}")
print(f"📄 CSV Summary → {DATA_DIR / 'DE_Data.csv'}")
