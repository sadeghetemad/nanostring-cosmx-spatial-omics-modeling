#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare Treated vs Untreated per Cell Type (Auto Test + Violin Plot)
Each page shows one cell type comparison instead of metabolic tasks.

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

TASK_TREATED_FILE = DATA_DIR / "Metabolic_Tasks_Treated.csv"
TASK_UNTREATED_FILE = DATA_DIR / "Metabolic_Tasks_Untreated.csv"

# ---------------- LOAD DATA ----------------
task_treated = pd.read_csv(TASK_TREATED_FILE)
task_untreated = pd.read_csv(TASK_UNTREATED_FILE)

# Check that Cell_type column exists
if "Cell_type" not in task_treated.columns or "Cell_type" not in task_untreated.columns:
    raise ValueError("❌ Column 'Cell_type' not found in input CSV files.")

# ---------------- COMPUTE MEAN SCORE PER CELL ----------------
# Drop non-feature columns
drop_cols = ["Subject_ID", "Treatment_Status", "Cell_ID", "Cell_type"]
task_cols = [c for c in task_treated.columns if c not in drop_cols]

# Compute mean metabolic activity per cell (average of all tasks)
task_treated["Mean_Metabolic_Activity"] = task_treated[task_cols].mean(axis=1)
task_untreated["Mean_Metabolic_Activity"] = task_untreated[task_cols].mean(axis=1)

# ---------------- FIND COMMON CELL TYPES ----------------
celltypes_treated = set(task_treated["Cell_type"].unique())
celltypes_untreated = set(task_untreated["Cell_type"].unique())
common_celltypes = sorted(celltypes_treated.intersection(celltypes_untreated))

# ---------------- PREPARE OUTPUT FILES ----------------
pdf_all = PdfPages(PDF_DIR / "Violin_All_CellTypes.pdf")
pdf_sig = PdfPages(PDF_DIR / "Violin_Significant_CellTypes.pdf")
significant_cells = []

# ---------------- HELPER FUNCTIONS ----------------
def is_normal(data, alpha=0.05):
    """Check normality with Shapiro–Wilk test."""
    if len(data) < 3:
        return False
    stat, p = shapiro(data)
    return p > alpha  # True if normal

# ---------------- MAIN LOOP ----------------
for cell_type in common_celltypes:
    treated_vals = task_treated.loc[task_treated["Cell_type"] == cell_type, "Mean_Metabolic_Activity"].dropna().values
    untreated_vals = task_untreated.loc[task_untreated["Cell_type"] == cell_type, "Mean_Metabolic_Activity"].dropna().values

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

    plt.ylabel("Mean Metabolic Task Score", fontsize=11)
    plt.xlabel("")
    plt.xticks(fontsize=10)

    # --- Title ---
    ax.set_title(f"{cell_type}", fontsize=9, loc="center", pad=15, wrap=True)

    # --- Add test result info ---
    sig_text = "p < 0.001" if p_val < 0.001 else f"p = {p_val:.3f}"
    color = "darkred" if p_val < 0.05 else "black"
    test_text = f"{test_name} • {sig_text}"

    ax.text(0.5, ymax * 1.17, test_text,
            ha="center", va="bottom",
            fontsize=9, color=color, fontweight="bold")

    # --- Add means below violins ---
    mean_t = treated_vals.mean()
    std_t = treated_vals.std()
    mean_u = untreated_vals.mean()
    std_u = untreated_vals.std()
    ax.text(0, ymin * 0.97, f"{mean_t:.2f} ± {std_t:.2f}", ha="center", va="top", fontsize=7, color="#c23b22")
    ax.text(1, ymin * 0.97, f"{mean_u:.2f} ± {std_u:.2f}", ha="center", va="top", fontsize=7, color="#2255c2")

    plt.tight_layout()
    pdf_all.savefig()
    if p_val < 0.05:
        pdf_sig.savefig()
        significant_cells.append((cell_type, test_name, p_val, normal_treated, normal_untreated, mean_t, std_t, mean_u, std_u))
    plt.close()
    gc.collect()

pdf_all.close()
pdf_sig.close()

# ---------------- SAVE SUMMARY CSV ----------------
summary_df = pd.DataFrame(significant_cells,
    columns=["Cell_type", "Test_Type", "p_value", "Normal_Treated", "Normal_Untreated", "Mean_Treated", "STD_Treated", "Mean_Untreated", "STD_Untreated"]
)
summary_df.to_csv(DATA_DIR / "DE_CellType.csv", index=False)

print("✅ PDFs saved successfully.")
print(f"📘 All Cell Types → {PDF_DIR / 'Violin_All_CellTypes.pdf'}")
print(f"📕 Significant Cell Types → {PDF_DIR / 'Violin_Significant_CellTypes.pdf'}")
print(f"📄 CSV Summary → {DATA_DIR / 'DE_CellType.csv'}")