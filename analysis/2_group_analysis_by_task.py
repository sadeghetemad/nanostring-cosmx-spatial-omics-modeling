#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate Clean Violin Plots (Auto Statistical Test, No Warnings)
Each page compares Treated vs Untreated for one metabolic task.

Author: Sadegh Etemad
"""

import os, gc, warnings
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

# ---------------- MORAN SUMMARY (MEAN + STD) ----------------
moran_summary = (
    moran_df.groupby(["Treatment_Status", "Task"], as_index=False)
    .agg(Mean_MoranI=("I", "mean"), STD_MoranI=("I", "std"))
)

# ---------------- FIND COMMON TASKS ----------------
drop_cols = ["Subject_ID", "Treatment_Status", "Cell_ID", "Cell_type"]
tasks_treated = set(task_treated.drop(columns=drop_cols, errors="ignore").columns)
tasks_untreated = set(task_untreated.drop(columns=drop_cols, errors="ignore").columns)
common_tasks = sorted(tasks_treated.intersection(tasks_untreated))

# ---------------- PREPARE OUTPUT FILES ----------------
pdf_all = PdfPages(PDF_DIR / "Violin_All_Tasks.pdf")
pdf_sig = PdfPages(PDF_DIR / "Violin_Significant_Tasks.pdf")
significant_tasks = []

# ---------------- HELPER FUNCTION ----------------
def is_normal(data, alpha=0.05):
    """Check normality with Shapiro–Wilk test."""
    if len(data) < 3:
        return False
    _, p = shapiro(data)
    return p > alpha

# ---------------- MAIN LOOP ----------------
for task in common_tasks:
    treated_vals = task_treated[task].dropna().values
    untreated_vals = task_untreated[task].dropna().values
    if len(treated_vals) < 5 or len(untreated_vals) < 5:
        continue

    normal_treated = is_normal(treated_vals)
    normal_untreated = is_normal(untreated_vals)

    if normal_treated and normal_untreated:
        test_name = "t-test"
        _, p_val = ttest_ind(treated_vals, untreated_vals, equal_var=False)
    else:
        test_name = "Mann–Whitney U"
        _, p_val = mannwhitneyu(treated_vals, untreated_vals, alternative="two-sided")

    df_plot = pd.DataFrame({
        "Score": np.concatenate([treated_vals, untreated_vals]),
        "Group": ["Treated"] * len(treated_vals) + ["Untreated"] * len(untreated_vals)
    })

    plt.figure(figsize=(6, 6))
    ax = sns.violinplot(data=df_plot, x="Group", y="Score", inner="box",
                        cut=0, linewidth=1.2, palette=["#c23b22", "#2255c2"])

    ymax, ymin = df_plot["Score"].max(), df_plot["Score"].min()
    plt.ylim(ymin * 0.95, ymax * 1.25)
    plt.ylabel("Metabolic Task Score", fontsize=11)
    plt.xlabel("")
    plt.xticks(fontsize=10)
    ax.set_title(task, fontsize=8, loc="center", pad=12, wrap=True)

    sig_text = "p < 0.001" if p_val < 0.001 else f"p = {p_val:.3f}"
    color = "darkred" if p_val < 0.05 else "black"
    ax.text(0.5, ymax * 1.18, f"{test_name} • {sig_text}",
            ha="center", va="bottom", fontsize=7, color=color, fontweight="bold")

    for group in ["Treated", "Untreated"]:
        moran_stats = moran_summary.loc[
            (moran_summary["Treatment_Status"] == group) &
            (moran_summary["Task"] == task),
            ["Mean_MoranI", "STD_MoranI"]
        ]
        if not moran_stats.empty:
            xpos = 0 if group == "Treated" else 1
            mean_I = moran_stats["Mean_MoranI"].values[0]
            var_I = moran_stats["STD_MoranI"].values[0]
            ax.text(xpos, ymax * 1.03,
                    f"Moran's I = {mean_I:.3f} ± {var_I:.3f}",
                    ha="center", va="bottom", fontsize=8,
                    color="darkred", fontweight="medium")

    plt.tight_layout()
    pdf_all.savefig()
    if p_val < 0.05:
        pdf_sig.savefig()
        significant_tasks.append((task, test_name, p_val, normal_treated, normal_untreated))
    plt.close()
    gc.collect()

pdf_all.close()
pdf_sig.close()

# ---------------- CREATE ANNOTATED DATASETS ----------------
summary_df = pd.DataFrame(
    significant_tasks,
    columns=["Task", "Test_Type", "p_value", "Normal_Treated", "Normal_Untreated"]
)

task_info_df = pd.read_csv(OUTPUT_DIR / "Task_Info_with_CRC_binary.csv")
summary_df = summary_df.merge(task_info_df, on="Task", how="left")

summary_df = summary_df.sort_values("p_value", ascending=True).reset_index(drop=True)
top20_df = summary_df.head(20).copy()

summary_df.to_csv(DATA_DIR / "DE_Data_Annotated.csv", index=False)
top20_df.to_csv(DATA_DIR / "DE_Top20_Annotated.csv", index=False)

# ---------------- GENERATE PDF FOR TOP 20 ----------------
PDF_TOP20 = PDF_DIR / "Violin_Top20_Tasks.pdf"
pdf_top20 = PdfPages(PDF_TOP20)
top20_tasks = top20_df["Task"].tolist()

print(f"\n📘 Generating violin plots for Top 20 Tasks ({len(top20_tasks)} tasks)...")

for _, row in top20_df.iterrows():
    task = row["Task"]
    test_name = row["Test_Type"]
    p_val = row["p_value"]
    treated_vals = task_treated[task].dropna().values
    untreated_vals = task_untreated[task].dropna().values
    if len(treated_vals) < 5 or len(untreated_vals) < 5:
        continue

    df_plot = pd.DataFrame({
        "Score": np.concatenate([treated_vals, untreated_vals]),
        "Group": ["Treated"] * len(treated_vals) + ["Untreated"] * len(untreated_vals)
    })

    plt.figure(figsize=(6, 6))
    ax = sns.violinplot(data=df_plot, x="Group", y="Score", inner="box",
                        cut=0, linewidth=1.2, palette=["#c23b22", "#2255c2"])
    ymax, ymin = df_plot["Score"].max(), df_plot["Score"].min()
    plt.ylim(ymin * 0.95, ymax * 1.25)
    plt.ylabel("Metabolic Task Score", fontsize=11)
    plt.xlabel("")
    plt.xticks(fontsize=10)
    ax.set_title(task, fontsize=9, loc="center", pad=10, wrap=True)
    ax.text(0.5, ymax * 1.18, f"{test_name} • p = {p_val:.4f}",
            ha="center", va="bottom", fontsize=8, color="darkred", fontweight="bold")

    for group in ["Treated", "Untreated"]:
        moran_stats = moran_summary.loc[
            (moran_summary["Treatment_Status"] == group) &
            (moran_summary["Task"] == task),
            ["Mean_MoranI", "STD_MoranI"]
        ]
        if not moran_stats.empty:
            xpos = 0 if group == "Treated" else 1
            mean_I = moran_stats["Mean_MoranI"].values[0]
            var_I = moran_stats["STD_MoranI"].values[0]
            ax.text(xpos, ymax * 1.03,
                    f"Moran's I = {mean_I:.3f} ± {var_I:.3f}",
                    ha="center", va="bottom", fontsize=8,
                    color="darkred", fontweight="medium")

    plt.tight_layout()
    pdf_top20.savefig()
    plt.close()
    gc.collect()

pdf_top20.close()
print(f"✅ Top 20 Tasks PDF saved → {PDF_TOP20}")

# ---------------- SUMMARY ----------------
print("\n✅ Annotated datasets saved successfully.")
print(f"📄 All Annotated → {DATA_DIR / 'DE_Data_Annotated.csv'}")
print(f"🏆 Top 20 Annotated → {DATA_DIR / 'DE_Top20_Annotated.csv'}\n")
print("🏆 Top 20 Tasks with Most Significant Differences (Annotated):")
print(top20_df[["Task", "p_value", "System", "Subsystem", "Upregulated_in_CRC"]]
      .to_string(index=False))