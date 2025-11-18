import warnings
warnings.filterwarnings("ignore")

import scanpy as sc
import matplotlib.pyplot as plt
import sccellfie
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path
import pandas as pd

# ============================
# CONFIGURATION
# ============================
BASE_DIR = Path().resolve()
RESULTS_DIR = BASE_DIR / "results"
OUTPUT_DIR = BASE_DIR / "analysis"
DATA_DIR = OUTPUT_DIR / "data"
PDF_DIR = OUTPUT_DIR / "pdf"

PDF_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)

TOP20_TASKS_FILE = DATA_DIR / "DE_Top20_Annotated.csv"
CELLTYPE_FILE = DATA_DIR / "DE_CellType.csv"

spot_size = 200
cmap = "YlGnBu"

# ============================
# LOAD TASKS & CELL TYPES
# ============================
print("Loading task and cell type lists...")

top20_df = pd.read_csv(TOP20_TASKS_FILE)
top_tasks = top20_df["Task"].dropna().unique().tolist()[:2]

celltype_df = pd.read_csv(CELLTYPE_FILE)
cell_types = celltype_df["Cell_type"].dropna().unique().tolist()

print(f"Top 20 tasks: {top_tasks}")
print(f"Cell types: {cell_types}\n")

# ============================
# DEFINE LOG FUNCTION
# ============================
def log(msg: str):
    print(msg, flush=True)

# ============================
# LOAD DATA AND SEPARATE GROUPS
# ============================
treated_patients = []
untreated_patients = []

for folder in sorted(RESULTS_DIR.iterdir()):
    if not folder.is_dir():
        continue

    files = list(folder.glob("*metabolic_tasks*.h5ad"))
    if not files:
        continue

    f = files[0]
    patient_id = folder.name
    log(f"🧠 Loading patient {patient_id} ...")

    try:
        adata = sc.read_h5ad(f)
    except Exception as e:
        log(f"❌ Error reading {patient_id}: {e}")
        continue

    if "Treatment_Status" not in adata.obs.columns:
        log(f"⚠️ No 'Treatment_Status' column in {patient_id}, skipping.")
        continue

    # Determine patient type
    status_values = adata.obs["Treatment_Status"].str.lower().unique().tolist()
    if "treated" in status_values:
        treated_patients.append((patient_id, adata))
    elif "untreated" in status_values:
        untreated_patients.append((patient_id, adata))
    else:
        log(f"⚠️ Unknown Treatment_Status in {patient_id}: {status_values}")

log(f"\nTotal Treated Patients: {len(treated_patients)}")
log(f"Total Untreated Patients: {len(untreated_patients)}\n")

# ============================
# DEFINE PLOT FUNCTION
# ============================
def plot_group_pdf(group_name, patients_list):
    pdf_path = PDF_DIR / f"{group_name.upper()}_AllPatients_Top20Tasks.pdf"
    log(f"📄 Creating {pdf_path} ...")

    with PdfPages(pdf_path) as pdf:
        for task in top_tasks:
            for patient_id, adata in patients_list:
                # if not hasattr(adata, "metabolic_tasks"):
                #     log(f"⚠️ No 'metabolic_tasks' for {patient_id}, skipping.")
                #     continue

                if task not in adata.var_names:
                    log(f"⚠️ Task '{task}' not found in {patient_id}, skipping.")
                    continue

                # ---- Plot both maps side-by-side ----
                log(f"   Plotting {patient_id} - {group_name} - {task}")
                fig, axes = plt.subplots(1, 2, figsize=(12, 6))

                adata = adata[adata.obs["cell_type"].isin(cell_types)].copy()

                try:
                    # Left: Cell type map
                    sccellfie.plotting.plot_spatial(
                        adata,
                        keys=["cell_type"],
                        spot_size=spot_size,
                        use_raw=False,
                        ncols=1,
                        hspace=0.3,
                        save=False,
                        ax=axes[0]
                    )
                    axes[0].set_title("Cell Type", fontsize=12)

                    # Right: Metabolic task map
                    sccellfie.plotting.plot_spatial(
                        adata,
                        keys=[task],
                        spot_size=spot_size,
                        cmap=cmap,
                        use_raw=False,
                        ncols=1,
                        hspace=0.3,
                        vmin=0,
                        save=False,
                        ax=axes[1]
                    )
                    axes[1].set_title(f"Task: {task}", fontsize=12)

                    fig.suptitle(f"{patient_id} ({group_name.capitalize()})", fontsize=14)
                    fig.tight_layout()
                    pdf.savefig(fig)
                    plt.close(fig)

                except Exception as e:
                    log(f"   ⚠️ Error plotting {patient_id}-{task}: {e}")

    log(f"✅ Saved {pdf_path}\n")

# ============================
# CREATE PDFs
# ============================
if treated_patients:
    treated_patients = ['92604']
    plot_group_pdf("treated", treated_patients)

if untreated_patients:
    untreated_patients = ['1209622B']
    plot_group_pdf("untreated", untreated_patients)

print("🎯 All done successfully.")