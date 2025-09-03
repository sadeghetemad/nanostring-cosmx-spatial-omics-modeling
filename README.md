# Multimodal analysis of Nanostring CosMX spatial transcriptomics and fluxomics.
 A collection of Jupyter notebooks for multimodal analysis and modeling of Nanostring CosMX spatial transcriptomics and fluxomics data.



Each notebook demonstrates a specific task or workflow.  
This repo is designed as a reference / tutorial / research workflow rather than a single Python package.

---

## 📚 Notebooks

- `01_tma_data_analysis.ipynb` → Read and analyze NanoString tissue microarray (TMA) data
- `02_tma_preprocessing.ipynb` → Perform preprocessing steps on NanoString data
- `03_apply_sc_cellfie.ipynb` → Apply scCellFie for single-cell functional inference to extract reaction- and metabolic-task–level data
- `04_cell_detecotr.ipynb` → Train and Test a binary classifier to detect cell status (treated vs. untreated)


---

## 📜 R Script

- `scripts/preprocess_data.R`  
  This script loads the original NanoString TMA data and converts it into a format that can be read by Python notebooks.  

### Run the R script

```bash
Rscript scripts/preprocess_data.R
```

---

## 📦 Installation

Clone this repository:

```bash
git clone https://github.com/Occhipinti-Lab/nanostring-cosmx-spatial-omics-modeling.git
cd nanostring-cosmx-spatial-omics-modeling
```

Install the dependencies (adjust if you’re using conda/venv):

```bash
pip install -r requirements.txt
```

Or create a conda environment:

```bash
conda create -n myenv python=3.10
conda activate myenv
pip install -r requirements.txt
```

---

## ▶️ Usage

Start Jupyter Lab or Notebook:

```bash
jupyter lab
```

Then open any notebook (`.ipynb`) you want to explore.

---

## 🛠 Requirements

- Python 3.9+  
- Jupyter Lab / Notebook  
- See `requirements.txt` for Python libraries

---

## 🤝 Contributing

Feel free to open issues or pull requests if you’d like to add new notebooks or improve existing ones.

---

## 📄 License

[MIT](LICENSE) © 2025 Ochipinti-Lab