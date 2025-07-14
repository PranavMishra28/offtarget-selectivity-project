# 🧪 Off-Target Selectivity Evaluation Pipeline

This project implements a multi-step pipeline to evaluate the **off-target risks**, **selectivity**, and **binding promiscuity** of small molecules using deep learning, cheminformatics, and public bioactivity databases (ChEMBL, UniProt).

> ✅ Built for drug discovery, chemical biology, and lead optimization workflows.

---

## 🚀 Pipeline Overview

1. **NEBULA** – Generate 3D analogs from a SMILES (via RDKit)
2. **SPARROW** – Estimate synthetic feasibility (Lipinski & complexity-based)
3. **ChEMBL Mining** – Predict empirical targets from real-world bioactivity (MoA & fallback)
4. **Structure Modeling** – Simulate AlphaFold–PLIP binding mock scores
5. **IMPACT Estimator** – Compute selectivity index and safety flag
6. **Metrics & Outputs** – Generate heatmaps, JSON summaries, prioritization scores
7. **Dashboard (Optional)** – Streamlit visualization layer (coming soon)

---

## 📦 Project Structure

```
offtarget_selectivity/
│
├── 01_nebula/                   # 3D generation
├── 02_sparrow/                  # Synthetic prioritization
├── 03_empirical_binding/        # ChEMBL-based off-target prediction
├── 04_structure_modeling/       # Structure-based risk estimation
├── 05_impact_risk/              # Final selectivity and decision logic
├── output_formatter/            # CSV, JSON, heatmap generation
├── target_panel/                # UniProt panel construction
├── model_wrappers/              # DeepPurpose integration
├── dashboard/                   # (Optional) Streamlit interface
├── run.py                       # 🔁 Main runner script
└── requirements.txt             # 📦 Dependencies
```

---

## ⚙️ Setup Instructions

```bash
# 1. Clone the repo
git clone https://github.com/yourusername/offtarget-selectivity.git
cd offtarget-selectivity

# 2. Create a virtual environment
python -m venv .venv
source .venv/Scripts/activate   # On Windows
# or
source .venv/bin/activate       # On macOS/Linux

# 3. Install dependencies
pip install -r requirements.txt

# 4. Download GTEx expression dataset
mkdir -p data/gtex
curl -L \
  "https://zenodo.org/record/4073231/files/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.tsv" \
  -o data/gtex/GTEx_median_tpm.tsv

gzip -f data/gtex/GTEx_median_tpm.tsv
```

---

## ▶️ How to Run

```bash
# Run full pipeline on aspirin (COX-1 as primary target)
python run.py
```

- Input: SMILES string in `run.py`
- Output: Summary metrics in `output/`, plots, JSON, and CSV rankings

---

## 🧠 Features

- DeepPurpose binding predictions (MPNN-CNN)
- Selectivity scoring (on-target vs off-target)
- Promiscuity and synthesis risk flagging
- UniProt panel from curated protein classes (GPCR, Kinase, etc.)
- ChEMBL integration (via official Python client)
- Fallback strategies for robust offline execution

---

## 📊 Example Output

- `offtarget_report.csv` – Ranked list of off-targets
- `offtarget_report.json` – SMILES, predictions, metrics
- `offtarget_report_heatmap.png` – Visual score heatmap
- `impact_summary.json` – Final decision: Reject / Watch / Synthesize

---

## 📚 Dependencies

Major libraries include:

- `rdkit`
- `pandas`, `numpy`, `matplotlib`, `seaborn`
- `DeepPurpose`
- `chembl_webresource_client`
- `requests`
- `scikit-learn`

---

## 💡 Credits

Developed as part of a research and engineering effort for evaluating the off-target landscape of small molecules using **hybrid empirical + predictive** approaches.

Inspired by academic tools like SEA, SwissTargetPrediction, PLIP, and AlphaFold.

---

## 📄 License

MIT License (c) 2025 — Feel free to fork, cite, and contribute!

---

## 🤝 Contributing

Pull requests welcome! Please open an issue first for any major changes.
