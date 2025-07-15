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
7. **Dashboard** – Aggregated summary of all risk metrics, expression scores, and final recommendations
8. **Conflict & Redundancy Resolution (Optional Task)** - conflict resolution & redundancy checking logic between AlphaFold vs Empirical
   models

---

## ⚙️ Setup Instructions

### 1. Clone the repo

```bash
git clone https://github.com/yourusername/offtarget-selectivity.git
cd offtarget-selectivity
```

### 2. Create a virtual environment

```bash
python -m venv .venv
source .venv/Scripts/activate   # On Windows
# or
source .venv/bin/activate       # On macOS/Linux
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

### 4. Download GTEx expression dataset

```bash
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

## 📊 Output Files

### 🧪 Task 1: 🌌 NEBULA (Library Generation)

- `/nebula/generated_library.sdf`
  → SDF file of structurally similar analogs generated from input compound

### 🧬 Task 2: 🧪 SPARROW (Synthetic Accessibility)

- `/sparrow/triaged_library.csv`
  → Molecules ranked by synthetic difficulty and prioritization score

### 📈 Task 3: Empirical Binding (ChEMBL)

- `/empirical_binding/offtarget_predictions.json`
  → Empirically associated off-target UniProt IDs and confidence scores

### 🧠 Task 4: Structural Binding Risk (AlphaFold + PLIP)

- `/structure_modeling/binding_risk.json`
  → Structural risk estimate (mocked or real docking results per UniProt ID)

### ⚖️ Task 5: Selectivity & Safety (IMPACT)

- `/impact_risk/impact_summary.json`
  → Combined risk summary with decision flag: **Synthesize**, **Watch**, or **Reject**

### 🧬 Task 6: Expression-Aware Risk Reweighting (GTEx)

- `/expression_filter/tissue_weighted_risk.json`
  → Risk scores adjusted based on average GTEx tissue expression

### 📊 Task 7: Final Dashboard Summary

- `/final_dashboard/compound_summary.json`
  → Aggregated summary of all risk metrics, expression scores, and final recommendations

### 🔀 Optional Task 8: Conflict & Redundancy Resolution

- `/conflict_resolution/conflict_summary.json`
  → Flags for targets with model disagreement (`conflict`), similarity (`redundant`), or agreement

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
