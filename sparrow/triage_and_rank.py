# offtarget_selectivity/02_sparrow/triage_and_rank.py

import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors


def compute_prioritization(smiles_list, output_path="sparrow/ranked_candidates.csv"):
    scored = []

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"⚠️ Skipping invalid SMILES: {smiles}")
            continue

        try:
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            score = 1 / (mw * (logp + 1))  # Simplified "ease" score
            scored.append((smiles, mw, logp, score))
        except Exception as e:
            print(f"⚠️ Error computing score for {smiles}: {e}")
            continue

    if not scored:
        print("❌ No valid molecules to prioritize.")
        return

    df = pd.DataFrame(scored, columns=["smiles", "mol_wt", "logp", "prioritization_score"])
    df.sort_values("prioritization_score", ascending=False, inplace=True)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df.to_csv(output_path, index=False)
    print(f"✅ Ranked molecules saved to {output_path}")
