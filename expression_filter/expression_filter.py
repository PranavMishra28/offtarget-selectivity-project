# offtarget_selectivity/06_expression_filter/expression_filter.py

import pandas as pd
import os
import json
from scipy.stats import percentileofscore

def apply_expression_weighting(
    combined_scores: dict,
    gtex_filepath="data/gtex/GTEx_median_tpm.tsv.gz",
    output_path="expression_filter/tissue_weighted_risk.json"
):
    '''
    Reweights risk scores based on average tissue expression from GTEx data.
    GTEx data expected in TSV format with genes as index and tissues as columns.
    '''
    print(f"📊 Loading GTEx data from {gtex_filepath}...")
    gtex = pd.read_csv(gtex_filepath, sep='\\t', index_col=0)

    # Calculate average expression across tissues
    gtex['avg_tpm'] = gtex.iloc[:, 1:].mean(axis=1)

    results = {}

    all_expr_values = gtex['avg_tpm'].values
    for uid, score in combined_scores.items():
        expr = gtex.loc[uid, 'avg_tpm'] if uid in gtex.index else 0.0
        expr_weight = min(expr / 100.0, 1.0)  # normalize and cap
        weighted_score = round(score * (1 + expr_weight), 4)
        percentile = round(percentileofscore(all_expr_values, expr), 2)

        results[uid] = {
            "original_score": score,
            "avg_expression_tpm": round(expr, 2),
            "avg_expression_percentile": percentile,
            "tissue_weighted_score": weighted_score
        }

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(results, f, indent=4)

    print(f"✅ Saved tissue-weighted risk scores to {output_path}")
    return results
