# offtarget_selectivity/06_expression_filter/expression_filter.py

import os
import json
import random


def get_mock_expression(uniprot_id: str, tissue: str = "brain") -> float:
    """
    Mock expression percentile (0 to 1) for given protein in a tissue.
    Replace with real GTEx or HPA data query.
    """
    random.seed(hash(uniprot_id + tissue) % 100000)
    return round(random.uniform(0.1, 1.0), 2)


def apply_expression_weighting(
    off_target_scores: dict,
    tissues: list = ["brain", "heart", "liver"],
    output_path: str = "expression_filter/tissue_weighted_risk.json"
):
    """
    Multiplies each off-target risk score by its expression percentile in high-risk tissues.

    Args:
        off_target_scores (dict): {uniprot_id: risk_score}
        tissues (list): Tissues to consider
        output_path (str): Save final tissue-weighted risk
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    weighted_scores = {}
    for uid, base_score in off_target_scores.items():
        # Get average expression across specified tissues
        expr_scores = [get_mock_expression(uid, tissue) for tissue in tissues]
        avg_expr = round(sum(expr_scores) / len(expr_scores), 3)
        weighted_score = round(base_score * avg_expr, 4)

        weighted_scores[uid] = {
            "original_score": base_score,
            "avg_expression_percentile": avg_expr,
            "tissue_weighted_score": weighted_score
        }

    with open(output_path, "w") as f:
        json.dump(weighted_scores, f, indent=4)

    print(f"✅ Tissue-weighted scores saved to {output_path}")
    return weighted_scores
