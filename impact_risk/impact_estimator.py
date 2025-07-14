# offtarget_selectivity/05_impact_risk/impact_estimator.py

import os
import json


def estimate_impact_risk(
    empirical_scores: dict,
    structural_scores: dict,
    on_target_id: str,
    output_path="impact_risk/impact_summary.json",
):
    """
    Combines structural and empirical scores to estimate selectivity and risk.

    Returns:
        A JSON file containing:
            - selectivity_score
            - risky_offtargets
            - decision_flag
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Merge structural + empirical scores by averaging
    combined_scores = {}
    all_uids = set(empirical_scores) | set(structural_scores)
    for uid in all_uids:
        emp = empirical_scores.get(uid, 0.0)
        struct = structural_scores.get(uid, 0.0)
        combined_scores[uid] = round((emp + struct) / 2, 4)

    # Check if primary target is present
    on_score = combined_scores.get(on_target_id)
    if on_score is None:
        print(f"⚠️ Warning: Primary target {on_target_id} not found in combined scores.")

        summary = {
            "on_target_id": on_target_id,
            "on_target_score": None,
            "avg_off_target_score": None,
            "selectivity_score": None,
            "risky_offtargets": [],
            "decision_flag": "Reject",
            "note": "Primary target not found in predictions.",
        }

        with open(output_path, "w") as f:
            json.dump(summary, f, indent=4)

        return summary

    # Calculate average off-target score
    off_scores = [v for k, v in combined_scores.items() if k != on_target_id]
    avg_off_target = round(sum(off_scores) / len(off_scores), 4) if off_scores else None

    # Selectivity index
    selectivity_score = (
        round(on_score / avg_off_target, 4)
        if avg_off_target and avg_off_target != 0
        else None
    )

    # Risk ranking
    risky_offtargets = sorted(
        [(k, v) for k, v in combined_scores.items() if k != on_target_id],
        key=lambda x: x[1],
        reverse=True,
    )[:10]

    # Decision logic
    if selectivity_score and selectivity_score > 2:
        decision = "Synthesize"
    elif selectivity_score and 1 <= selectivity_score <= 2:
        decision = "Watch"
    else:
        decision = "Reject"

    summary = {
        "on_target_id": on_target_id,
        "on_target_score": on_score,
        "avg_off_target_score": avg_off_target,
        "selectivity_score": selectivity_score,
        "risky_offtargets": [
            {"uniprot_id": uid, "combined_score": score}
            for uid, score in risky_offtargets
        ],
        "decision_flag": decision,
    }

    with open(output_path, "w") as f:
        json.dump(summary, f, indent=4)

    print(f"✅ Impact summary saved to {output_path}")
    return summary
