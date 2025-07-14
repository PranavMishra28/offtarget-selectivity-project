# offtarget_selectivity/final_dashboard/generate_dashboard.py

import os
import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def generate_final_dashboard(
    impact_path="impact_risk/impact_summary.json",
    expression_path="expression_filter/tissue_weighted_risk.json",
    output_dir="final_dashboard"
):
    os.makedirs(output_dir, exist_ok=True)

    # Load inputs
    with open(impact_path, "r") as f:
        impact = json.load(f)
    with open(expression_path, "r") as f:
        expression = json.load(f)

    # Combine into summary
    summary = {
        "compound_on_target": impact["on_target_id"],
        "on_target_score": impact["on_target_score"],
        "avg_off_target_score": impact["avg_off_target_score"],
        "selectivity_score": impact["selectivity_score"],
        "decision_flag": impact["decision_flag"],
        "top_risky_offtargets": impact["risky_offtargets"],
        "expression_weighted_risks": expression
    }

    # Save full JSON
    json_path = os.path.join(output_dir, "compound_summary.json")
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=4)
    print(f"✅ Saved summary to {json_path}")

    # Save CSV (flatten expression risks)
    flat_data = []
    for uid, info in expression.items():
        flat_data.append({
            "uniprot_id": uid,
            "original_score": info["original_score"],
            "expression_weighted_score": info["tissue_weighted_score"],
            "avg_expression": info.get("avg_expression_tpm") or info.get("avg_expression_percentile", "N/A")
        })
    df = pd.DataFrame(flat_data)
    df.to_csv(os.path.join(output_dir, "compound_summary.csv"), index=False)
    print(f"✅ Saved summary table to compound_summary.csv")

    # Optional Plot: Top 10 expression-weighted risks
    df_sorted = df.sort_values("expression_weighted_score", ascending=False).head(10)
    plt.figure(figsize=(10, 5))
    sns.barplot(data=df_sorted, x="uniprot_id", y="expression_weighted_score", palette="rocket")
    plt.title("Top 10 Tissue-Weighted Off-Target Risks")
    plt.ylabel("Tissue-Weighted Risk Score")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt_path = os.path.join(output_dir, "selectivity_dashboard.png")
    plt.savefig(plt_path)
    plt.close()
    print(f"✅ Dashboard image saved to {plt_path}")
