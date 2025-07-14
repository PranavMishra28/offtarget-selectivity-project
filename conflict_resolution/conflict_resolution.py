# offtarget_selectivity/08_conflict_resolution/conflict_resolution.py

import json
import os


def resolve_model_conflicts(empirical_path, structural_path, output_path="conflict_resolution/conflict_summary.json", threshold=0.3):
    '''
    Compares empirical vs. structural scores for each target and flags conflicts or redundancy.
    Accepts structural JSON either as:
      - { "P08172": 0.711, ... }
      - { "uniprot_id": "P08172", "binding_score": 0.711, ... }
    '''
    with open(empirical_path) as f1:
        empirical = json.load(f1)

    with open(structural_path) as f2:
        structural_raw = json.load(f2)

    # Normalize structural input to {uid: score} format
    if isinstance(structural_raw, dict) and "uniprot_id" in structural_raw:
        uid = structural_raw["uniprot_id"]
        score = structural_raw.get("binding_score", 0)
        structural = {uid: score}
    else:
        structural = structural_raw

    resolved = {}
    all_keys = set(empirical) | set(structural)

    for uid in all_keys:
        emp_score = empirical.get(uid, 0)
        struct_score = structural.get(uid, 0)

        # Convert to float to avoid str/int mismatch
        try:
            emp_score = float(emp_score)
            struct_score = float(struct_score)
        except ValueError:
            continue

        delta = abs(emp_score - struct_score)

        if delta >= threshold:
            label = "conflict"
        elif delta < 0.1:
            label = "redundant"
        else:
            label = "agreeing"

        resolved[uid] = {
            "empirical_score": emp_score,
            "structural_score": struct_score,
            "delta": round(delta, 4),
            "relation": label
        }

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(resolved, f, indent=4)

    print(f"✅ Conflict summary saved to {output_path}")
    return resolved
