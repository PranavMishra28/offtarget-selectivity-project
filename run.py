# run.py

from rdkit import Chem
from sparrow.triage_and_rank import compute_prioritization
from nebula.generate_library import generate_library
from empirical_binding.empirical_binding import get_empirical_offtargets
from structure_modeling.structure_modeling import mock_structure_binding_analysis
from impact_risk.impact_estimator import estimate_impact_risk
from expression_filter.expression_filter import apply_expression_weighting
from final_dashboard.generate_dashboard import generate_final_dashboard
from conflict_resolution.conflict_resolution import resolve_model_conflicts



def run_offtarget_selectivity_pipeline(
    smiles: str,
    primary_uniprot_id: str
):
    print("\n🎯 [STEP 1] Generating 3D Molecule Library (NEBULA)")
    generate_library(smiles)

    print("\n🔬 [STEP 2] Prioritizing Candidates by Synthesis Difficulty (SPARROW)")
    sdf_path = "nebula/generated_library.sdf"
    suppl = Chem.SDMolSupplier(sdf_path)
    smiles_list = [Chem.MolToSmiles(mol) for mol in suppl if mol is not None]
    compute_prioritization(smiles_list)

    print("\n🧠 [STEP 3] Predicting Empirical Off-Targets (SwissTargetPrediction)")
    empirical_scores = get_empirical_offtargets(smiles)

    print("\n🏗 [STEP 4] Predicting Structural Binding Risks (Mock AlphaFold + PLIP)")
    structural_scores = {}
    for uid in empirical_scores.keys():
        result = mock_structure_binding_analysis(smiles, uid)
        structural_scores[uid] = result["binding_score"]

    print("\n⚖️ [STEP 5] Estimating Selectivity & Safety (IMPACT)")
    impact_summary = estimate_impact_risk(empirical_scores, structural_scores, primary_uniprot_id)

    risky_targets = {
        entry["uniprot_id"]: entry["combined_score"]
        for entry in impact_summary["risky_offtargets"]
    }

    print("\n🧬 [STEP 6] Reweighting Risks by Tissue Expression (GTEx/HPA)")
    apply_expression_weighting(risky_targets)

    print("\n📊 [STEP 7] Generating Final Summary Dashboard")
    generate_final_dashboard()

    print("\n⚕️ [STEP 8] Optional Task: Conflict & Redundancy Resolution")
    conflict_summary = resolve_model_conflicts(
    empirical_path="empirical_binding/offtarget_predictions.json",
    structural_path="structure_modeling/binding_risk.json"
)
    if conflict_summary:
        print(conflict_summary)
    print("\n✅ All steps complete. Check respective folders for results.")


if __name__ == "__main__":
    # Example compound: Aspirin
    compound_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
    primary_target_uniprot = "P33259"  # COX-1

    run_offtarget_selectivity_pipeline(
        smiles=compound_smiles,
        primary_uniprot_id=primary_target_uniprot
    )
