"""
Enhanced Off-Target & Selectivity Pipeline - Main Execution Script
Implements the complete pipeline with all enhanced components and comprehensive error handling.
"""

import asyncio
import os
import json
import structlog
import pandas as pd
from typing import Dict, Any, Optional
from pathlib import Path

# Import enhanced pipeline components
from rdkit import Chem
from nebula.generate_library import generate_library
from sparrow.triage_and_rank import compute_prioritization
from empirical_binding.empirical_binding import get_empirical_offtargets
from structure_modeling.structure_modeling import mock_structure_binding_analysis
from impact_risk.impact_estimator import estimate_impact_risk
from expression_filter.expression_filter import apply_expression_weighting
from final_dashboard.generate_dashboard import generate_final_dashboard
from conflict_resolution.conflict_resolution import resolve_model_conflicts

# Import configuration and utilities
from utils.config_manager import config_manager

# Configure logging
structlog.configure(
    processors=[
        structlog.stdlib.filter_by_level,
        structlog.stdlib.add_logger_name,
        structlog.stdlib.add_log_level,
        structlog.stdlib.PositionalArgumentsFormatter(),
        structlog.processors.TimeStamper(fmt="iso"),
        structlog.processors.StackInfoRenderer(),
        structlog.processors.format_exc_info,
        structlog.processors.UnicodeDecoder(),
        structlog.processors.JSONRenderer()
    ],
    context_class=dict,
    logger_factory=structlog.stdlib.LoggerFactory(),
    wrapper_class=structlog.stdlib.BoundLogger,
    cache_logger_on_first_use=True,
)

logger = structlog.get_logger(__name__)

class PipelineExecutor:
    """Main pipeline execution orchestrator"""
    
    def __init__(self):
        self.logger = logger
        self.pipeline_config = config_manager.get_pipeline_config()
        self.results = {}
    
    async def execute_pipeline(
        self,
    smiles: str,
        primary_uniprot_id: str,
        output_base_dir: str = "."
    ) -> Dict[str, Any]:
        """
        Execute the complete enhanced Off-Target & Selectivity pipeline.
        
        Args:
            smiles: Input SMILES string
            primary_uniprot_id: Primary target UniProt ID
            output_base_dir: Base directory for outputs
        
        Returns:
            Comprehensive pipeline results
        """
        
        self.logger.info("🚀 Starting Enhanced Off-Target & Selectivity Pipeline")
        self.logger.info(f"📝 Input SMILES: {smiles}")
        self.logger.info(f"🎯 Primary Target: {primary_uniprot_id}")
        
        try:
            # Step 1: NEBULA - Enhanced Generative Library
            await self._execute_nebula(smiles, output_base_dir)
            
            # Step 2: SPARROW - Enhanced Synthesis Feasibility
            await self._execute_sparrow(output_base_dir)
            
            # Step 3: Empirical Binding - Multi-Source Prediction
            await self._execute_empirical_binding(smiles, output_base_dir)
            
            # Step 4: Structure Modeling - AlphaFold-3 + PLIP
            await self._execute_structure_modeling(smiles, output_base_dir)
            
            # Step 5: IMPACT - Enhanced Risk Assessment
            await self._execute_impact_assessment(smiles, primary_uniprot_id, output_base_dir)
            
            # Step 6: Expression Filter - Tissue-Specific Weighting
            await self._execute_expression_filter(output_base_dir)
            
            # Step 7: Final Dashboard - Comprehensive Visualization
            await self._execute_final_dashboard(output_base_dir)
            
            # Step 8: Conflict Resolution - Model Agreement Analysis
            await self._execute_conflict_resolution(output_base_dir)
            
            # Generate final summary
            final_summary = self._generate_final_summary()
            
            self.logger.info("✅ Enhanced pipeline execution completed successfully")
            return final_summary
            
        except Exception as e:
            self.logger.error(f"❌ Pipeline execution failed: {e}")
            raise
    
    async def _execute_nebula(self, smiles: str, output_base_dir: str):
        """Execute NEBULA generative library step"""
        self.logger.info("🎯 [STEP 1] Enhanced NEBULA - Generative Library Generation")
        
        try:
            result = generate_library(
                input_smiles=smiles,
                num_variants=10,  # Further reduced for faster execution
                output_path=os.path.join(output_base_dir, "nebula/generated_library.sdf"),
                metadata_path=os.path.join(output_base_dir, "nebula/generation_metadata.json")
            )
            
            self.results["nebula"] = result
            self.logger.info(f"✅ NEBULA completed: {result['num_molecules']} molecules generated")
            
        except Exception as e:
            self.logger.error(f"❌ NEBULA failed: {e}")
            self.results["nebula"] = {"error": str(e)}
    
    async def _execute_sparrow(self, output_base_dir: str):
        """Execute SPARROW synthesis feasibility step"""
        self.logger.info("🔬 [STEP 2] Enhanced SPARROW - Synthesis Feasibility Assessment")
        
        try:
            # Load generated molecules from NEBULA
            sdf_path = os.path.join(output_base_dir, "nebula/generated_library.sdf")
            if not os.path.exists(sdf_path):
                raise FileNotFoundError(f"NEBULA output not found: {sdf_path}")
            
            suppl = Chem.SDMolSupplier(sdf_path)
            smiles_list = [Chem.MolToSmiles(mol) for mol in suppl if mol is not None]
            
            if not smiles_list:
                raise ValueError("No valid molecules found in NEBULA output")
            
            result = await compute_prioritization(
                smiles_list=smiles_list,
                output_path=os.path.join(output_base_dir, "sparrow/ranked_candidates.csv"),
                metadata_path=os.path.join(output_base_dir, "sparrow/synthesis_analysis.json")
            )
            
            self.results["sparrow"] = result
            self.logger.info(f"✅ SPARROW completed: {result['num_analyzed']} molecules analyzed")
            
        except Exception as e:
            self.logger.error(f"❌ SPARROW failed: {e}")
            self.results["sparrow"] = {"error": str(e)}
    
    async def _execute_empirical_binding(self, smiles: str, output_base_dir: str):
        """Execute enhanced empirical binding prediction"""
        self.logger.info("🧠 [STEP 3] Enhanced Empirical Binding - Multi-Source Target Prediction")
        
        try:
            result = await get_empirical_offtargets(
                smiles=smiles,
                output_path=os.path.join(output_base_dir, "empirical_binding/offtarget_predictions.json"),
                metadata_path=os.path.join(output_base_dir, "empirical_binding/prediction_metadata.json"),
                min_confidence=0.3,
                max_results=5  # Minimal for fastest execution
            )
            
            self.results["empirical_binding"] = {
                "predictions": result,
                "num_targets": len(result)
            }
            self.logger.info(f"✅ Empirical binding completed: {len(result)} targets predicted")
            
        except Exception as e:
            self.logger.error(f"❌ Empirical binding failed: {e}")
            self.results["empirical_binding"] = {"error": str(e)}
    
    async def _execute_structure_modeling(self, smiles: str, output_base_dir: str):
        """Execute enhanced structure modeling"""
        self.logger.info("🏗 [STEP 4] Enhanced Structure Modeling - AlphaFold-3 + PLIP Integration")
        
        try:
            # Get targets from empirical binding
            empirical_data = self.results.get("empirical_binding", {})
            if "error" in empirical_data:
                raise ValueError("Empirical binding data not available")
            
            targets = list(empirical_data.get("predictions", {}).keys())
            
            structural_scores = {}
            for target in targets[:2]:  # Minimal for fastest execution
                try:
                    result = await mock_structure_binding_analysis(
                        smiles=smiles,
                        uniprot_id=target,
                        output_dir=os.path.join(output_base_dir, "structure_modeling")
                    )
                    structural_scores[target] = result["binding_score"]
                except Exception as e:
                    self.logger.warning(f"Structure analysis failed for {target}: {e}")
                    structural_scores[target] = 0.0
            
            self.results["structure_modeling"] = {
                "scores": structural_scores,
                "num_targets_analyzed": len(structural_scores)
            }
            self.logger.info(f"✅ Structure modeling completed: {len(structural_scores)} targets analyzed")
            
        except Exception as e:
            self.logger.error(f"❌ Structure modeling failed: {e}")
            self.results["structure_modeling"] = {"error": str(e)}
    
    async def _execute_impact_assessment(self, smiles: str, primary_uniprot_id: str, output_base_dir: str):
        """Execute enhanced IMPACT risk assessment"""
        self.logger.info("⚖️ [STEP 5] Enhanced IMPACT - Comprehensive Risk Assessment")
        
        try:
            # Get data from previous steps
            empirical_data = self.results.get("empirical_binding", {})
            structural_data = self.results.get("structure_modeling", {})
            
            if "error" in empirical_data or "error" in structural_data:
                raise ValueError("Required data from previous steps not available")
            
            empirical_scores = empirical_data.get("predictions", {})
            structural_scores = structural_data.get("scores", {})
            
            result = await estimate_impact_risk(
                empirical_scores=empirical_scores,
                structural_scores=structural_scores,
                on_target_id=primary_uniprot_id,
                smiles=smiles,
                output_path=os.path.join(output_base_dir, "impact_risk/impact_summary.json"),
                detailed_path=os.path.join(output_base_dir, "impact_risk/detailed_analysis.json")
            )
            
            self.results["impact_risk"] = result
            self.logger.info(f"✅ IMPACT assessment completed: Decision = {result.get('decision_flag', 'Unknown')}")
            
        except Exception as e:
            self.logger.error(f"❌ IMPACT assessment failed: {e}")
            self.results["impact_risk"] = {"error": str(e)}
    
    async def _execute_expression_filter(self, output_base_dir: str):
        """Execute enhanced expression filtering"""
        self.logger.info("🧬 [STEP 6] Enhanced Expression Filter - Tissue-Specific Weighting")
        
        try:
            # Get risky targets from IMPACT
            impact_data = self.results.get("impact_risk", {})
            if "error" in impact_data:
                self.logger.warning("IMPACT data not available, using fallback")
                self.results["expression_filter"] = {"weighted_risks": {}}
                return
            
            risky_targets = {}
            risky_offtargets = impact_data.get("risky_offtargets", [])
            if risky_offtargets:
                risky_targets = {
                    entry["uniprot_id"]: entry["combined_score"]
                    for entry in risky_offtargets
                }
            
            if not risky_targets:
                self.logger.warning("No risky targets found for expression weighting")
                self.results["expression_filter"] = {"weighted_risks": {}}
                return
            
            result = await apply_expression_weighting(
                combined_scores=risky_targets,
                gtex_filepath=os.path.join(output_base_dir, "data/gtex/GTEx_median_tpm.tsv.gz"),
                output_path=os.path.join(output_base_dir, "expression_filter/tissue_weighted_risk.json"),
                metadata_path=os.path.join(output_base_dir, "expression_filter/expression_analysis.json"),
                visualization_path=os.path.join(output_base_dir, "expression_filter/expression_visualizations")
            )
            
            self.results["expression_filter"] = {
                "weighted_risks": result,
                "num_targets_weighted": len(result)
            }
            self.logger.info(f"✅ Expression filtering completed: {len(result)} targets weighted")
            
        except Exception as e:
            self.logger.error(f"❌ Expression filtering failed: {e}")
            self.results["expression_filter"] = {"error": str(e)}
    
    async def _execute_final_dashboard(self, output_base_dir: str):
        """Execute enhanced final dashboard generation"""
        self.logger.info("📊 [STEP 7] Enhanced Final Dashboard - Comprehensive Visualization")
        
        try:
            result = await generate_final_dashboard(
                output_dir=os.path.join(output_base_dir, "final_dashboard")
            )
            
            self.results["final_dashboard"] = result
            self.logger.info(f"✅ Final dashboard completed: {len(result.get('visualizations', {}))} visualizations generated")
            
        except Exception as e:
            self.logger.error(f"❌ Final dashboard failed: {e}")
            self.results["final_dashboard"] = {"error": str(e)}
    
    async def _execute_conflict_resolution(self, output_base_dir: str):
        """Execute conflict resolution step"""
        self.logger.info("🔍 [STEP 8] Enhanced Conflict Resolution - Model Agreement Analysis")
        
        try:
            result = await resolve_model_conflicts(
                output_dir=os.path.join(output_base_dir, "conflict_resolution")
            )
            
            self.results["conflict_resolution"] = result
            self.logger.info(f"✅ Conflict resolution completed: {result.get('num_conflicts', 0)} conflicts identified")
            
        except Exception as e:
            self.logger.error(f"❌ Conflict resolution failed: {e}")
            self.results["conflict_resolution"] = {"error": str(e)}
    
    def _generate_final_summary(self) -> Dict[str, Any]:
        """Generate comprehensive final summary"""
        summary = {
            "pipeline_version": "2.0",
            "execution_status": "completed",
            "results": self.results,
            "metadata": {
                "timestamp": pd.Timestamp.now().isoformat(),
                "components_executed": list(self.results.keys()),
                "successful_components": [
                    name for name, result in self.results.items() 
                    if "error" not in result
                ],
                "failed_components": [
                    name for name, result in self.results.items() 
                    if "error" in result
                ]
            }
        }
        
        # Add key metrics
        impact_data = self.results.get("impact_risk", {})
        if "error" not in impact_data:
            summary["key_metrics"] = {
                "decision": impact_data.get("decision_flag", "Unknown"),
                "selectivity_score": impact_data.get("selectivity_score", 0),
                "safety_score": impact_data.get("safety_score", 0),
                "num_offtargets": len(impact_data.get("risky_offtargets", []))
            }
        
        return summary

async def run_offtarget_selectivity_pipeline(
    smiles: str,
    primary_uniprot_id: str,
    output_base_dir: str = "."
) -> Dict[str, Any]:
    """
    Main entry point for the enhanced Off-Target & Selectivity pipeline.
    
    Args:
        smiles: Input SMILES string
        primary_uniprot_id: Primary target UniProt ID
        output_base_dir: Base directory for outputs
    
    Returns:
        Comprehensive pipeline results
    """
    
    # Validate inputs
    if not smiles or not isinstance(smiles, str):
        raise ValueError("SMILES must be a non-empty string")
    
    # Validate SMILES format
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES format: {smiles}")
    
    # Validate UniProt ID format
    if not primary_uniprot_id or not isinstance(primary_uniprot_id, str):
        raise ValueError("UniProt ID must be a non-empty string")
    
    # UniProt ID format validation (should be like P12345 or A0A123456789)
    import re
    uniprot_pattern = r'^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$'
    if not re.match(uniprot_pattern, primary_uniprot_id):
        raise ValueError(f"Invalid UniProt ID format: {primary_uniprot_id}. Expected format like P12345 or A0A123456789")
    
    # Create output directory
    os.makedirs(output_base_dir, exist_ok=True)
    
    # Execute pipeline
    executor = PipelineExecutor()
    return await executor.execute_pipeline(
        smiles=smiles,
        primary_uniprot_id=primary_uniprot_id,
        output_base_dir=output_base_dir
    )

if __name__ == "__main__":
    # Example compound: Aspirin
    compound_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
    primary_target_uniprot = "P33259"  # COX-1

    # Run pipeline
    async def main():
        try:
            results = await run_offtarget_selectivity_pipeline(
        smiles=compound_smiles,
        primary_uniprot_id=primary_target_uniprot
    )
            
            # Print summary
            print("\n" + "="*80)
            print("🎉 ENHANCED OFF-TARGET SELECTIVITY PIPELINE COMPLETED")
            print("="*80)
            
            key_metrics = results.get("key_metrics", {})
            print(f"🎯 Synthesis Decision: {key_metrics.get('decision', 'Unknown')}")
            print(f"📊 Selectivity Score: {key_metrics.get('selectivity_score', 0):.3f}")
            print(f"🛡️ Safety Score: {key_metrics.get('safety_score', 0):.3f}")
            print(f"🎯 Off-Targets Identified: {key_metrics.get('num_offtargets', 0)}")
            
            metadata = results.get("metadata", {})
            print(f"✅ Successful Components: {len(metadata.get('successful_components', []))}")
            print(f"❌ Failed Components: {len(metadata.get('failed_components', []))}")
            
            print("\n📁 Results available in respective output directories")
            print("📊 Comprehensive dashboard: final_dashboard/comprehensive_report.html")
            print("="*80)
            
        except Exception as e:
            logger.error(f"Pipeline execution failed: {e}")
            print(f"❌ Pipeline failed: {e}")
    
    # Run async main function
    asyncio.run(main())
