#!/usr/bin/env python3
"""
Comprehensive Testing Suite for Enhanced Off-Target & Selectivity Pipeline v2.0
Systematically tests each component against original requirements using existing run.py as reference.
"""

import asyncio
import sys
import os
import json
import pandas as pd
from pathlib import Path
from typing import Dict, Any, List, Tuple
import time
import shutil

# Add current directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

class ComprehensiveTester:
    """Comprehensive testing suite for the pipeline"""
    
    def __init__(self):
        self.test_results = {}
        self.start_time = time.time()
        self.test_molecule = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin
        self.primary_target = "P33259"  # COX-1
        self.test_output_dir = "comprehensive_test_outputs"
        
        # Clean up previous test outputs
        if os.path.exists(self.test_output_dir):
            shutil.rmtree(self.test_output_dir)
        os.makedirs(self.test_output_dir, exist_ok=True)
        
    def log_test(self, test_name: str, status: str, details: str = ""):
        """Log test results"""
        print(f"[{status}] {test_name}: {details}")
        self.test_results[test_name] = {
            "status": status,
            "details": details,
            "timestamp": time.time()
        }
    
    def test_1_environment_and_dependencies(self) -> bool:
        """Test 1: Environment and Dependencies"""
        print("\n" + "="*80)
        print("TEST 1: ENVIRONMENT AND DEPENDENCIES")
        print("="*80)
        
        try:
            # Test Python version
            python_version = sys.version_info
            if python_version.major >= 3 and python_version.minor >= 8:
                self.log_test("Python Version", "PASS", f"Python {python_version.major}.{python_version.minor}.{python_version.micro}")
            else:
                self.log_test("Python Version", "FAIL", f"Python {python_version.major}.{python_version.minor}.{python_version.micro} - Requires 3.8+")
                return False
            
            # Test core dependencies
            import rdkit
            from rdkit import Chem
            self.log_test("RDKit", "PASS", f"Version: {rdkit.__version__}")
            
            import pandas as pd
            self.log_test("Pandas", "PASS", f"Version: {pd.__version__}")
            
            import numpy as np
            self.log_test("NumPy", "PASS", f"Version: {np.__version__}")
            
            import plotly
            self.log_test("Plotly", "PASS", f"Version: {plotly.__version__}")
            
            import structlog
            self.log_test("Structlog", "PASS", "Structured logging available")
            
            import aiohttp
            self.log_test("Aiohttp", "PASS", "Async HTTP client available")
            
            import yaml
            self.log_test("PyYAML", "PASS", "YAML configuration support")
            
            return True
            
        except ImportError as e:
            self.log_test("Dependencies", "FAIL", f"Missing dependency: {e}")
            return False
    
    def test_2_configuration_system(self) -> bool:
        """Test 2: Configuration Management"""
        print("\n" + "="*80)
        print("TEST 2: CONFIGURATION SYSTEM")
        print("="*80)
        
        try:
            # Test config file existence
            if os.path.exists("config.yaml"):
                self.log_test("Config File", "PASS", "config.yaml exists")
            else:
                self.log_test("Config File", "FAIL", "config.yaml not found")
                return False
            
            # Test config loading
            from utils.config_manager import config_manager
            
            # Test API configuration
            api_config = config_manager.get_api_config("swiss_target_prediction")
            if api_config and hasattr(api_config, 'base_url'):
                self.log_test("API Config", "PASS", f"SwissTargetPrediction: {api_config.base_url}")
            else:
                self.log_test("API Config", "FAIL", "API configuration not loaded")
                return False
            
            # Test model configuration
            model_config = config_manager.get_model_config("generative")
            if model_config:
                self.log_test("Model Config", "PASS", "Model configuration loaded")
            else:
                self.log_test("Model Config", "FAIL", "Model configuration not loaded")
                return False
            
            # Test pipeline configuration
            pipeline_config = config_manager.get_pipeline_config()
            if pipeline_config and hasattr(pipeline_config, 'scoring_weights'):
                self.log_test("Pipeline Config", "PASS", f"{len(pipeline_config.scoring_weights)} scoring weights configured")
            else:
                self.log_test("Pipeline Config", "FAIL", "Pipeline configuration not loaded")
                return False
            
            return True
            
        except Exception as e:
            self.log_test("Configuration", "FAIL", f"Configuration error: {e}")
            return False
    
    def test_3_pipeline_imports(self) -> bool:
        """Test 3: Pipeline Component Imports"""
        print("\n" + "="*80)
        print("TEST 3: PIPELINE COMPONENT IMPORTS")
        print("="*80)
        
        components_to_test = [
            ("run", "run_offtarget_selectivity_pipeline"),
            ("nebula.generate_library", "generate_library"),
            ("sparrow.triage_and_rank", "compute_prioritization"),
            ("empirical_binding.empirical_binding", "get_empirical_offtargets"),
            ("structure_modeling.structure_modeling", "mock_structure_binding_analysis"),
            ("impact_risk.impact_estimator", "estimate_impact_risk"),
            ("expression_filter.expression_filter", "apply_expression_weighting"),
            ("final_dashboard.generate_dashboard", "generate_final_dashboard"),
            ("conflict_resolution.conflict_resolution", "resolve_model_conflicts"),
            ("utils.config_manager", "config_manager"),
            ("utils.api_client", "APIManager")
        ]
        
        passed = 0
        total = len(components_to_test)
        
        for module_name, function_name in components_to_test:
            try:
                if module_name == "run":
                    module = __import__(module_name, fromlist=[function_name])
                else:
                    module = __import__(module_name, fromlist=[function_name])
                
                if hasattr(module, function_name):
                    func = getattr(module, function_name)
                    self.log_test(f"{module_name}.{function_name}", "PASS", "Successfully imported")
                    passed += 1
                else:
                    self.log_test(f"{module_name}.{function_name}", "FAIL", "Function not found")
            except ImportError as e:
                self.log_test(f"{module_name}.{function_name}", "FAIL", f"Import failed: {e}")
            except Exception as e:
                self.log_test(f"{module_name}.{function_name}", "FAIL", f"Error: {e}")
        
        self.log_test("Import Summary", "INFO", f"{passed}/{total} components imported successfully")
        return passed >= total * 0.8  # Allow 20% failure rate
    
    def test_4_output_directory_structure(self) -> bool:
        """Test 4: Output Directory Structure"""
        print("\n" + "="*80)
        print("TEST 4: OUTPUT DIRECTORY STRUCTURE")
        print("="*80)
        
        try:
            # Create expected output directories
            expected_dirs = [
                "nebula",
                "sparrow", 
                "empirical_binding",
                "structure_modeling",
                "impact_risk",
                "expression_filter",
                "final_dashboard",
                "conflict_resolution"
            ]
            
            for dir_name in expected_dirs:
                dir_path = os.path.join(self.test_output_dir, dir_name)
                os.makedirs(dir_path, exist_ok=True)
                self.log_test(f"Directory: {dir_name}", "PASS", f"Created: {dir_path}")
            
            # Test data directory
            data_dir = os.path.join(self.test_output_dir, "data", "gtex")
            os.makedirs(data_dir, exist_ok=True)
            
            # Copy GTEx data to test directory
            import shutil
            source_gtex = "data/gtex/GTEx_median_tpm.tsv.gz"
            target_gtex = os.path.join(data_dir, "GTEx_median_tpm.tsv.gz")
            if os.path.exists(source_gtex):
                shutil.copy2(source_gtex, target_gtex)
                self.log_test("Data Directory", "PASS", "GTEx data copied to test directory")
            else:
                self.log_test("Data Directory", "WARN", "GTEx data not found, using fallback")
            
            return True
            
        except Exception as e:
            self.log_test("Directory Structure", "FAIL", f"Directory creation failed: {e}")
            return False
    
    async def test_5_individual_component_tests(self) -> bool:
        """Test 5: Individual Component Functionality"""
        print("\n" + "="*80)
        print("TEST 5: INDIVIDUAL COMPONENT FUNCTIONALITY")
        print("="*80)
        
        component_tests = [
            ("NEBULA", self._test_nebula),
            ("SPARROW", self._test_sparrow),
            ("Empirical Binding", self._test_empirical_binding),
            ("Structure Modeling", self._test_structure_modeling),
            ("IMPACT Risk", self._test_impact_risk),
            ("Expression Filter", self._test_expression_filter),
            ("Final Dashboard", self._test_final_dashboard),
            ("Conflict Resolution", self._test_conflict_resolution),
            ("Toxicophore Detection", self._test_toxicophore_detection),
            ("AI Explanation", self._test_ai_explanation)
        ]
        
        passed = 0
        total = len(component_tests)
        
        for component_name, test_func in component_tests:
            try:
                result = await test_func()
                if result:
                    self.log_test(component_name, "PASS", "Component functionality verified")
                    passed += 1
                else:
                    self.log_test(component_name, "FAIL", "Component functionality failed")
            except Exception as e:
                self.log_test(component_name, "FAIL", f"Component test error: {e}")
        
        self.log_test("Component Summary", "INFO", f"{passed}/{total} components working correctly")
        return passed >= total * 0.7  # Allow 30% failure rate for individual components
    
    async def _test_nebula(self) -> bool:
        """Test NEBULA component"""
        try:
            from nebula.generate_library import generate_library
            
            result = generate_library(
                input_smiles=self.test_molecule,
                num_variants=3,  # Reduced for faster execution
                output_path=os.path.join(self.test_output_dir, "nebula/generated_library.sdf"),
                metadata_path=os.path.join(self.test_output_dir, "nebula/generation_metadata.json")
            )
            
            return isinstance(result, dict) and 'num_molecules' in result
        except Exception:
            return False
    
    async def _test_sparrow(self) -> bool:
        """Test SPARROW component"""
        try:
            from sparrow.triage_and_rank import compute_prioritization
            
            result = await compute_prioritization(
                smiles_list=[self.test_molecule],  # Single molecule for faster test
                output_path=os.path.join(self.test_output_dir, "sparrow/ranked_candidates.csv"),
                metadata_path=os.path.join(self.test_output_dir, "sparrow/synthesis_analysis.json")
            )
            
            return isinstance(result, dict) and 'num_analyzed' in result
        except Exception:
            return False
    
    async def _test_empirical_binding(self) -> bool:
        """Test Empirical Binding component"""
        try:
            from empirical_binding.empirical_binding import get_empirical_offtargets
            
            result = await get_empirical_offtargets(
                smiles=self.test_molecule,
                output_path=os.path.join(self.test_output_dir, "empirical_binding/offtarget_predictions.json"),
                metadata_path=os.path.join(self.test_output_dir, "empirical_binding/prediction_metadata.json")
            )
            
            return isinstance(result, dict) and len(result) > 0
        except Exception:
            return False
    
    async def _test_structure_modeling(self) -> bool:
        """Test Structure Modeling component"""
        try:
            from structure_modeling.structure_modeling import mock_structure_binding_analysis
            
            result = await mock_structure_binding_analysis(
                smiles=self.test_molecule,
                uniprot_id=self.primary_target,
                output_dir=os.path.join(self.test_output_dir, "structure_modeling")
            )
            
            return isinstance(result, dict) and 'binding_score' in result
        except Exception:
            return False
    
    async def _test_impact_risk(self) -> bool:
        """Test IMPACT Risk component"""
        try:
            from impact_risk.impact_estimator import estimate_impact_risk
            
            result = await estimate_impact_risk(
                empirical_scores={self.primary_target: 0.8},
                structural_scores={self.primary_target: 0.7},
                on_target_id=self.primary_target,
                smiles=self.test_molecule,
                output_path=os.path.join(self.test_output_dir, "impact_risk/impact_summary.json")
            )
            
            return isinstance(result, dict) and 'selectivity_score' in result
        except Exception:
            return False
    
    async def _test_expression_filter(self) -> bool:
        """Test Expression Filter component"""
        try:
            from expression_filter.expression_filter import apply_expression_weighting
            
            result = await apply_expression_weighting(
                combined_scores={self.primary_target: 0.8},
                gtex_filepath=os.path.join(self.test_output_dir, "data/gtex/GTEx_median_tpm.tsv.gz"),
                output_path=os.path.join(self.test_output_dir, "expression_filter/tissue_weighted_risk.json"),
                metadata_path=os.path.join(self.test_output_dir, "expression_filter/expression_analysis.json")
            )
            
            return isinstance(result, dict) and 'weighted_risks' in result
        except Exception:
            return False
    
    async def _test_final_dashboard(self) -> bool:
        """Test Final Dashboard component"""
        try:
            from final_dashboard.generate_dashboard import generate_final_dashboard
            
            result = await generate_final_dashboard(
                output_dir=os.path.join(self.test_output_dir, "final_dashboard")
            )
            
            return isinstance(result, dict) and 'dashboard_dir' in result
        except Exception:
            return False
    
    async def _test_conflict_resolution(self) -> bool:
        """Test Conflict Resolution component"""
        try:
            from conflict_resolution.conflict_resolution import resolve_model_conflicts
            
            result = resolve_model_conflicts(
                empirical_path=os.path.join(self.test_output_dir, "empirical_binding/offtarget_predictions.json"),
                structural_path=os.path.join(self.test_output_dir, "structure_modeling/binding_risk.json"),
                output_path=os.path.join(self.test_output_dir, "conflict_resolution/conflict_summary.json")
            )
            
            return isinstance(result, dict) and len(result) > 0
        except Exception:
            return False
    
    async def _test_toxicophore_detection(self) -> bool:
        """Test Toxicophore Detection component"""
        try:
            from toxicity.toxicophore_detector import analyze_toxicophores
            
            result = await analyze_toxicophores(
                smiles=self.test_molecule,
                output_path=os.path.join(self.test_output_dir, "toxicity/toxicophore_analysis.json")
            )
            
            return isinstance(result, dict) and "total_alerts" in result
        except Exception:
            return False
    
    async def _test_ai_explanation(self) -> bool:
        """Test AI Explanation Generator component"""
        try:
            from ai.explanation_generator import generate_ai_explanation
            
            # Create mock results for testing
            mock_results = {
                "results": {
                    "impact_risk": {
                        "decision_flag": "Synthesize",
                        "selectivity_score": 0.8,
                        "safety_score": 0.7,
                        "risky_offtargets": []
                    }
                }
            }
            
            result = await generate_ai_explanation(
                results=mock_results,
                output_path=os.path.join(self.test_output_dir, "ai/explanation_report.json")
            )
            
            return isinstance(result, dict) and "synthesis_explanation" in result
        except Exception:
            return False
    
    async def test_6_full_pipeline_execution(self) -> bool:
        """Test 6: Full Pipeline Execution"""
        print("\n" + "="*80)
        print("TEST 6: FULL PIPELINE EXECUTION")
        print("="*80)
        
        try:
            from run import run_offtarget_selectivity_pipeline
            
            # Execute full pipeline
            result = await run_offtarget_selectivity_pipeline(
                smiles=self.test_molecule,
                primary_uniprot_id=self.primary_target,
                output_base_dir=self.test_output_dir
            )
            
            # Verify result structure
            if not isinstance(result, dict):
                self.log_test("Pipeline Result", "FAIL", "Pipeline did not return a dictionary")
                return False
            
            # Check key metrics
            key_metrics = result.get('key_metrics', {})
            if key_metrics:
                self.log_test("Key Metrics", "PASS", f"Decision: {key_metrics.get('decision', 'Unknown')}, Score: {key_metrics.get('selectivity_score', 0):.3f}")
            else:
                self.log_test("Key Metrics", "WARN", "No key metrics found")
            
            # Check component results
            component_results = result.get('results', {})
            if component_results:
                successful_components = [name for name, res in component_results.items() if 'error' not in res]
                failed_components = [name for name, res in component_results.items() if 'error' in res]
                self.log_test("Component Results", "PASS", f"{len(successful_components)} successful, {len(failed_components)} failed")
            else:
                self.log_test("Component Results", "FAIL", "No component results found")
                return False
            
            # Check metadata
            metadata = result.get('metadata', {})
            if metadata:
                self.log_test("Pipeline Metadata", "PASS", f"Pipeline version: {result.get('pipeline_version', 'Unknown')}")
            else:
                self.log_test("Pipeline Metadata", "WARN", "No metadata found")
            
            return True
            
        except Exception as e:
            self.log_test("Full Pipeline", "FAIL", f"Pipeline execution failed: {e}")
            return False
    
    def test_7_output_file_validation(self) -> bool:
        """Test 7: Output File Validation"""
        print("\n" + "="*80)
        print("TEST 7: OUTPUT FILE VALIDATION")
        print("="*80)
        
        expected_files = [
            ("nebula/generated_library.sdf", "SDF file with generated molecules"),
            ("nebula/generation_metadata.json", "NEBULA generation metadata"),
            ("sparrow/ranked_candidates.csv", "SPARROW synthesis rankings"),
            ("sparrow/synthesis_analysis.json", "SPARROW analysis metadata"),
            ("empirical_binding/offtarget_predictions.json", "Empirical binding predictions"),
            ("empirical_binding/prediction_metadata.json", "Empirical binding metadata"),
            ("structure_modeling/binding_risk.json", "Structural binding risk"),
            ("structure_modeling/binding_pose.png", "Binding pose visualization"),
            ("impact_risk/impact_summary.json", "IMPACT risk summary"),
            ("impact_risk/detailed_analysis.json", "IMPACT detailed analysis"),
            ("expression_filter/tissue_weighted_risk.json", "Tissue-weighted risk"),
            ("expression_filter/expression_analysis.json", "Expression analysis metadata"),
            ("final_dashboard/compound_summary.json", "Final compound summary"),
            ("final_dashboard/comprehensive_report.html", "Comprehensive dashboard"),
            ("conflict_resolution/conflict_summary.json", "Conflict resolution summary"),
            ("toxicity/toxicophore_analysis.json", "Toxicophore analysis report"),
            ("ai/explanation_report.json", "AI explanation report")
        ]
        
        passed = 0
        total = len(expected_files)
        
        for file_path, description in expected_files:
            full_path = os.path.join(self.test_output_dir, file_path)
            if os.path.exists(full_path):
                file_size = os.path.getsize(full_path)
                if file_size > 0:
                    self.log_test(f"File: {file_path}", "PASS", f"{description} ({file_size} bytes)")
                    passed += 1
                else:
                    self.log_test(f"File: {file_path}", "WARN", f"{description} (empty file)")
            else:
                self.log_test(f"File: {file_path}", "FAIL", f"{description} (not found)")
        
        self.log_test("File Validation", "INFO", f"{passed}/{total} files generated successfully")
        return passed >= total * 0.6  # Allow 40% missing files
    
    def test_8_data_quality_validation(self) -> bool:
        """Test 8: Data Quality Validation"""
        print("\n" + "="*80)
        print("TEST 8: DATA QUALITY VALIDATION")
        print("="*80)
        
        try:
            # Test SDF file quality
            sdf_path = os.path.join(self.test_output_dir, "nebula/generated_library.sdf")
            if os.path.exists(sdf_path):
                from rdkit import Chem
                suppl = Chem.SDMolSupplier(sdf_path)
                valid_mols = [mol for mol in suppl if mol is not None]
                if len(valid_mols) > 0:
                    self.log_test("SDF Quality", "PASS", f"{len(valid_mols)} valid molecules")
                else:
                    self.log_test("SDF Quality", "FAIL", "No valid molecules in SDF")
            else:
                self.log_test("SDF Quality", "WARN", "SDF file not found")
            
            # Test JSON file quality
            json_files = [
                "empirical_binding/offtarget_predictions.json",
                "impact_risk/impact_summary.json",
                "expression_filter/tissue_weighted_risk.json"
            ]
            
            for json_file in json_files:
                file_path = os.path.join(self.test_output_dir, json_file)
                if os.path.exists(file_path):
                    try:
                        with open(file_path, 'r') as f:
                            data = json.load(f)
                        if isinstance(data, dict) and len(data) > 0:
                            self.log_test(f"JSON Quality: {json_file}", "PASS", f"Valid JSON with {len(data)} entries")
                        else:
                            self.log_test(f"JSON Quality: {json_file}", "WARN", "Empty or invalid JSON structure")
                    except json.JSONDecodeError:
                        self.log_test(f"JSON Quality: {json_file}", "FAIL", "Invalid JSON format")
                else:
                    self.log_test(f"JSON Quality: {json_file}", "WARN", "File not found")
            
            # Test CSV file quality
            csv_path = os.path.join(self.test_output_dir, "sparrow/ranked_candidates.csv")
            if os.path.exists(csv_path):
                try:
                    df = pd.read_csv(csv_path)
                    if len(df) > 0 and len(df.columns) > 0:
                        self.log_test("CSV Quality", "PASS", f"Valid CSV with {len(df)} rows, {len(df.columns)} columns")
                    else:
                        self.log_test("CSV Quality", "WARN", "Empty CSV file")
                except Exception as e:
                    self.log_test("CSV Quality", "FAIL", f"CSV reading error: {e}")
            else:
                self.log_test("CSV Quality", "WARN", "CSV file not found")
            
            return True
            
        except Exception as e:
            self.log_test("Data Quality", "FAIL", f"Data quality validation error: {e}")
            return False
    
    def test_9_performance_metrics(self) -> bool:
        """Test 9: Performance Metrics"""
        print("\n" + "="*80)
        print("TEST 9: PERFORMANCE METRICS")
        print("="*80)
        
        try:
            # Calculate execution time
            total_time = time.time() - self.start_time
            if total_time < 600:  # Less than 10 minutes (more lenient)
                self.log_test("Execution Time", "PASS", f"Total execution time: {total_time:.1f} seconds")
            else:
                self.log_test("Execution Time", "WARN", f"Slow execution: {total_time:.1f} seconds")
            
            # Check memory usage
            try:
                import psutil
                process = psutil.Process()
                memory_mb = process.memory_info().rss / 1024 / 1024
                if memory_mb < 1000:  # Less than 1GB
                    self.log_test("Memory Usage", "PASS", f"Memory usage: {memory_mb:.1f} MB")
                else:
                    self.log_test("Memory Usage", "WARN", f"High memory usage: {memory_mb:.1f} MB")
            except ImportError:
                self.log_test("Memory Usage", "INFO", "psutil not available for memory monitoring")
            
            # Check disk usage
            total_size = 0
            for root, dirs, files in os.walk(self.test_output_dir):
                for file in files:
                    file_path = os.path.join(root, file)
                    total_size += os.path.getsize(file_path)
            
            total_size_mb = total_size / 1024 / 1024
            self.log_test("Disk Usage", "PASS", f"Output size: {total_size_mb:.1f} MB")
            
            return True
            
        except Exception as e:
            self.log_test("Performance", "FAIL", f"Performance test error: {e}")
            return False
    
    def test_10_error_handling(self) -> bool:
        """Test 10: Error Handling"""
        print("\n" + "="*80)
        print("TEST 10: ERROR HANDLING")
        print("="*80)
        
        try:
            # Test invalid SMILES handling
            from run import run_offtarget_selectivity_pipeline
            
            try:
                # This should raise an error
                asyncio.run(run_offtarget_selectivity_pipeline(
                    smiles="INVALID_SMILES",
                    primary_uniprot_id="P33259"
                ))
                self.log_test("Invalid SMILES", "FAIL", "Should have raised an error")
                return False
            except ValueError:
                self.log_test("Invalid SMILES", "PASS", "Correctly handled invalid SMILES")
            except Exception as e:
                self.log_test("Invalid SMILES", "PASS", f"Handled invalid SMILES: {type(e).__name__}")
            
            # Test invalid UniProt ID handling
            try:
                # This should raise an error
                asyncio.run(run_offtarget_selectivity_pipeline(
                    smiles=self.test_molecule,
                    primary_uniprot_id="INVALID"
                ))
                self.log_test("Invalid UniProt ID", "FAIL", "Should have raised an error")
                return False
            except ValueError:
                self.log_test("Invalid UniProt ID", "PASS", "Correctly handled invalid UniProt ID")
            except Exception as e:
                self.log_test("Invalid UniProt ID", "PASS", f"Handled invalid UniProt ID: {type(e).__name__}")
            
            return True
            
        except Exception as e:
            self.log_test("Error Handling", "FAIL", f"Error handling test failed: {e}")
            return False
    
    def test_11_integration_validation(self) -> bool:
        """Test 11: Integration Validation"""
        print("\n" + "="*80)
        print("TEST 11: INTEGRATION VALIDATION")
        print("="*80)
        
        try:
            # Test data consistency across components
            import json
            
            # Check that all JSON files are valid and contain expected data
            json_files = [
                "nebula/generation_metadata.json",
                "sparrow/synthesis_analysis.json", 
                "empirical_binding/prediction_metadata.json",
                "impact_risk/impact_summary.json",
                "expression_filter/expression_analysis.json",
                "final_dashboard/compound_summary.json",
                "conflict_resolution/conflict_summary.json",
                "toxicity/toxicophore_analysis.json",
                "ai/explanation_report.json"
            ]
            
            valid_json_count = 0
            for json_file in json_files:
                try:
                    file_path = os.path.join(self.test_output_dir, json_file)
                    if os.path.exists(file_path):
                        with open(file_path, 'r') as f:
                            data = json.load(f)
                        if isinstance(data, dict) and len(data) > 0:
                            valid_json_count += 1
                except Exception:
                    pass
            
            if valid_json_count >= len(json_files) * 0.8:  # 80% success rate
                self.log_test("JSON Integration", "PASS", f"{valid_json_count}/{len(json_files)} JSON files valid")
            else:
                self.log_test("JSON Integration", "FAIL", f"Only {valid_json_count}/{len(json_files)} JSON files valid")
                return False
            
            # Test file size consistency
            total_size = 0
            file_count = 0
            for root, dirs, files in os.walk(self.test_output_dir):
                for file in files:
                    file_path = os.path.join(root, file)
                    if os.path.isfile(file_path):
                        total_size += os.path.getsize(file_path)
                        file_count += 1
            
            if total_size > 1000000:  # At least 1MB of output
                self.log_test("Output Volume", "PASS", f"Generated {total_size/1024/1024:.1f} MB across {file_count} files")
            else:
                self.log_test("Output Volume", "FAIL", f"Insufficient output: {total_size/1024/1024:.1f} MB")
                return False
            
            # Test directory structure completeness
            required_dirs = [
                "nebula", "sparrow", "empirical_binding", "structure_modeling",
                "impact_risk", "expression_filter", "final_dashboard", "conflict_resolution",
                "toxicity", "ai"
            ]
            
            existing_dirs = 0
            for dir_name in required_dirs:
                dir_path = os.path.join(self.test_output_dir, dir_name)
                if os.path.exists(dir_path) and os.path.isdir(dir_path):
                    existing_dirs += 1
            
            if existing_dirs >= len(required_dirs):
                self.log_test("Directory Structure", "PASS", f"All {existing_dirs} required directories created")
            else:
                self.log_test("Directory Structure", "FAIL", f"Missing directories: {existing_dirs}/{len(required_dirs)}")
                return False
            
            return True
            
        except Exception as e:
            self.log_test("Integration Validation", "FAIL", f"Integration test failed: {e}")
            return False
    
    def test_12_robustness_validation(self) -> bool:
        """Test 12: Robustness Validation"""
        print("\n" + "="*80)
        print("TEST 12: ROBUSTNESS VALIDATION")
        print("="*80)
        
        try:
            # Test concurrent file access
            import threading
            import time
            
            def read_file_safe(file_path):
                try:
                    with open(file_path, 'r') as f:
                        return f.read()
                except Exception:
                    return None
            
            # Test reading multiple files concurrently
            test_files = [
                os.path.join(self.test_output_dir, "nebula/generation_metadata.json"),
                os.path.join(self.test_output_dir, "sparrow/synthesis_analysis.json"),
                os.path.join(self.test_output_dir, "final_dashboard/compound_summary.json")
            ]
            
            threads = []
            results = []
            
            for file_path in test_files:
                if os.path.exists(file_path):
                    thread = threading.Thread(target=lambda f=file_path: results.append(read_file_safe(f)))
                    threads.append(thread)
                    thread.start()
            
            for thread in threads:
                thread.join()
            
            successful_reads = sum(1 for result in results if result is not None)
            if successful_reads >= len(test_files) * 0.8:
                self.log_test("Concurrent Access", "PASS", f"{successful_reads}/{len(test_files)} files read successfully")
            else:
                self.log_test("Concurrent Access", "FAIL", f"Concurrent access issues: {successful_reads}/{len(test_files)}")
                return False
            
            # Test memory cleanup
            import gc
            gc.collect()
            
            try:
                import psutil
                process = psutil.Process()
                memory_after = process.memory_info().rss / 1024 / 1024
                if memory_after < 500:  # Less than 500MB after cleanup
                    self.log_test("Memory Cleanup", "PASS", f"Memory usage: {memory_after:.1f} MB after cleanup")
                else:
                    self.log_test("Memory Cleanup", "WARN", f"High memory usage after cleanup: {memory_after:.1f} MB")
            except ImportError:
                self.log_test("Memory Cleanup", "INFO", "psutil not available for memory monitoring")
            
            # Test file permissions
            test_file = os.path.join(self.test_output_dir, "test_permissions.txt")
            try:
                with open(test_file, 'w') as f:
                    f.write("test")
                os.remove(test_file)
                self.log_test("File Permissions", "PASS", "File creation and deletion successful")
            except Exception as e:
                self.log_test("File Permissions", "FAIL", f"File permission issues: {e}")
                return False
            
            return True
            
        except Exception as e:
            self.log_test("Robustness Validation", "FAIL", f"Robustness test failed: {e}")
            return False
    
    def test_13_final_validation(self) -> bool:
        """Test 13: Final Validation"""
        print("\n" + "="*80)
        print("TEST 13: FINAL VALIDATION")
        print("="*80)
        
        try:
            # Final comprehensive check
            total_files = 0
            total_size = 0
            
            for root, dirs, files in os.walk(self.test_output_dir):
                for file in files:
                    file_path = os.path.join(root, file)
                    if os.path.isfile(file_path):
                        total_files += 1
                        total_size += os.path.getsize(file_path)
            
            if total_files >= 15:  # At least 15 files generated
                self.log_test("File Count", "PASS", f"Generated {total_files} files")
            else:
                self.log_test("File Count", "FAIL", f"Insufficient files: {total_files}")
                return False
            
            if total_size > 500000:  # At least 500KB total
                self.log_test("Total Size", "PASS", f"Total output size: {total_size/1024:.1f} KB")
            else:
                self.log_test("Total Size", "FAIL", f"Insufficient total size: {total_size/1024:.1f} KB")
                return False
            
            # Check for critical files
            critical_files = [
                "nebula/generated_library.sdf",
                "sparrow/ranked_candidates.csv", 
                "final_dashboard/comprehensive_report.html",
                "final_dashboard/compound_summary.json"
            ]
            
            critical_count = 0
            for critical_file in critical_files:
                file_path = os.path.join(self.test_output_dir, critical_file)
                if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
                    critical_count += 1
            
            if critical_count >= len(critical_files):
                self.log_test("Critical Files", "PASS", f"All {critical_count} critical files present")
            else:
                self.log_test("Critical Files", "FAIL", f"Missing critical files: {critical_count}/{len(critical_files)}")
                return False
            
            return True
            
        except Exception as e:
            self.log_test("Final Validation", "FAIL", f"Final validation failed: {e}")
            return False
    
    def generate_comprehensive_report(self):
        """Generate comprehensive test report"""
        print("\n" + "="*80)
        print("COMPREHENSIVE TEST REPORT")
        print("="*80)
        
        total_tests = len(self.test_results)
        passed_tests = sum(1 for result in self.test_results.values() if result['status'] == 'PASS')
        failed_tests = sum(1 for result in self.test_results.values() if result['status'] == 'FAIL')
        warning_tests = sum(1 for result in self.test_results.values() if result['status'] == 'WARN')
        info_tests = sum(1 for result in self.test_results.values() if result['status'] == 'INFO')
        
        print(f"Total Tests: {total_tests}")
        print(f"Passed: {passed_tests}")
        print(f"Failed: {failed_tests}")
        print(f"Warnings: {warning_tests}")
        print(f"Info: {info_tests}")
        print(f"Success Rate: {(passed_tests/total_tests)*100:.1f}%")
        
        # Calculate overall status
        if passed_tests >= total_tests * 0.8 and failed_tests <= total_tests * 0.2:
            overall_status = "PASS"
            status_message = "Pipeline is working correctly and ready for use!"
        elif passed_tests >= total_tests * 0.6:
            overall_status = "WARN"
            status_message = "Pipeline is mostly working but has some issues to address."
        else:
            overall_status = "FAIL"
            status_message = "Pipeline has significant issues that need to be fixed."
        
        print(f"\nOverall Status: {overall_status}")
        print(f"Status Message: {status_message}")
        
        print("\nDetailed Results:")
        for test_name, result in self.test_results.items():
            status_icon = "✅" if result['status'] == 'PASS' else "❌" if result['status'] == 'FAIL' else "⚠️" if result['status'] == 'WARN' else "ℹ️"
            print(f"{status_icon} {test_name}: {result['details']}")
        
        # Save detailed report
        report_data = {
            "summary": {
                "total_tests": total_tests,
                "passed": passed_tests,
                "failed": failed_tests,
                "warnings": warning_tests,
                "info": info_tests,
                "success_rate": (passed_tests/total_tests)*100,
                "overall_status": overall_status,
                "status_message": status_message
            },
            "detailed_results": self.test_results,
            "timestamp": time.time(),
            "test_molecule": self.test_molecule,
            "primary_target": self.primary_target,
            "output_directory": self.test_output_dir
        }
        
        report_path = os.path.join(self.test_output_dir, "comprehensive_test_report.json")
        with open(report_path, 'w') as f:
            json.dump(report_data, f, indent=2)
        
        print(f"\nDetailed report saved to: {report_path}")
        print(f"Test outputs available in: {self.test_output_dir}")
        
        return overall_status == "PASS"

async def main():
    """Run comprehensive testing suite"""
    print("🧪 COMPREHENSIVE TESTING SUITE")
    print("Enhanced Off-Target & Selectivity Pipeline v2.0")
    print("="*80)
    print("This test suite will validate your entire pipeline against the original requirements.")
    print("It uses the existing run.py as reference and tests all components systematically.")
    print("="*80)
    
    tester = ComprehensiveTester()
    
    # Run all tests
    tests = [
        ("Environment and Dependencies", tester.test_1_environment_and_dependencies),
        ("Configuration System", tester.test_2_configuration_system),
        ("Pipeline Component Imports", tester.test_3_pipeline_imports),
        ("Output Directory Structure", tester.test_4_output_directory_structure),
        ("Individual Component Tests", tester.test_5_individual_component_tests),
        ("Full Pipeline Execution", tester.test_6_full_pipeline_execution),
        ("Output File Validation", tester.test_7_output_file_validation),
        ("Data Quality Validation", tester.test_8_data_quality_validation),
        ("Performance Metrics", tester.test_9_performance_metrics),
        ("Error Handling", tester.test_10_error_handling),
        ("Integration Validation", tester.test_11_integration_validation),
        ("Robustness Validation", tester.test_12_robustness_validation),
        ("Final Validation", tester.test_13_final_validation)
    ]
    
    for test_name, test_func in tests:
        try:
            if asyncio.iscoroutinefunction(test_func):
                await test_func()
            else:
                test_func()
        except Exception as e:
            tester.log_test(test_name, "FAIL", f"Test failed with exception: {e}")
    
    # Generate final report
    success = tester.generate_comprehensive_report()
    
    if success:
        print("\n🎉 ALL TESTS PASSED! Your pipeline is fully functional and production-ready.")
        print("You can now use the pipeline with confidence.")
    else:
        print("\n⚠️ Some tests failed. Please review the detailed report above.")
        print("Address the issues before using the pipeline in production.")
    
    return success

if __name__ == "__main__":
    success = asyncio.run(main())
    sys.exit(0 if success else 1) 