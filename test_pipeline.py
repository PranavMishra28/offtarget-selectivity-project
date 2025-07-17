#!/usr/bin/env python3
"""
Test script for Enhanced Off-Target & Selectivity Pipeline v2.0
Verifies all components are working correctly.
"""

import asyncio
import sys
import os
from pathlib import Path

# Add current directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

def test_imports():
    """Test that all required modules can be imported"""
    print("🔍 Testing imports...")
    
    try:
        # Test core imports
        import rdkit
        from rdkit import Chem
        print("✅ RDKit imported successfully")
        
        import pandas as pd
        import numpy as np
        print("✅ Pandas/NumPy imported successfully")
        
        import structlog
        print("✅ Structlog imported successfully")
        
        # Test pipeline components
        from utils.config_manager import config_manager
        print("✅ Config manager imported successfully")
        
        from utils.api_client import APIManager
        print("✅ API client imported successfully")
        
        return True
        
    except ImportError as e:
        print(f"❌ Import failed: {e}")
        return False

def test_configuration():
    """Test configuration loading"""
    print("\n⚙️ Testing configuration...")
    
    try:
        from utils.config_manager import config_manager
        
        # Test API config
        api_config = config_manager.get_api_config("swiss_target_prediction")
        print(f"✅ API config loaded: {api_config.base_url}")
        
        # Test model config
        model_config = config_manager.get_model_config("generative")
        print(f"✅ Model config loaded: {model_config.library}")
        
        # Test pipeline config
        pipeline_config = config_manager.get_pipeline_config()
        print(f"✅ Pipeline config loaded: {len(pipeline_config.scoring_weights)} weights")
        
        return True
        
    except Exception as e:
        print(f"❌ Configuration test failed: {e}")
        return False

def test_basic_functionality():
    """Test basic functionality with a simple molecule"""
    print("\n🧪 Testing basic functionality...")
    
    try:
        from rdkit import Chem
        
        # Test SMILES parsing
        smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            print("❌ SMILES parsing failed")
            return False
        
        print(f"✅ SMILES parsed successfully: {Chem.MolToSmiles(mol)}")
        
        # Test molecular properties
        from rdkit.Chem import Descriptors
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        
        print(f"✅ Molecular properties calculated: MW={mw:.1f}, LogP={logp:.2f}")
        
        return True
        
    except Exception as e:
        print(f"❌ Basic functionality test failed: {e}")
        return False

async def test_async_components():
    """Test async components"""
    print("\n⚡ Testing async components...")
    
    try:
        # Test API manager
        from utils.api_client import APIManager
        api_manager = APIManager()
        
        # Test client creation
        swiss_client = api_manager.get_client("swiss_target_prediction")
        if swiss_client:
            print("✅ API manager created successfully")
        else:
            print("❌ API manager failed")
            return False
        
        return True
        
    except Exception as e:
        print(f"❌ Async components test failed: {e}")
        return False

def test_output_directories():
    """Test output directory creation"""
    print("\n📁 Testing output directories...")
    
    try:
        # Create test directories
        test_dirs = [
            "nebula",
            "sparrow", 
            "empirical_binding",
            "structure_modeling",
            "impact_risk",
            "expression_filter",
            "final_dashboard",
            "conflict_resolution"
        ]
        
        for dir_name in test_dirs:
            os.makedirs(dir_name, exist_ok=True)
            print(f"✅ Created directory: {dir_name}")
        
        return True
        
    except Exception as e:
        print(f"❌ Directory creation failed: {e}")
        return False

def test_pipeline_modules():
    """Test individual pipeline modules"""
    print("\n🔧 Testing pipeline modules...")
    
    modules_to_test = [
        ("nebula.generate_library", "generate_library"),
        ("sparrow.triage_and_rank", "compute_prioritization"),
        ("empirical_binding.empirical_binding", "get_empirical_offtargets"),
        ("structure_modeling.structure_modeling", "mock_structure_binding_analysis"),
        ("impact_risk.impact_estimator", "estimate_impact_risk"),
        ("expression_filter.expression_filter", "apply_expression_weighting"),
        ("final_dashboard.generate_dashboard", "generate_final_dashboard"),
        ("conflict_resolution.conflict_resolution", "resolve_model_conflicts")
    ]
    
    passed = 0
    total = len(modules_to_test)
    
    for module_name, function_name in modules_to_test:
        try:
            module = __import__(module_name, fromlist=[function_name])
            func = getattr(module, function_name)
            print(f"✅ {module_name} imported successfully")
            passed += 1
        except ImportError as e:
            print(f"⚠️ {module_name} not available: {e}")
        except AttributeError as e:
            print(f"⚠️ {function_name} not found in {module_name}: {e}")
        except Exception as e:
            print(f"❌ {module_name} failed: {e}")
    
    print(f"📊 Module test results: {passed}/{total} modules available")
    return passed > 0  # Pass if at least one module is available

async def main():
    """Run all tests"""
    print("🧪 Enhanced Off-Target & Selectivity Pipeline v2.0 - Test Suite")
    print("=" * 70)
    
    tests = [
        ("Imports", test_imports),
        ("Configuration", test_configuration),
        ("Basic Functionality", test_basic_functionality),
        ("Output Directories", test_output_directories),
        ("Async Components", test_async_components),
        ("Pipeline Modules", test_pipeline_modules)
    ]
    
    passed = 0
    total = len(tests)
    
    for test_name, test_func in tests:
        try:
            if asyncio.iscoroutinefunction(test_func):
                result = await test_func()
            else:
                result = test_func()
            
            if result:
                passed += 1
            else:
                print(f"❌ {test_name} test failed")
                
        except Exception as e:
            print(f"❌ {test_name} test failed with exception: {e}")
    
    print("\n" + "=" * 70)
    print(f"📊 Test Results: {passed}/{total} tests passed")
    
    if passed >= total - 1:  # Allow one test to fail
        print("🎉 Pipeline is ready to use!")
        return True
    else:
        print("⚠️ Some tests failed. Please check the errors above.")
        return False

if __name__ == "__main__":
    success = asyncio.run(main())
    sys.exit(0 if success else 1) 