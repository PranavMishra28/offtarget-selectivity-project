"""
Enhanced IMPACT - Advanced Selectivity & Safety Assessment
Implements comprehensive selectivity scoring, ion channel liability assessment, and safety profiling.
"""

import os
import json
import asyncio
from typing import Dict, Any, List, Optional, Tuple
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import structlog
from tqdm import tqdm

# Import configuration
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.config_manager import config_manager

class IonChannelPredictor:
    """Ion channel liability prediction (hERG, Nav, Ca2+)"""
    
    def __init__(self):
        self.logger = structlog.get_logger(__name__)
        self.risk_threshold = 0.3  # Configurable threshold
    
    def predict_ion_channel_risk(self, smiles: str) -> Dict[str, Any]:
        """Predict ion channel liability risks"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return self._empty_ion_channel_results()
        
        # Calculate molecular descriptors for ion channel prediction
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        aromatic_rings = Descriptors.NumAromaticRings(mol)
        rotatable = Descriptors.NumRotatableBonds(mol)
        
        # Predict hERG liability
        herg_risk = self._predict_herg_risk(mol, mw, logp, tpsa, hbd, hba, aromatic_rings)
        
        # Predict Nav liability
        nav_risk = self._predict_nav_risk(mol, mw, logp, tpsa, hbd, hba)
        
        # Predict Ca2+ channel liability
        ca_risk = self._predict_ca_risk(mol, mw, logp, tpsa, hbd, hba, rotatable)
        
        # Calculate overall ion channel risk
        overall_risk = max(herg_risk, nav_risk, ca_risk)
        
        return {
            "herg_risk": herg_risk,
            "nav_risk": nav_risk,
            "ca_risk": ca_risk,
            "overall_ion_channel_risk": overall_risk,
            "risk_level": self._get_risk_level(overall_risk),
            "descriptors": {
                "molecular_weight": mw,
                "logp": logp,
                "tpsa": tpsa,
                "hbd": hbd,
                "hba": hba,
                "aromatic_rings": aromatic_rings,
                "rotatable_bonds": rotatable
            }
        }
    
    def _predict_herg_risk(self, mol: Chem.Mol, mw: float, logp: float, tpsa: float, 
                          hbd: int, hba: int, aromatic_rings: int) -> float:
        """Predict hERG liability risk"""
        # hERG risk factors based on literature
        risk_factors = []
        
        # Lipophilicity (high logP increases hERG risk)
        if logp > 3.0:
            risk_factors.append(min(1.0, (logp - 3.0) / 2.0))
        
        # Molecular weight (larger molecules have higher hERG risk)
        if mw > 400:
            risk_factors.append(min(1.0, (mw - 400) / 200))
        
        # Aromatic rings (more rings increase hERG risk)
        if aromatic_rings > 2:
            risk_factors.append(min(1.0, (aromatic_rings - 2) / 3))
        
        # TPSA (lower TPSA increases hERG risk)
        if tpsa < 90:
            risk_factors.append(min(1.0, (90 - tpsa) / 90))
        
        # Calculate overall hERG risk
        if risk_factors:
            herg_risk = np.mean(risk_factors)
        else:
            herg_risk = 0.1  # Low baseline risk
        
        return min(1.0, herg_risk)
    
    def _predict_nav_risk(self, mol: Chem.Mol, mw: float, logp: float, tpsa: float, 
                         hbd: int, hba: int) -> float:
        """Predict Nav channel liability risk"""
        risk_factors = []
        
        # Similar factors to hERG but with different weights
        if logp > 2.5:
            risk_factors.append(min(1.0, (logp - 2.5) / 2.5))
        
        if mw > 350:
            risk_factors.append(min(1.0, (mw - 350) / 250))
        
        if tpsa < 100:
            risk_factors.append(min(1.0, (100 - tpsa) / 100))
        
        # Nav-specific: high HBD can increase risk
        if hbd > 3:
            risk_factors.append(min(1.0, (hbd - 3) / 5))
        
        if risk_factors:
            nav_risk = np.mean(risk_factors) * 0.8  # Slightly lower than hERG
        else:
            nav_risk = 0.08
        
        return min(1.0, nav_risk)
    
    def _predict_ca_risk(self, mol: Chem.Mol, mw: float, logp: float, tpsa: float, 
                        hbd: int, hba: int, rotatable: int) -> float:
        """Predict Ca2+ channel liability risk"""
        risk_factors = []
        
        # Ca2+ channel risk factors
        if logp > 2.0:
            risk_factors.append(min(1.0, (logp - 2.0) / 3.0))
        
        if mw > 300:
            risk_factors.append(min(1.0, (mw - 300) / 300))
        
        # High rotatable bonds can increase Ca2+ channel risk
        if rotatable > 5:
            risk_factors.append(min(1.0, (rotatable - 5) / 10))
        
        if risk_factors:
            ca_risk = np.mean(risk_factors) * 0.6  # Lower than hERG/Nav
        else:
            ca_risk = 0.05
        
        return min(1.0, ca_risk)
    
    def _get_risk_level(self, risk_score: float) -> str:
        """Get risk level classification"""
        if risk_score >= 0.7:
            return "High"
        elif risk_score >= 0.4:
            return "Medium"
        else:
            return "Low"
    
    def _empty_ion_channel_results(self) -> Dict[str, Any]:
        """Empty results for invalid molecules"""
        return {
            "herg_risk": 0.0,
            "nav_risk": 0.0,
            "ca_risk": 0.0,
            "overall_ion_channel_risk": 0.0,
            "risk_level": "Unknown",
            "descriptors": {}
        }

class SelectivityCalculator:
    """Advanced selectivity calculation and analysis"""
    
    def __init__(self):
        self.logger = structlog.get_logger(__name__)
    
    def calculate_selectivity_score(
        self, 
        on_target_score: float, 
        off_target_scores: Dict[str, float],
        on_target_id: str
    ) -> Dict[str, Any]:
        """Calculate comprehensive selectivity metrics"""
        
        if not off_target_scores:
            return {
                "selectivity_score": 1.0,
                "selectivity_index": float('inf'),
                "off_target_count": 0,
                "risk_distribution": {},
                "selectivity_classification": "Unknown"
            }
        
        # Remove on-target from off-target scores
        off_target_only = {k: v for k, v in off_target_scores.items() if k != on_target_id}
        
        if not off_target_only:
            return {
                "selectivity_score": 1.0,
                "selectivity_index": float('inf'),
                "off_target_count": 0,
                "risk_distribution": {},
                "selectivity_classification": "Unknown"
            }
        
        # Calculate various selectivity metrics
        off_target_values = list(off_target_only.values())
        
        # Selectivity index (on-target / average off-target)
        avg_off_target = np.mean(off_target_values)
        selectivity_index = on_target_score / avg_off_target if avg_off_target > 0 else float('inf')
        
        # Selectivity score (normalized)
        max_off_target = max(off_target_values)
        selectivity_score = max(0.0, min(1.0, (on_target_score - max_off_target) / on_target_score)) if on_target_score > 0 else 0.0
        
        # Risk distribution analysis
        risk_distribution = self._analyze_risk_distribution(off_target_values)
        
        # Selectivity classification
        selectivity_classification = self._classify_selectivity(selectivity_index, selectivity_score)
        
        return {
            "selectivity_score": round(selectivity_score, 4),
            "selectivity_index": round(selectivity_index, 4) if selectivity_index != float('inf') else float('inf'),
            "off_target_count": len(off_target_only),
            "avg_off_target_score": round(avg_off_target, 4),
            "max_off_target_score": round(max_off_target, 4),
            "risk_distribution": risk_distribution,
            "selectivity_classification": selectivity_classification
        }
    
    def _analyze_risk_distribution(self, off_target_scores: List[float]) -> Dict[str, Any]:
        """Analyze distribution of off-target scores"""
        if not off_target_scores:
            return {}
        
        scores = np.array(off_target_scores)
        
        return {
            "mean": round(float(np.mean(scores)), 4),
            "median": round(float(np.median(scores)), 4),
            "std": round(float(np.std(scores)), 4),
            "min": round(float(np.min(scores)), 4),
            "max": round(float(np.max(scores)), 4),
            "quartiles": {
                "q25": round(float(np.percentile(scores, 25)), 4),
                "q50": round(float(np.percentile(scores, 50)), 4),
                "q75": round(float(np.percentile(scores, 75)), 4)
            },
            "high_risk_count": int(np.sum(scores > 0.7)),
            "medium_risk_count": int(np.sum((scores > 0.4) & (scores <= 0.7))),
            "low_risk_count": int(np.sum(scores <= 0.4))
        }
    
    def _classify_selectivity(self, selectivity_index: float, selectivity_score: float) -> str:
        """Classify selectivity based on metrics"""
        if selectivity_index == float('inf'):
            return "Perfect"
        elif selectivity_index > 10:
            return "Excellent"
        elif selectivity_index > 5:
            return "Good"
        elif selectivity_index > 2:
            return "Moderate"
        elif selectivity_index > 1:
            return "Poor"
        else:
            return "Very Poor"

class SafetyProfiler:
    """Comprehensive safety profiling"""
    
    def __init__(self):
        self.logger = structlog.get_logger(__name__)
    
    def profile_safety(self, smiles: str, ion_channel_risks: Dict[str, Any]) -> Dict[str, Any]:
        """Comprehensive safety profiling"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return self._empty_safety_profile()
        
        # Calculate safety-related molecular properties
        safety_properties = self._calculate_safety_properties(mol)
        
        # Assess various safety risks
        safety_risks = {
            "ion_channel": ion_channel_risks,
            "metabolic": self._assess_metabolic_risk(mol),
            "genotoxicity": self._assess_genotoxicity_risk(mol),
            "hepatotoxicity": self._assess_hepatotoxicity_risk(mol),
            "cardiotoxicity": self._assess_cardiotoxicity_risk(mol, ion_channel_risks)
        }
        
        # Calculate overall safety score
        overall_safety = self._calculate_overall_safety(safety_risks)
        
        return {
            "safety_properties": safety_properties,
            "safety_risks": safety_risks,
            "overall_safety_score": overall_safety,
            "safety_classification": self._classify_safety(overall_safety)
        }
    
    def _calculate_safety_properties(self, mol: Chem.Mol) -> Dict[str, Any]:
        """Calculate safety-related molecular properties"""
        return {
            "molecular_weight": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "tpsa": Descriptors.TPSA(mol),
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol),
            "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
            "rings": Descriptors.RingCount(mol),
            "aromatic_rings": Descriptors.NumAromaticRings(mol),
            "sp3_fraction": self._calculate_sp3_fraction(mol),
            "lipinski_violations": self._count_lipinski_violations(mol),
            "veber_violations": self._count_veber_violations(mol)
        }
    
    def _assess_metabolic_risk(self, mol: Chem.Mol) -> Dict[str, Any]:
        """Assess metabolic stability risk"""
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        
        # Metabolic risk factors
        risk_factors = []
        
        if mw > 500:
            risk_factors.append(min(1.0, (mw - 500) / 200))
        
        if logp > 3:
            risk_factors.append(min(1.0, (logp - 3) / 2))
        
        if tpsa < 90:
            risk_factors.append(min(1.0, (90 - tpsa) / 90))
        
        risk_score = np.mean(risk_factors) if risk_factors else 0.1
        
        return {
            "risk_score": round(risk_score, 4),
            "risk_level": "High" if risk_score > 0.6 else "Medium" if risk_score > 0.3 else "Low"
        }
    
    def _assess_genotoxicity_risk(self, mol: Chem.Mol) -> Dict[str, Any]:
        """Assess genotoxicity risk"""
        # Simplified genotoxicity assessment
        aromatic_rings = Descriptors.NumAromaticRings(mol)
        
        risk_score = min(1.0, aromatic_rings * 0.2)
        
        return {
            "risk_score": round(risk_score, 4),
            "risk_level": "High" if risk_score > 0.6 else "Medium" if risk_score > 0.3 else "Low"
        }
    
    def _assess_hepatotoxicity_risk(self, mol: Chem.Mol) -> Dict[str, Any]:
        """Assess hepatotoxicity risk"""
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        
        risk_factors = []
        
        if mw > 400:
            risk_factors.append(min(1.0, (mw - 400) / 200))
        
        if logp > 3.5:
            risk_factors.append(min(1.0, (logp - 3.5) / 1.5))
        
        risk_score = np.mean(risk_factors) if risk_factors else 0.1
        
        return {
            "risk_score": round(risk_score, 4),
            "risk_level": "High" if risk_score > 0.6 else "Medium" if risk_score > 0.3 else "Low"
        }
    
    def _assess_cardiotoxicity_risk(self, mol: Chem.Mol, ion_channel_risks: Dict[str, Any]) -> Dict[str, Any]:
        """Assess cardiotoxicity risk"""
        # Combine ion channel risks with structural factors
        ion_risk = ion_channel_risks.get("overall_ion_channel_risk", 0.0)
        
        # Additional structural factors
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        
        structural_risk = 0.0
        if mw > 450:
            structural_risk += 0.2
        if logp > 4:
            structural_risk += 0.3
        
        # Combined risk
        combined_risk = max(ion_risk, structural_risk)
        
        return {
            "risk_score": round(combined_risk, 4),
            "risk_level": "High" if combined_risk > 0.6 else "Medium" if combined_risk > 0.3 else "Low",
            "ion_channel_contribution": round(ion_risk, 4),
            "structural_contribution": round(structural_risk, 4)
        }
    
    def _calculate_overall_safety(self, safety_risks: Dict[str, Any]) -> float:
        """Calculate overall safety score"""
        risk_scores = []
        
        for risk_type, risk_data in safety_risks.items():
            if isinstance(risk_data, dict) and "risk_score" in risk_data:
                risk_scores.append(risk_data["risk_score"])
        
        if not risk_scores:
            return 0.0
        
        # Use maximum risk as overall safety concern
        max_risk = max(risk_scores)
        safety_score = 1.0 - max_risk  # Convert risk to safety
        
        return round(safety_score, 4)
    
    def _classify_safety(self, safety_score: float) -> str:
        """Classify safety based on score"""
        if safety_score >= 0.8:
            return "Excellent"
        elif safety_score >= 0.6:
            return "Good"
        elif safety_score >= 0.4:
            return "Moderate"
        elif safety_score >= 0.2:
            return "Poor"
        else:
            return "Very Poor"
    
    def _count_lipinski_violations(self, mol: Chem.Mol) -> int:
        """Count Lipinski rule violations"""
        violations = 0
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        
        if mw > 500:
            violations += 1
        if logp > 5:
            violations += 1
        if hbd > 5:
            violations += 1
        if hba > 10:
            violations += 1
        
        return violations
    
    def _count_veber_violations(self, mol: Chem.Mol) -> int:
        """Count Veber rule violations"""
        violations = 0
        rotatable = Descriptors.NumRotatableBonds(mol)
        tpsa = Descriptors.TPSA(mol)
        
        if rotatable > 10:
            violations += 1
        if tpsa > 140:
            violations += 1
        
        return violations
    
    def _calculate_sp3_fraction(self, mol: Chem.Mol) -> float:
        """Calculate sp3 fraction manually since FractionCsp3 may not be available"""
        try:
            # Try to use the built-in function first
            return rdMolDescriptors.FractionCsp3(mol)
        except AttributeError:
            # Fallback calculation
            total_carbons = 0
            sp3_carbons = 0
            
            for atom in mol.GetAtoms():
                if atom.GetSymbol() == 'C':
                    total_carbons += 1
                    # Check if it's sp3 hybridized (4 single bonds)
                    if atom.GetTotalNumHs() + atom.GetDegree() == 4:
                        sp3_carbons += 1
            
            if total_carbons == 0:
                return 0.0
            
            return sp3_carbons / total_carbons
    
    def _empty_safety_profile(self) -> Dict[str, Any]:
        """Empty safety profile for invalid molecules"""
        return {
            "safety_properties": {},
            "safety_risks": {},
            "overall_safety_score": 0.0,
            "safety_classification": "Unknown"
        }

class DecisionEngine:
    """Advanced decision engine with configurable logic"""
    
    def __init__(self):
        self.logger = structlog.get_logger(__name__)
        # Get configuration
        config = config_manager.get_pipeline_config()
        scoring_weights = config.get("scoring_weights", {})
        thresholds = config.get("thresholds", {})
    
    def make_decision(
        self, 
        selectivity_metrics: Dict[str, Any],
        safety_profile: Dict[str, Any],
        ion_channel_risks: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Make synthesis decision based on comprehensive analysis"""
        
        # Extract key metrics
        selectivity_score = selectivity_metrics.get("selectivity_score", 0.0)
        selectivity_index = selectivity_metrics.get("selectivity_index", 0.0)
        safety_score = safety_profile.get("overall_safety_score", 0.0)
        ion_channel_risk = ion_channel_risks.get("overall_ion_channel_risk", 0.0)
        
        # Get decision thresholds
        thresholds = config_manager.get_thresholds()
        
        # Calculate decision score
        decision_score = self._calculate_decision_score(
            selectivity_score, selectivity_index, safety_score, ion_channel_risk
        )
        
        # Make decision
        decision = self._classify_decision(decision_score, thresholds)
        
        # Generate recommendations
        recommendations = self._generate_recommendations(
            decision, selectivity_metrics, safety_profile, ion_channel_risks
        )
        
        return {
            "decision": decision,
            "decision_score": round(decision_score, 4),
            "recommendations": recommendations,
            "key_factors": {
                "selectivity_score": selectivity_score,
                "selectivity_index": selectivity_index,
                "safety_score": safety_score,
                "ion_channel_risk": ion_channel_risk
            }
        }
    
    def _calculate_decision_score(
        self, 
        selectivity_score: float, 
        selectivity_index: float, 
        safety_score: float, 
        ion_channel_risk: float
    ) -> float:
        """Calculate overall decision score"""
        
        # Normalize selectivity index (handle infinity)
        if selectivity_index == float('inf'):
            norm_selectivity_index = 1.0
        else:
            norm_selectivity_index = min(1.0, selectivity_index / 10.0)
        
        # Weighted combination
        weights = {
            "selectivity_score": 0.3,
            "selectivity_index": 0.3,
            "safety_score": 0.25,
            "ion_channel_safety": 0.15
        }
        
        ion_channel_safety = 1.0 - ion_channel_risk
        
        decision_score = (
            weights["selectivity_score"] * selectivity_score +
            weights["selectivity_index"] * norm_selectivity_index +
            weights["safety_score"] * safety_score +
            weights["ion_channel_safety"] * ion_channel_safety
        )
        
        return min(1.0, max(0.0, decision_score))
    
    def _classify_decision(self, decision_score: float, thresholds: Dict[str, float]) -> str:
        """Classify decision based on score and thresholds"""
        synthesize_threshold = thresholds.get("synthesize", 2.0)
        watch_threshold = thresholds.get("watch", 1.0)
        
        if decision_score >= 0.8:
            return "Synthesize"
        elif decision_score >= 0.6:
            return "Watch"
        elif decision_score >= 0.4:
            return "Modify"
        else:
            return "Reject"
    
    def _generate_recommendations(
        self, 
        decision: str, 
        selectivity_metrics: Dict[str, Any],
        safety_profile: Dict[str, Any],
        ion_channel_risks: Dict[str, Any]
    ) -> List[str]:
        """Generate specific recommendations"""
        recommendations = []
        
        if decision == "Synthesize":
            recommendations.append("Proceed with synthesis and testing")
            recommendations.append("Monitor for any unexpected off-target effects")
        
        elif decision == "Watch":
            recommendations.append("Proceed with caution and enhanced monitoring")
            recommendations.append("Consider structural modifications to improve selectivity")
        
        elif decision == "Modify":
            recommendations.append("Structural modifications required before synthesis")
            if selectivity_metrics.get("selectivity_score", 0) < 0.5:
                recommendations.append("Focus on improving target selectivity")
            if safety_profile.get("overall_safety_score", 0) < 0.6:
                recommendations.append("Address safety concerns through structural optimization")
        
        else:  # Reject
            recommendations.append("Do not proceed with current compound")
            recommendations.append("Consider alternative chemical scaffolds")
        
        # Add specific recommendations based on analysis
        if ion_channel_risks.get("herg_risk", 0) > 0.5:
            recommendations.append("High hERG risk - consider structural modifications to reduce ion channel binding")
        
        if selectivity_metrics.get("off_target_count", 0) > 10:
            recommendations.append("High number of off-targets - consider more selective analogs")
        
        return recommendations

async def estimate_impact_risk(
    empirical_scores: Dict[str, float],
    structural_scores: Dict[str, float],
    on_target_id: str,
    smiles: str,
    output_path: str = "impact_risk/impact_summary.json",
    detailed_path: str = "impact_risk/detailed_analysis.json"
) -> Dict[str, Any]:
    """
    Enhanced impact risk assessment with comprehensive selectivity and safety analysis.
    
    Args:
        empirical_scores: Empirical binding scores
        structural_scores: Structural binding scores
        on_target_id: Primary target UniProt ID
        smiles: Compound SMILES
        output_path: Path for summary output
        detailed_path: Path for detailed analysis

    Returns:
        Comprehensive impact assessment results
    """
    
    logger = structlog.get_logger(__name__)
    logger.info(f"Starting enhanced IMPACT risk assessment for {on_target_id}")
    
    # Initialize analyzers
    ion_channel_predictor = IonChannelPredictor()
    selectivity_calculator = SelectivityCalculator()
    safety_profiler = SafetyProfiler()
    decision_engine = DecisionEngine()
    
    # Step 1: Ion channel liability assessment
    logger.info("Assessing ion channel liability...")
    ion_channel_risks = ion_channel_predictor.predict_ion_channel_risk(smiles)
    
    # Step 2: Combine empirical and structural scores
    combined_scores = {}
    all_targets = set(empirical_scores.keys()) | set(structural_scores.keys())
    
    for target in all_targets:
        emp_score = empirical_scores.get(target, 0.0)
        struct_score = structural_scores.get(target, 0.0)
        
        # Weighted combination (configurable)
        combined_score = 0.6 * emp_score + 0.4 * struct_score
        combined_scores[target] = round(combined_score, 4)
    
    # Step 3: Calculate selectivity metrics
    logger.info("Calculating selectivity metrics...")
    on_target_score = combined_scores.get(on_target_id, 0.0)
    off_target_scores = {k: v for k, v in combined_scores.items() if k != on_target_id}
    
    selectivity_metrics = selectivity_calculator.calculate_selectivity_score(
        on_target_score, off_target_scores, on_target_id
    )
    
    # Step 4: Safety profiling
    logger.info("Profiling safety characteristics...")
    safety_profile = safety_profiler.profile_safety(smiles, ion_channel_risks)
    
    # Step 5: Decision making
    logger.info("Making synthesis decision...")
    decision_result = decision_engine.make_decision(
        selectivity_metrics, safety_profile, ion_channel_risks
    )
    
    # Compile comprehensive results
    comprehensive_results = {
        "compound_info": {
            "smiles": smiles,
            "on_target_id": on_target_id,
            "on_target_score": on_target_score
        },
        "selectivity_analysis": {
            "metrics": selectivity_metrics,
            "off_target_scores": off_target_scores,
            "combined_scores": combined_scores
        },
        "safety_analysis": {
            "ion_channel_risks": ion_channel_risks,
            "safety_profile": safety_profile
        },
        "decision": decision_result,
        "metadata": {
            "analysis_timestamp": pd.Timestamp.now().isoformat(),
            "analysis_method": "enhanced_impact"
        }
    }
    
    # Save detailed analysis
    os.makedirs(os.path.dirname(detailed_path), exist_ok=True)
    with open(detailed_path, 'w') as f:
        json.dump(comprehensive_results, f, indent=2, default=str)
    
    # Create summary for legacy compatibility
    summary = {
        "on_target_id": on_target_id,
        "on_target_score": on_target_score,
        "avg_off_target_score": selectivity_metrics.get("avg_off_target_score", 0.0),
        "selectivity_score": selectivity_metrics.get("selectivity_score", 0.0),
        "selectivity_index": selectivity_metrics.get("selectivity_index", 0.0),
        "decision_flag": decision_result["decision"],
        "risky_offtargets": [
            {"uniprot_id": uid, "combined_score": score}
            for uid, score in sorted(off_target_scores.items(), key=lambda x: x[1], reverse=True)[:10]
        ],
        "ion_channel_risks": ion_channel_risks,
        "safety_score": safety_profile.get("overall_safety_score", 0.0),
        "recommendations": decision_result.get("recommendations", [])
    }
    
    # Save summary
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(summary, f, indent=2)
    
    logger.info(f"✅ Enhanced IMPACT risk assessment complete")
    logger.info(f"📁 Summary: {output_path}")
    logger.info(f"📄 Detailed analysis: {detailed_path}")
    logger.info(f"🎯 Decision: {decision_result['decision']}")
    
    return summary
