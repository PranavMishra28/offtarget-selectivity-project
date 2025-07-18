"""
Toxicophore Detection Module
Identifies known toxic structural patterns and alerts for potential safety concerns.
"""

import os
import json
from typing import Dict, Any, List, Optional, Tuple
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np
import structlog
from dataclasses import dataclass

@dataclass
class ToxicophoreAlert:
    """Represents a toxicophore alert"""
    name: str
    description: str
    severity: str  # "High", "Medium", "Low"
    smarts_pattern: str
    matched_atoms: List[int]
    risk_score: float
    recommendations: List[str]

class ToxicophoreDetector:
    """Advanced toxicophore detection system"""
    
    def __init__(self):
        self.logger = structlog.get_logger(__name__)
        self.toxicophores = self._load_toxicophore_patterns()
        
    def _load_toxicophore_patterns(self) -> Dict[str, Dict[str, Any]]:
        """Load toxicophore patterns from database"""
        return {
            "hERG_blockers": {
                "patterns": [
                    "c1ccc(cc1)-c2ccccc2",  # Biphenyl
                    "c1ccc(cc1)-c2cccnc2",   # Phenylpyridine
                    "c1ccc(cc1)-c2ccncc2",   # Phenylpyrimidine
                    "c1ccc(cc1)-c2cccnc2",   # Phenylquinoline
                ],
                "description": "hERG potassium channel blockers",
                "severity": "High",
                "recommendations": [
                    "Consider structural modifications to reduce lipophilicity",
                    "Add polar groups to reduce membrane permeability",
                    "Monitor for QT prolongation in preclinical studies"
                ]
            },
            "genotoxic_alert": {
                "patterns": [
                    "c1ccc(cc1)N(=O)=O",     # Nitroaromatic
                    "c1ccc(cc1)N=Nc2ccccc2",  # Azo compounds
                    "c1ccc(cc1)C(=O)N=O",     # Nitroso compounds
                    "c1ccc(cc1)C(=O)Cl",      # Acyl halides
                ],
                "description": "Potential genotoxic compounds",
                "severity": "High",
                "recommendations": [
                    "Perform Ames test for mutagenicity",
                    "Consider alternative scaffolds",
                    "Add detoxifying substituents"
                ]
            },
            "reactive_electrophiles": {
                "patterns": [
                    "C(=O)Cl",               # Acid chlorides
                    "C(=O)Br",               # Acid bromides
                    "C(=O)I",                # Acid iodides
                    "C(=O)OC(=O)",           # Anhydrides
                    "C(=O)SC(=O)",           # Thioanhydrides
                ],
                "description": "Reactive electrophilic compounds",
                "severity": "High",
                "recommendations": [
                    "Avoid reactive electrophiles in drug candidates",
                    "Consider prodrug approaches",
                    "Use protective groups in synthesis"
                ]
            },
            "metabolic_alert": {
                "patterns": [
                    "c1ccc(cc1)OC",          # Aryl ethers
                    "c1ccc(cc1)SC",          # Aryl thioethers
                    "c1ccc(cc1)NC",          # Aryl amines
                    "c1ccc(cc1)CC",          # Aryl alkyl
                ],
                "description": "Metabolic liability alerts",
                "severity": "Medium",
                "recommendations": [
                    "Monitor metabolic stability",
                    "Consider blocking metabolic sites",
                    "Perform microsomal stability studies"
                ]
            },
            "phospholipidosis": {
                "patterns": [
                    "c1ccc(cc1)CCCCCC",      # Long alkyl chains
                    "c1ccc(cc1)CCCCCCCC",    # Very long alkyl chains
                    "c1ccc(cc1)CCCCCCCCCC",  # Extremely long alkyl chains
                ],
                "description": "Phospholipidosis inducers",
                "severity": "Medium",
                "recommendations": [
                    "Reduce lipophilicity",
                    "Add polar groups",
                    "Monitor for phospholipidosis in preclinical studies"
                ]
            },
            "hepatotoxicity": {
                "patterns": [
                    "c1ccc(cc1)C(=O)c2ccccc2",  # Benzophenone
                    "c1ccc(cc1)C(=O)c2cccnc2",  # Aryl ketones
                    "c1ccc(cc1)C(=O)OC",        # Aryl esters
                    "c1ccc(cc1)C(=O)SC",        # Aryl thioesters
                ],
                "description": "Hepatotoxicity alerts",
                "severity": "Medium",
                "recommendations": [
                    "Monitor liver function in preclinical studies",
                    "Consider alternative scaffolds",
                    "Add hepatoprotective groups"
                ]
            }
        }
    
    def detect_toxicophores(self, smiles: str) -> List[ToxicophoreAlert]:
        """Detect toxicophores in a compound"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return []
        
        alerts = []
        
        for category, data in self.toxicophores.items():
            for pattern in data["patterns"]:
                matches = self._find_pattern_matches(mol, pattern)
                if matches:
                    alert = ToxicophoreAlert(
                        name=category,
                        description=data["description"],
                        severity=data["severity"],
                        smarts_pattern=pattern,
                        matched_atoms=matches,
                        risk_score=self._calculate_risk_score(mol, pattern, data["severity"]),
                        recommendations=data["recommendations"]
                    )
                    alerts.append(alert)
        
        return alerts
    
    def _find_pattern_matches(self, mol: Chem.Mol, pattern: str) -> List[int]:
        """Find atoms matching a SMARTS pattern"""
        try:
            pattern_mol = Chem.MolFromSmarts(pattern)
            if pattern_mol is None:
                return []
            
            matches = mol.GetSubstructMatches(pattern_mol)
            if matches:
                # Return atom indices from first match
                return list(matches[0])
            return []
        except Exception as e:
            self.logger.warning(f"Pattern matching failed for {pattern}: {e}")
            return []
    
    def _calculate_risk_score(self, mol: Chem.Mol, pattern: str, severity: str) -> float:
        """Calculate risk score for a toxicophore"""
        base_scores = {"High": 0.8, "Medium": 0.5, "Low": 0.2}
        base_score = base_scores.get(severity, 0.5)
        
        # Adjust based on molecular properties
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        
        # Higher risk for larger, more lipophilic compounds
        mw_factor = min(1.0, mw / 500)
        logp_factor = min(1.0, max(0, logp) / 5)
        
        risk_score = base_score * (0.6 + 0.2 * mw_factor + 0.2 * logp_factor)
        return min(1.0, risk_score)
    
    def generate_toxicophore_report(self, smiles: str) -> Dict[str, Any]:
        """Generate comprehensive toxicophore report"""
        alerts = self.detect_toxicophores(smiles)
        
        # Calculate overall risk
        if alerts:
            max_risk = max(alert.risk_score for alert in alerts)
            avg_risk = np.mean([alert.risk_score for alert in alerts])
            high_severity_count = len([a for a in alerts if a.severity == "High"])
        else:
            max_risk = 0.0
            avg_risk = 0.0
            high_severity_count = 0
        
        # Generate recommendations
        all_recommendations = []
        for alert in alerts:
            all_recommendations.extend(alert.recommendations)
        unique_recommendations = list(set(all_recommendations))
        
        report = {
            "compound_smiles": smiles,
            "total_alerts": len(alerts),
            "high_severity_alerts": high_severity_count,
            "max_risk_score": max_risk,
            "average_risk_score": avg_risk,
            "overall_assessment": self._get_overall_assessment(max_risk, high_severity_count),
            "alerts": [
                {
                    "name": alert.name,
                    "description": alert.description,
                    "severity": alert.severity,
                    "risk_score": alert.risk_score,
                    "recommendations": alert.recommendations
                }
                for alert in alerts
            ],
            "recommendations": unique_recommendations,
            "metadata": {
                "analysis_timestamp": pd.Timestamp.now().isoformat(),
                "version": "1.0"
            }
        }
        
        return report
    
    def _get_overall_assessment(self, max_risk: float, high_severity_count: int) -> str:
        """Get overall assessment based on alerts"""
        if max_risk > 0.8 or high_severity_count > 2:
            return "High Risk - Significant toxicophore alerts detected"
        elif max_risk > 0.5 or high_severity_count > 0:
            return "Medium Risk - Some toxicophore alerts detected"
        else:
            return "Low Risk - No significant toxicophore alerts"
    
    def suggest_modifications(self, smiles: str) -> List[str]:
        """Suggest structural modifications to reduce toxicity"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ["Invalid SMILES structure"]
        
        suggestions = []
        
        # Check for specific toxicophores and suggest modifications
        alerts = self.detect_toxicophores(smiles)
        
        for alert in alerts:
            if "hERG_blockers" in alert.name:
                suggestions.append("Add polar substituents to reduce lipophilicity")
                suggestions.append("Introduce hydrogen bond acceptors")
                suggestions.append("Reduce aromatic ring count")
            
            elif "genotoxic_alert" in alert.name:
                suggestions.append("Replace nitro groups with alternative electron-withdrawing groups")
                suggestions.append("Add detoxifying substituents")
                suggestions.append("Consider alternative scaffolds")
            
            elif "reactive_electrophiles" in alert.name:
                suggestions.append("Replace reactive groups with stable alternatives")
                suggestions.append("Use prodrug approaches")
                suggestions.append("Add protective groups")
        
        # General suggestions based on molecular properties
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        
        if mw > 500:
            suggestions.append("Reduce molecular weight by removing non-essential groups")
        
        if logp > 3:
            suggestions.append("Add polar groups to reduce lipophilicity")
        
        if hbd > 5:
            suggestions.append("Reduce hydrogen bond donors")
        
        if hba > 10:
            suggestions.append("Reduce hydrogen bond acceptors")
        
        return list(set(suggestions))  # Remove duplicates

async def analyze_toxicophores(
    smiles: str,
    output_path: str = "toxicity/toxicophore_analysis.json"
) -> Dict[str, Any]:
    """
    Analyze compound for toxicophores and generate report.
    
    Args:
        smiles: Compound SMILES string
        output_path: Path for output file
    
    Returns:
        Comprehensive toxicophore analysis report
    """
    
    logger = structlog.get_logger(__name__)
    logger.info(f"Starting toxicophore analysis for compound")
    
    detector = ToxicophoreDetector()
    report = detector.generate_toxicophore_report(smiles)
    
    # Save report
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(report, f, indent=2)
    
    logger.info(f"✅ Toxicophore analysis complete")
    logger.info(f"📁 Report: {output_path}")
    logger.info(f"🚨 Alerts: {report['total_alerts']} detected")
    
    return report 