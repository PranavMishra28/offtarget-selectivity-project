"""
Enhanced SPARROW - Advanced Synthesis Feasibility Assessment
Implements route planning, feasibility scoring, and synthetic accessibility using multiple APIs.
"""

import os
import json
import asyncio
from typing import List, Dict, Any, Optional, Tuple
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import AllChem
import structlog
from tqdm import tqdm

# Import configuration and API clients
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.config_manager import config_manager
from utils.api_client import api_manager

class SynthesisRoute:
    """Represents a synthesis route with scoring"""
    
    def __init__(self, route_data: Dict[str, Any]):
        self.steps = route_data.get("steps", [])
        self.total_cost = route_data.get("total_cost", 0.0)
        self.num_steps = len(self.steps)
        self.success_probability = route_data.get("success_probability", 0.5)
        self.reagent_availability = route_data.get("reagent_availability", 0.5)
        self.complexity_score = self._calculate_complexity()
        
    def _calculate_complexity(self) -> float:
        """Calculate route complexity score"""
        if not self.steps:
            return 1.0
        
        # Factors: number of steps, step complexity, yield expectations
        step_complexities = []
        for step in self.steps:
            complexity = 0.5  # Base complexity
            
            # Add complexity for different reaction types
            reaction_type = step.get("reaction_type", "").lower()
            if "c-c coupling" in reaction_type:
                complexity += 0.3
            elif "reduction" in reaction_type:
                complexity += 0.1
            elif "oxidation" in reaction_type:
                complexity += 0.2
            elif "protection" in reaction_type:
                complexity += 0.2
            elif "deprotection" in reaction_type:
                complexity += 0.1
            
            step_complexities.append(complexity)
        
        return np.mean(step_complexities) if step_complexities else 1.0
    
    def get_feasibility_score(self) -> float:
        """Calculate overall feasibility score"""
        # Weighted combination of factors
        weights = {
            "steps": 0.3,
            "cost": 0.2,
            "success_prob": 0.25,
            "availability": 0.15,
            "complexity": 0.1
        }
        
        # Normalize factors (lower is better for most)
        step_score = max(0, 1 - (self.num_steps - 1) / 10)  # Penalize more steps
        cost_score = max(0, 1 - self.total_cost / 1000)  # Penalize high cost
        complexity_score = max(0, 1 - self.complexity_score)
        
        feasibility = (
            weights["steps"] * step_score +
            weights["cost"] * cost_score +
            weights["success_prob"] * self.success_probability +
            weights["availability"] * self.reagent_availability +
            weights["complexity"] * complexity_score
        )
        
        return min(1.0, max(0.0, feasibility))

class SynthesisAnalyzer:
    """Analyzes synthesis feasibility using multiple approaches"""
    
    def __init__(self):
        self.logger = structlog.get_logger(__name__)
        self.askcos_client = api_manager.get_client("askcos")
        self.ibm_rxn_client = api_manager.get_client("ibm_rxn")
        
    async def analyze_synthesis(self, smiles: str) -> Dict[str, Any]:
        """Analyze synthesis feasibility using multiple APIs"""
        results = {
            "smiles": smiles,
            "routes": [],
            "feasibility_score": 0.0,
            "recommendation": "Unknown",
            "metadata": {}
        }
        
        # Try ASKCOS first
        if self.askcos_client:
            try:
                askcos_result = await self.askcos_client.get_retrosynthesis_routes(smiles)
                if askcos_result.success:
                    results["routes"].extend(self._parse_askcos_routes(askcos_result.data))
                    results["metadata"]["askcos_available"] = True
                else:
                    results["metadata"]["askcos_error"] = askcos_result.error
            except Exception as e:
                self.logger.warning(f"ASKCOS analysis failed: {e}")
                results["metadata"]["askcos_error"] = str(e)
        
        # Try IBM RXN as backup
        if self.ibm_rxn_client and not results["routes"]:
            try:
                ibm_result = await self.ibm_rxn_client.get_retrosynthesis_routes(smiles)
                if ibm_result.success:
                    results["routes"].extend(self._parse_ibm_routes(ibm_result.data))
                    results["metadata"]["ibm_rxn_available"] = True
                else:
                    results["metadata"]["ibm_rxn_error"] = ibm_result.error
            except Exception as e:
                self.logger.warning(f"IBM RXN analysis failed: {e}")
                results["metadata"]["ibm_rxn_error"] = str(e)
        
        # Fallback to rule-based scoring
        if not results["routes"]:
            results["routes"] = self._rule_based_routes(smiles)
            results["metadata"]["rule_based"] = True
        
        # Calculate overall feasibility
        if results["routes"]:
            route_scores = [route.get_feasibility_score() for route in results["routes"]]
            results["feasibility_score"] = max(route_scores)  # Best route score
            results["recommendation"] = self._get_recommendation(results["feasibility_score"])
        
        return results
    
    def _parse_askcos_routes(self, data: Dict[str, Any]) -> List[SynthesisRoute]:
        """Parse ASKCOS route data"""
        routes = []
        try:
            trees = data.get("trees", [])
            for tree in trees:
                route_data = {
                    "steps": [],
                    "total_cost": tree.get("ppg", 0.0),
                    "success_probability": 0.7,  # Default
                    "reagent_availability": 0.8   # Default
                }
                
                # Parse tree structure to extract steps
                if "children" in tree:
                    for child in tree["children"]:
                        step = {
                            "reaction_type": child.get("reaction_name", "Unknown"),
                            "reagents": child.get("reagents", []),
                            "yield": child.get("yield", 0.8)
                        }
                        route_data["steps"].append(step)
                
                routes.append(SynthesisRoute(route_data))
        except Exception as e:
            self.logger.error(f"Failed to parse ASKCOS routes: {e}")
        
        return routes
    
    def _parse_ibm_routes(self, data: Dict[str, Any]) -> List[SynthesisRoute]:
        """Parse IBM RXN route data"""
        routes = []
        try:
            # IBM RXN format parsing
            route_data = {
                "steps": [],
                "total_cost": data.get("cost", 0.0),
                "success_probability": 0.6,
                "reagent_availability": 0.7
            }
            
            # Parse steps if available
            if "steps" in data:
                for step in data["steps"]:
                    route_data["steps"].append({
                        "reaction_type": step.get("type", "Unknown"),
                        "reagents": step.get("reagents", []),
                        "yield": step.get("yield", 0.8)
                    })
            
            routes.append(SynthesisRoute(route_data))
        except Exception as e:
            self.logger.error(f"Failed to parse IBM RXN routes: {e}")
        
        return routes
    
    def _rule_based_routes(self, smiles: str) -> List[SynthesisRoute]:
        """Generate rule-based synthesis routes when APIs are unavailable"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return []
        
        # Calculate molecular complexity
        complexity = self._calculate_molecular_complexity(mol)
        
        # Estimate synthesis difficulty based on molecular properties
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        rings = Descriptors.RingCount(mol)
        rotatable = Descriptors.NumRotatableBonds(mol)
        
        # Estimate number of steps based on complexity
        estimated_steps = max(1, int(complexity * 3))
        
        # Estimate cost based on molecular weight and complexity
        estimated_cost = mw * complexity * 10
        
        # Calculate success probability based on properties
        success_prob = max(0.1, min(0.9, 1 - complexity * 0.5))
        
        route_data = {
            "steps": [{"reaction_type": "Estimated", "reagents": [], "yield": 0.8}] * estimated_steps,
            "total_cost": estimated_cost,
            "success_probability": success_prob,
            "reagent_availability": 0.8
        }
        
        return [SynthesisRoute(route_data)]
    
    def _calculate_molecular_complexity(self, mol: Chem.Mol) -> float:
        """Calculate molecular complexity score"""
        complexity = 0.0
        
        # Ring complexity
        rings = Descriptors.RingCount(mol)
        aromatic_rings = Descriptors.NumAromaticRings(mol)
        complexity += rings * 0.1 + aromatic_rings * 0.2
        
        # Functional group complexity
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        complexity += (hbd + hba) * 0.05
        
        # Stereochemistry
        chiral_centers = len(Chem.FindMolChiralCenters(mol))
        complexity += chiral_centers * 0.3
        
        # Size complexity
        mw = Descriptors.MolWt(mol)
        complexity += mw / 1000
        
        return min(1.0, complexity)
    
    def _get_recommendation(self, feasibility_score: float) -> str:
        """Get synthesis recommendation based on feasibility score"""
        if feasibility_score >= 0.8:
            return "Highly Feasible"
        elif feasibility_score >= 0.6:
            return "Moderately Feasible"
        elif feasibility_score >= 0.4:
            return "Challenging"
        else:
            return "Very Difficult"

class CheminformaticsScorer:
    """Advanced cheminformatics-based synthesis scoring"""
    
    def __init__(self):
        self.logger = structlog.get_logger(__name__)
    
    def calculate_synthesis_score(self, smiles: str) -> Dict[str, float]:
        """Calculate comprehensive synthesis accessibility score"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"synthesis_score": 0.0, "complexity": 1.0}
        
        # Calculate various molecular properties
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rotatable = Descriptors.NumRotatableBonds(mol)
        rings = Descriptors.RingCount(mol)
        aromatic_rings = Descriptors.NumAromaticRings(mol)
        tpsa = Descriptors.TPSA(mol)
        
        # Calculate complexity factors
        complexity_factors = {
            "size": min(1.0, mw / 500),  # Normalize by typical drug MW
            "lipinski_violations": self._count_lipinski_violations(mol),
            "ring_complexity": rings * 0.1 + aromatic_rings * 0.2,
            "functional_groups": (hbd + hba) * 0.05,
            "flexibility": rotatable * 0.02,
            "polarity": tpsa / 200  # Normalize by typical TPSA
        }
        
        # Calculate overall complexity
        complexity = sum(complexity_factors.values()) / len(complexity_factors)
        
        # Calculate synthesis score (inverse of complexity)
        synthesis_score = max(0.0, 1.0 - complexity)
        
        return {
            "synthesis_score": synthesis_score,
            "complexity": complexity,
            "factors": complexity_factors
        }
    
    def _count_lipinski_violations(self, mol: Chem.Mol) -> float:
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
        
        return violations / 4.0  # Normalize to 0-1

async def compute_prioritization(
    smiles_list: List[str], 
    output_path: str = "sparrow/ranked_candidates.csv",
    metadata_path: str = "sparrow/synthesis_analysis.json"
) -> Dict[str, Any]:
    """
    Enhanced synthesis prioritization with multiple analysis methods.
    
    Args:
        smiles_list: List of SMILES strings to analyze
        output_path: Path for CSV output
        metadata_path: Path for detailed analysis metadata
    """
    
    logger = structlog.get_logger(__name__)
    logger.info(f"Starting enhanced SPARROW analysis for {len(smiles_list)} molecules")
    
    # Initialize analyzers
    synthesis_analyzer = SynthesisAnalyzer()
    cheminformatics_scorer = CheminformaticsScorer()
    
    # Analyze each molecule
    results = []
    detailed_analysis = {}
    
    for i, smiles in enumerate(tqdm(smiles_list, desc="Analyzing synthesis feasibility")):
        try:
            # Get cheminformatics score
            chem_score = cheminformatics_scorer.calculate_synthesis_score(smiles)
            
            # Get synthesis analysis (async)
            synthesis_analysis = await synthesis_analyzer.analyze_synthesis(smiles)
            
            # Combine scores
            combined_score = (
                0.4 * synthesis_analysis["feasibility_score"] +
                0.6 * chem_score["synthesis_score"]
            )
            
            # Create result entry
            result = {
                "smiles": smiles,
                "synthesis_score": round(combined_score, 4),
                "feasibility_score": round(synthesis_analysis["feasibility_score"], 4),
                "cheminformatics_score": round(chem_score["synthesis_score"], 4),
                "complexity": round(chem_score["complexity"], 4),
                "recommendation": synthesis_analysis["recommendation"],
                "num_routes": len(synthesis_analysis["routes"]),
                "molecular_weight": round(Descriptors.MolWt(Chem.MolFromSmiles(smiles)), 2),
                "logp": round(Descriptors.MolLogP(Chem.MolFromSmiles(smiles)), 2)
            }
            
            results.append(result)
            # Convert SynthesisRoute objects to dictionaries for JSON serialization
            serializable_synthesis = {
                "smiles": synthesis_analysis["smiles"],
                "feasibility_score": synthesis_analysis["feasibility_score"],
                "recommendation": synthesis_analysis["recommendation"],
                "metadata": synthesis_analysis["metadata"],
                "routes": [
                    {
                        "steps": route.steps,
                        "total_cost": route.total_cost,
                        "num_steps": route.num_steps,
                        "success_probability": route.success_probability,
                        "reagent_availability": route.reagent_availability,
                        "complexity_score": route.complexity_score,
                        "feasibility_score": route.get_feasibility_score()
                    }
                    for route in synthesis_analysis["routes"]
                ]
            }
            
            detailed_analysis[smiles] = {
                "synthesis_analysis": serializable_synthesis,
                "cheminformatics_analysis": chem_score
            }
            
        except Exception as e:
            logger.error(f"Failed to analyze {smiles}: {e}")
            # Add fallback result
            results.append({
                "smiles": smiles,
                "synthesis_score": 0.0,
                "feasibility_score": 0.0,
                "cheminformatics_score": 0.0,
                "complexity": 1.0,
                "recommendation": "Analysis Failed",
                "num_routes": 0,
                "molecular_weight": 0.0,
                "logp": 0.0
            })
    
    # Sort by synthesis score (higher is better)
    results.sort(key=lambda x: x["synthesis_score"], reverse=True)
    
    # Save results
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df = pd.DataFrame(results)
    df.to_csv(output_path, index=False)
    
    # Save detailed analysis
    with open(metadata_path, 'w') as f:
        json.dump(detailed_analysis, f, indent=2)
    
    logger.info(f"✅ Enhanced SPARROW analysis complete")
    logger.info(f"📁 Results: {output_path}")
    logger.info(f"📄 Detailed analysis: {metadata_path}")
    
    return {
        "csv_path": output_path,
        "metadata_path": metadata_path,
        "num_analyzed": len(results),
        "top_recommendations": results[:5]
    }
