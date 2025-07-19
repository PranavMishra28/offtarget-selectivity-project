"""
Enhanced Conflict Resolution - Model Agreement Analysis
Implements comprehensive conflict resolution between different model predictions and provides consensus recommendations.
"""

import os
import json
import asyncio
from typing import List, Dict, Any, Optional, Tuple
import pandas as pd
import numpy as np
from rdkit import Chem
import structlog
from tqdm import tqdm

# Import configuration
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.config_manager import config_manager

class ModelPrediction:
    """Represents a prediction from a specific model"""
    
    def __init__(self, model_name: str, prediction_data: Dict[str, Any]):
        self.model_name = model_name
        self.prediction = prediction_data.get("prediction", "")
        self.confidence = prediction_data.get("confidence", 0.0)
        self.score = prediction_data.get("score", 0.0)
        self.metadata = prediction_data.get("metadata", {})
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            "model_name": self.model_name,
            "prediction": self.prediction,
            "confidence": self.confidence,
            "score": self.score,
            "metadata": self.metadata
        }

class ConflictAnalyzer:
    """Analyzes conflicts between model predictions"""
    
    def __init__(self):
        self.logger = structlog.get_logger(__name__)
    
    def analyze_conflicts(self, predictions: List[ModelPrediction]) -> Dict[str, Any]:
        """Analyze conflicts between model predictions"""
        try:
            if not predictions:
                return self._empty_conflict_analysis()
            
            # Group predictions by type
            prediction_groups = self._group_predictions(predictions)
            
            # Calculate agreement metrics
            agreement_metrics = self._calculate_agreement_metrics(prediction_groups)
            
            # Identify conflicts
            conflicts = self._identify_conflicts(prediction_groups)
            
            # Generate consensus
            consensus = self._generate_consensus(predictions, agreement_metrics)
            
            return {
                "predictions": [pred.to_dict() for pred in predictions],
                "prediction_groups": prediction_groups,
                "agreement_metrics": agreement_metrics,
                "conflicts": conflicts,
                "consensus": consensus,
                "analysis_timestamp": pd.Timestamp.now().isoformat()
            }
            
        except Exception as e:
            self.logger.error(f"Conflict analysis failed: {e}")
            return self._empty_conflict_analysis()
    
    def _group_predictions(self, predictions: List[ModelPrediction]) -> Dict[str, List[ModelPrediction]]:
        """Group predictions by type"""
        groups = {}
        
        for pred in predictions:
            pred_type = self._categorize_prediction(pred.prediction)
            if pred_type not in groups:
                groups[pred_type] = []
            groups[pred_type].append(pred)
        
        return groups
    
    def _categorize_prediction(self, prediction: str) -> str:
        """Categorize prediction type"""
        prediction_lower = prediction.lower()
        
        if "synthesize" in prediction_lower or "proceed" in prediction_lower:
            return "synthesize"
        elif "reject" in prediction_lower or "stop" in prediction_lower:
            return "reject"
        elif "modify" in prediction_lower or "optimize" in prediction_lower:
            return "modify"
        elif "caution" in prediction_lower or "watch" in prediction_lower:
            return "caution"
        else:
            return "unknown"
    
    def _calculate_agreement_metrics(self, prediction_groups: Dict[str, List[ModelPrediction]]) -> Dict[str, Any]:
        """Calculate agreement metrics between models"""
        metrics = {
            "total_predictions": sum(len(group) for group in prediction_groups.values()),
            "group_distribution": {k: len(v) for k, v in prediction_groups.items()},
            "agreement_score": 0.0,
            "consensus_type": "none",
            "confidence_range": {"min": 0.0, "max": 0.0, "mean": 0.0}
        }
        
        if not prediction_groups:
            return metrics
        
        # Calculate agreement score
        total_preds = metrics["total_predictions"]
        if total_preds > 0:
            # Find the most common prediction type
            max_group = max(prediction_groups.items(), key=lambda x: len(x[1]))
            consensus_count = len(max_group[1])
            metrics["agreement_score"] = consensus_count / total_preds
            metrics["consensus_type"] = max_group[0]
        
        # Calculate confidence statistics
        all_confidences = []
        for group in prediction_groups.values():
            all_confidences.extend([pred.confidence for pred in group])
        
        if all_confidences:
            metrics["confidence_range"] = {
                "min": min(all_confidences),
                "max": max(all_confidences),
                "mean": np.mean(all_confidences)
            }
        
        return metrics
    
    def _identify_conflicts(self, prediction_groups: Dict[str, List[ModelPrediction]]) -> List[Dict[str, Any]]:
        """Identify specific conflicts between models"""
        conflicts = []
        
        if len(prediction_groups) <= 1:
            return conflicts
        
        # Find conflicting groups
        group_names = list(prediction_groups.keys())
        for i in range(len(group_names)):
            for j in range(i + 1, len(group_names)):
                group1_name = group_names[i]
                group2_name = group_names[j]
                
                if self._are_conflicting(group1_name, group2_name):
                    conflict = {
                        "type": f"{group1_name}_vs_{group2_name}",
                        "group1": {
                            "type": group1_name,
                            "models": [pred.model_name for pred in prediction_groups[group1_name]],
                            "count": len(prediction_groups[group1_name])
                        },
                        "group2": {
                            "type": group2_name,
                            "models": [pred.model_name for pred in prediction_groups[group2_name]],
                            "count": len(prediction_groups[group2_name])
                        },
                        "severity": self._calculate_conflict_severity(
                            prediction_groups[group1_name],
                            prediction_groups[group2_name]
                        )
                    }
                    conflicts.append(conflict)
        
        return conflicts
    
    def _are_conflicting(self, type1: str, type2: str) -> bool:
        """Check if two prediction types are conflicting"""
        conflicting_pairs = [
            ("synthesize", "reject"),
            ("synthesize", "modify"),
            ("reject", "modify")
        ]
        
        return (type1, type2) in conflicting_pairs or (type2, type1) in conflicting_pairs
    
    def _calculate_conflict_severity(self, group1: List[ModelPrediction], group2: List[ModelPrediction]) -> str:
        """Calculate conflict severity"""
        total_models = len(group1) + len(group2)
        max_group_size = max(len(group1), len(group2))
        
        ratio = max_group_size / total_models
        
        if ratio >= 0.8:
            return "low"
        elif ratio >= 0.6:
            return "medium"
        else:
            return "high"
    
    def _generate_consensus(self, predictions: List[ModelPrediction], agreement_metrics: Dict[str, Any]) -> Dict[str, Any]:
        """Generate consensus recommendation"""
        consensus_type = agreement_metrics["consensus_type"]
        agreement_score = agreement_metrics["agreement_score"]
        
        # Determine consensus recommendation
        if agreement_score >= 0.8:
            recommendation = consensus_type
            confidence = "high"
        elif agreement_score >= 0.6:
            recommendation = consensus_type
            confidence = "medium"
        else:
            # No clear consensus, use weighted approach
            recommendation = self._weighted_consensus(predictions)
            confidence = "low"
        
        return {
            "recommendation": recommendation,
            "confidence": confidence,
            "agreement_score": agreement_score,
            "reasoning": self._generate_reasoning(predictions, consensus_type, agreement_score)
        }
    
    def _weighted_consensus(self, predictions: List[ModelPrediction]) -> str:
        """Generate weighted consensus when no clear majority"""
        # Weight by confidence
        weighted_scores = {}
        
        for pred in predictions:
            pred_type = self._categorize_prediction(pred.prediction)
            if pred_type not in weighted_scores:
                weighted_scores[pred_type] = 0.0
            weighted_scores[pred_type] += pred.confidence * pred.score
        
        if weighted_scores:
            return max(weighted_scores.items(), key=lambda x: x[1])[0]
        else:
            return "unknown"
    
    def _generate_reasoning(self, predictions: List[ModelPrediction], consensus_type: str, agreement_score: float) -> str:
        """Generate reasoning for consensus"""
        if agreement_score >= 0.8:
            return f"Strong consensus ({agreement_score:.1%}) among models for {consensus_type}"
        elif agreement_score >= 0.6:
            return f"Moderate consensus ({agreement_score:.1%}) among models for {consensus_type}"
        else:
            return f"Weak consensus ({agreement_score:.1%}), using weighted analysis"
    
    def _empty_conflict_analysis(self) -> Dict[str, Any]:
        """Return empty conflict analysis"""
        return {
            "predictions": [],
            "prediction_groups": {},
            "agreement_metrics": {
                "total_predictions": 0,
                "group_distribution": {},
                "agreement_score": 0.0,
                "consensus_type": "none",
                "confidence_range": {"min": 0.0, "max": 0.0, "mean": 0.0}
            },
            "conflicts": [],
            "consensus": {
                "recommendation": "unknown",
                "confidence": "low",
                "agreement_score": 0.0,
                "reasoning": "No predictions available"
            },
            "analysis_timestamp": pd.Timestamp.now().isoformat()
        }

class ConflictResolver:
    """Main conflict resolution orchestrator"""
    
    def __init__(self):
        self.logger = structlog.get_logger(__name__)
        self.analyzer = ConflictAnalyzer()
    
    async def resolve_model_conflicts(
        self,
        output_dir: str = "conflict_resolution"
    ) -> Dict[str, Any]:
        """Resolve conflicts between model predictions"""
        
        self.logger.info("🔍 Starting conflict resolution analysis")
        
        try:
            # Load predictions from different components
            predictions = await self._load_predictions()
            
            if not predictions:
                self.logger.warning("No predictions found for conflict resolution")
                return self._empty_resolution_result()
            
            # Analyze conflicts
            self.logger.info("🔄 Analyzing model conflicts...")
            conflict_analysis = self.analyzer.analyze_conflicts(predictions)
            
            # Generate resolution summary
            resolution_summary = self._generate_resolution_summary(conflict_analysis)
            
            # Save results
            self._save_results(conflict_analysis, resolution_summary, output_dir)
            
            self.logger.info(f"✅ Conflict resolution complete: {len(conflict_analysis['conflicts'])} conflicts identified")
            
            return resolution_summary
            
        except Exception as e:
            self.logger.error(f"❌ Conflict resolution failed: {e}")
            return self._empty_resolution_result()
    
    async def _load_predictions(self) -> List[ModelPrediction]:
        """Load predictions from different pipeline components"""
        predictions = []
        
        try:
            # Load IMPACT predictions
            impact_path = "impact_risk/impact_summary.json"
            if os.path.exists(impact_path):
                with open(impact_path, 'r') as f:
                    impact_data = json.load(f)
                
                if "decision_flag" in impact_data:
                    predictions.append(ModelPrediction(
                        "IMPACT",
                        {
                            "prediction": impact_data["decision_flag"],
                            "confidence": impact_data.get("confidence", 0.5),
                            "score": impact_data.get("selectivity_score", 0.0),
                            "metadata": {"component": "impact_risk"}
                        }
                    ))
            
            # Load SPARROW predictions
            sparrow_path = "sparrow/synthesis_analysis.json"
            if os.path.exists(sparrow_path):
                with open(sparrow_path, 'r') as f:
                    sparrow_data = json.load(f)
                
                if "summary_stats" in sparrow_data:
                    avg_score = sparrow_data["summary_stats"].get("avg_synthesis_score", 0.0)
                    if avg_score > 0.7:
                        prediction = "Synthesize"
                    elif avg_score > 0.4:
                        prediction = "Proceed with Caution"
                    else:
                        prediction = "Modify Required"
                    
                    predictions.append(ModelPrediction(
                        "SPARROW",
                        {
                            "prediction": prediction,
                            "confidence": avg_score,
                            "score": avg_score,
                            "metadata": {"component": "sparrow"}
                        }
                    ))
            
            # Load Expression Filter predictions
            expression_path = "expression_filter/tissue_weighted_risk.json"
            if os.path.exists(expression_path):
                with open(expression_path, 'r') as f:
                    expression_data = json.load(f)
                
                if "overall_risk" in expression_data:
                    risk_score = expression_data["overall_risk"]
                    if risk_score < 0.3:
                        prediction = "Synthesize"
                    elif risk_score < 0.6:
                        prediction = "Proceed with Caution"
                    else:
                        prediction = "Modify Required"
                    
                    predictions.append(ModelPrediction(
                        "Expression Filter",
                        {
                            "prediction": prediction,
                            "confidence": 1.0 - risk_score,
                            "score": 1.0 - risk_score,
                            "metadata": {"component": "expression_filter"}
                        }
                    ))
            
            # Load Structure Modeling predictions
            structure_path = "structure_modeling/binding_risk.json"
            if os.path.exists(structure_path):
                with open(structure_path, 'r') as f:
                    structure_data = json.load(f)
                
                if "binding_score" in structure_data:
                    binding_score = structure_data["binding_score"]
                    if binding_score < 0.3:
                        prediction = "Synthesize"
                    elif binding_score < 0.6:
                        prediction = "Proceed with Caution"
                    else:
                        prediction = "Modify Required"
                    
                    predictions.append(ModelPrediction(
                        "Structure Modeling",
                        {
                            "prediction": prediction,
                            "confidence": 1.0 - binding_score,
                            "score": 1.0 - binding_score,
                            "metadata": {"component": "structure_modeling"}
                        }
                    ))
            
        except Exception as e:
            self.logger.error(f"Failed to load predictions: {e}")
        
        return predictions
    
    def _generate_resolution_summary(self, conflict_analysis: Dict[str, Any]) -> Dict[str, Any]:
        """Generate resolution summary"""
        consensus = conflict_analysis["consensus"]
        conflicts = conflict_analysis["conflicts"]
        
        return {
            "final_recommendation": consensus["recommendation"],
            "confidence": consensus["confidence"],
            "agreement_score": consensus["agreement_score"],
            "num_conflicts": len(conflicts),
            "conflict_severity": self._assess_overall_conflict_severity(conflicts),
            "reasoning": consensus["reasoning"],
            "timestamp": pd.Timestamp.now().isoformat()
        }
    
    def _assess_overall_conflict_severity(self, conflicts: List[Dict[str, Any]]) -> str:
        """Assess overall conflict severity"""
        if not conflicts:
            return "none"
        
        severities = [conflict["severity"] for conflict in conflicts]
        
        if "high" in severities:
            return "high"
        elif "medium" in severities:
            return "medium"
        else:
            return "low"
    
    def _save_results(self, conflict_analysis: Dict[str, Any], resolution_summary: Dict[str, Any], output_dir: str):
        """Save conflict resolution results"""
        try:
            os.makedirs(output_dir, exist_ok=True)
            
            # Convert all objects to JSON-serializable format
            serializable_analysis = {}
            for key, value in conflict_analysis.items():
                if key == "predictions":
                    serializable_analysis[key] = [
                        pred.to_dict() if hasattr(pred, 'to_dict') else str(pred) 
                        for pred in value
                    ]
                elif key == "prediction_groups":
                    serializable_analysis[key] = {
                        group_name: [pred.to_dict() if hasattr(pred, 'to_dict') else str(pred) for pred in group_preds]
                        for group_name, group_preds in value.items()
                    }
                else:
                    serializable_analysis[key] = value
            
            # Save detailed analysis
            analysis_path = os.path.join(output_dir, "conflict_analysis.json")
            with open(analysis_path, 'w') as f:
                json.dump(serializable_analysis, f, indent=2, default=str)
            
            # Save resolution summary
            summary_path = os.path.join(output_dir, "conflict_summary.json")
            with open(summary_path, 'w') as f:
                json.dump(resolution_summary, f, indent=2, default=str)
            
        except Exception as e:
            self.logger.error(f"Failed to save conflict resolution results: {e}")
    
    def _empty_resolution_result(self) -> Dict[str, Any]:
        """Return empty resolution result"""
        return {
            "final_recommendation": "unknown",
            "confidence": "low",
            "agreement_score": 0.0,
            "num_conflicts": 0,
            "conflict_severity": "none",
            "reasoning": "No predictions available for conflict resolution",
            "timestamp": pd.Timestamp.now().isoformat()
        }

# Global resolver instance
_resolver = None

async def resolve_model_conflicts(
    output_dir: str = "conflict_resolution"
) -> Dict[str, Any]:
    """Main function to resolve model conflicts"""
    global _resolver
    
    if _resolver is None:
        _resolver = ConflictResolver()
    
    return await _resolver.resolve_model_conflicts(output_dir=output_dir)
