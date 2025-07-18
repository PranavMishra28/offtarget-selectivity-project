"""
AI-Powered Explanation Generator
Provides human-readable insights and rationale for synthesis decisions using LLM techniques.
"""

import os
import json
from typing import Dict, Any, List, Optional
import pandas as pd
import numpy as np
import structlog
from datetime import datetime
import re

class ExplanationGenerator:
    """AI-powered explanation generator for synthesis decisions"""
    
    def __init__(self):
        self.logger = structlog.get_logger(__name__)
        self.templates = self._load_explanation_templates()
        
    def _load_explanation_templates(self) -> Dict[str, Dict[str, Any]]:
        """Load explanation templates for different scenarios"""
        return {
            "synthesis_decision": {
                "synthesize": {
                    "title": "✅ Synthesis Recommended",
                    "summary": "This compound demonstrates excellent selectivity and safety profile, making it suitable for synthesis.",
                    "key_factors": [
                        "High selectivity score indicates minimal off-target binding",
                        "Favorable safety profile with low toxicity risk",
                        "Good synthesis feasibility score",
                        "Acceptable expression-weighted risk assessment"
                    ],
                    "next_steps": [
                        "Proceed with synthesis planning",
                        "Initiate preclinical safety studies",
                        "Prepare for IND-enabling studies"
                    ]
                },
                "watch": {
                    "title": "⚠️ Proceed with Caution",
                    "summary": "This compound shows moderate risk factors that require careful monitoring during development.",
                    "key_factors": [
                        "Moderate selectivity concerns identified",
                        "Some off-target binding predicted",
                        "Requires enhanced safety monitoring",
                        "Consider structural modifications"
                    ],
                    "next_steps": [
                        "Implement enhanced safety monitoring",
                        "Consider structural modifications",
                        "Perform additional selectivity studies"
                    ]
                },
                "modify": {
                    "title": "🔧 Structural Modifications Required",
                    "summary": "This compound requires structural modifications to improve selectivity or safety before synthesis.",
                    "key_factors": [
                        "Significant off-target binding predicted",
                        "Safety concerns identified",
                        "Synthesis feasibility may be challenging",
                        "Expression-weighted risks are elevated"
                    ],
                    "next_steps": [
                        "Design structural modifications",
                        "Focus on reducing off-target binding",
                        "Improve safety profile",
                        "Re-evaluate after modifications"
                    ]
                },
                "reject": {
                    "title": "❌ Synthesis Not Recommended",
                    "summary": "This compound has significant safety or selectivity issues that preclude synthesis.",
                    "key_factors": [
                        "High off-target binding risk",
                        "Significant safety concerns",
                        "Poor synthesis feasibility",
                        "Unacceptable expression-weighted risks"
                    ],
                    "next_steps": [
                        "Consider alternative scaffolds",
                        "Explore different chemical series",
                        "Re-evaluate target selection"
                    ]
                }
            },
            "risk_assessment": {
                "high_risk": "High risk compounds show significant off-target binding (>0.7) to multiple targets, particularly those expressed in critical tissues.",
                "medium_risk": "Medium risk compounds show moderate off-target binding (0.4-0.7) to a limited number of targets.",
                "low_risk": "Low risk compounds show minimal off-target binding (<0.4) and good selectivity profile."
            },
            "selectivity_analysis": {
                "excellent": "Excellent selectivity (>0.8) indicates strong binding to the primary target with minimal off-target interactions.",
                "good": "Good selectivity (0.6-0.8) shows preferential binding to the primary target with some off-target concerns.",
                "moderate": "Moderate selectivity (0.4-0.6) indicates potential off-target binding that requires attention.",
                "poor": "Poor selectivity (<0.4) shows significant off-target binding that may limit therapeutic utility."
            }
        }
    
    def generate_synthesis_explanation(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Generate comprehensive synthesis decision explanation"""
        impact_data = results.get("results", {}).get("impact_risk", {})
        decision = impact_data.get("decision_flag", "Unknown")
        
        # Get template for decision
        template = self.templates["synthesis_decision"].get(decision.lower(), 
                                                          self.templates["synthesis_decision"]["watch"])
        
        # Extract key metrics
        selectivity_score = impact_data.get("selectivity_score", 0)
        safety_score = impact_data.get("safety_score", 0)
        risky_targets = impact_data.get("risky_offtargets", [])
        
        # Generate detailed explanation
        explanation = {
            "decision": {
                "title": template["title"],
                "summary": template["summary"],
                "confidence": self._calculate_confidence(selectivity_score, safety_score, len(risky_targets))
            },
            "key_factors": template["key_factors"],
            "next_steps": template["next_steps"],
            "detailed_analysis": self._generate_detailed_analysis(results),
            "risk_breakdown": self._generate_risk_breakdown(results),
            "recommendations": self._generate_specific_recommendations(results),
            "metadata": {
                "generation_timestamp": pd.Timestamp.now().isoformat(),
                "version": "1.0"
            }
        }
        
        return explanation
    
    def _calculate_confidence(self, selectivity: float, safety: float, risky_targets: int) -> str:
        """Calculate confidence level in the decision"""
        avg_score = (selectivity + safety) / 2
        
        if avg_score > 0.8 and risky_targets <= 2:
            return "High"
        elif avg_score > 0.6 and risky_targets <= 5:
            return "Medium"
        else:
            return "Low"
    
    def _generate_detailed_analysis(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Generate detailed analysis of each component"""
        analysis = {}
        
        # Impact Risk Analysis
        impact_data = results.get("results", {}).get("impact_risk", {})
        if impact_data:
            analysis["impact_risk"] = {
                "selectivity_score": impact_data.get("selectivity_score", 0),
                "safety_score": impact_data.get("safety_score", 0),
                "risky_offtargets_count": len(impact_data.get("risky_offtargets", [])),
                "top_risky_targets": [
                    {"uniprot_id": t["uniprot_id"], "score": t["combined_score"]}
                    for t in impact_data.get("risky_offtargets", [])[:5]
                ]
            }
        
        # Expression Analysis
        expression_data = results.get("results", {}).get("expression_filter", {})
        if expression_data:
            weighted_risks = expression_data.get("weighted_risks", {})
            analysis["expression_analysis"] = {
                "tissues_analyzed": len(weighted_risks),
                "high_risk_tissues": [
                    tissue for tissue, score in weighted_risks.items() 
                    if score > 0.7
                ],
                "average_tissue_risk": np.mean(list(weighted_risks.values())) if weighted_risks else 0
            }
        
        # Synthesis Analysis
        sparrow_data = results.get("results", {}).get("sparrow", {})
        if sparrow_data:
            analysis["synthesis_analysis"] = {
                "compounds_analyzed": len(sparrow_data.get("detailed_analysis", {})),
                "feasible_compounds": len([
                    comp for comp, data in sparrow_data.get("detailed_analysis", {}).items()
                    if data.get("synthesis_analysis", {}).get("feasibility_score", 0) > 0.6
                ]),
                "average_feasibility": np.mean([
                    data.get("synthesis_analysis", {}).get("feasibility_score", 0)
                    for data in sparrow_data.get("detailed_analysis", {}).values()
                ]) if sparrow_data.get("detailed_analysis") else 0
            }
        
        return analysis
    
    def _generate_risk_breakdown(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Generate detailed risk breakdown"""
        impact_data = results.get("results", {}).get("impact_risk", {})
        risky_targets = impact_data.get("risky_offtargets", [])
        
        # Categorize risks by target class
        risk_categories = {
            "ion_channels": [],
            "enzymes": [],
            "receptors": [],
            "transporters": [],
            "other": []
        }
        
        for target in risky_targets:
            uniprot_id = target["uniprot_id"]
            score = target["combined_score"]
            
            # Simple categorization based on common patterns
            if any(pattern in uniprot_id.lower() for pattern in ["kcn", "scn", "cac"]):
                risk_categories["ion_channels"].append({"id": uniprot_id, "score": score})
            elif any(pattern in uniprot_id.lower() for pattern in ["kinase", "ase"]):
                risk_categories["enzymes"].append({"id": uniprot_id, "score": score})
            elif any(pattern in uniprot_id.lower() for pattern in ["gpr", "receptor"]):
                risk_categories["receptors"].append({"id": uniprot_id, "score": score})
            elif any(pattern in uniprot_id.lower() for pattern in ["slc", "abc"]):
                risk_categories["transporters"].append({"id": uniprot_id, "score": score})
            else:
                risk_categories["other"].append({"id": uniprot_id, "score": score})
        
        # Calculate category risks
        category_risks = {}
        for category, targets in risk_categories.items():
            if targets:
                category_risks[category] = {
                    "count": len(targets),
                    "max_score": max(t["score"] for t in targets),
                    "avg_score": np.mean([t["score"] for t in targets]),
                    "targets": targets
                }
            else:
                category_risks[category] = {
                    "count": 0,
                    "max_score": 0,
                    "avg_score": 0,
                    "targets": []
                }
        
        return category_risks
    
    def _generate_specific_recommendations(self, results: Dict[str, Any]) -> List[str]:
        """Generate specific, actionable recommendations"""
        recommendations = []
        
        impact_data = results.get("results", {}).get("impact_risk", {})
        selectivity_score = impact_data.get("selectivity_score", 0)
        safety_score = impact_data.get("safety_score", 0)
        risky_targets = impact_data.get("risky_offtargets", [])
        
        # Selectivity-based recommendations
        if selectivity_score < 0.5:
            recommendations.append("🔍 Perform additional selectivity profiling against key off-targets")
            recommendations.append("🧬 Consider structural modifications to improve target specificity")
        
        if selectivity_score < 0.3:
            recommendations.append("⚠️ Significant off-target binding predicted - consider alternative scaffolds")
        
        # Safety-based recommendations
        if safety_score < 0.6:
            recommendations.append("🛡️ Implement enhanced safety monitoring in preclinical studies")
            recommendations.append("🧪 Perform comprehensive toxicology studies")
        
        if safety_score < 0.4:
            recommendations.append("🚨 High safety risk - consider structural modifications or alternative approaches")
        
        # Off-target specific recommendations
        ion_channel_risks = [t for t in risky_targets if any(pattern in t["uniprot_id"].lower() 
                                                           for pattern in ["kcn", "scn", "cac"])]
        if ion_channel_risks:
            recommendations.append("💓 Monitor for cardiac effects - ion channel interactions detected")
        
        enzyme_risks = [t for t in risky_targets if any(pattern in t["uniprot_id"].lower() 
                                                       for pattern in ["kinase", "ase"])]
        if enzyme_risks:
            recommendations.append("⚗️ Monitor for metabolic interactions - enzyme off-targets detected")
        
        # Expression-based recommendations
        expression_data = results.get("results", {}).get("expression_filter", {})
        weighted_risks = expression_data.get("weighted_risks", {})
        high_risk_tissues = [tissue for tissue, score in weighted_risks.items() if score > 0.7]
        
        if "brain" in high_risk_tissues:
            recommendations.append("🧠 Monitor for CNS effects - high expression in brain tissues")
        
        if "heart" in high_risk_tissues:
            recommendations.append("❤️ Enhanced cardiac monitoring required - high expression in heart tissues")
        
        if "liver" in high_risk_tissues:
            recommendations.append("🫁 Monitor liver function - high expression in liver tissues")
        
        return recommendations
    
    def generate_ai_insights(self, results: Dict[str, Any]) -> Dict[str, Any]:
        """Generate AI-powered insights and observations"""
        insights = {
            "key_observations": [],
            "unusual_patterns": [],
            "opportunities": [],
            "warnings": []
        }
        
        impact_data = results.get("results", {}).get("impact_risk", {})
        selectivity_score = impact_data.get("selectivity_score", 0)
        safety_score = impact_data.get("safety_score", 0)
        risky_targets = impact_data.get("risky_offtargets", [])
        
        # Key observations
        if selectivity_score > 0.8:
            insights["key_observations"].append("🎯 Exceptional selectivity profile - minimal off-target binding predicted")
        
        if safety_score > 0.8:
            insights["key_observations"].append("🛡️ Excellent safety profile - low toxicity risk predicted")
        
        if len(risky_targets) <= 2:
            insights["key_observations"].append("✅ Limited off-target interactions - focused target profile")
        
        # Unusual patterns
        if selectivity_score > 0.9 and safety_score < 0.5:
            insights["unusual_patterns"].append("⚠️ Unusual pattern: High selectivity but low safety - investigate further")
        
        if len(risky_targets) > 10:
            insights["unusual_patterns"].append("🔍 High promiscuity detected - broad off-target binding profile")
        
        # Opportunities
        if selectivity_score > 0.7 and safety_score > 0.7:
            insights["opportunities"].append("🚀 Strong candidate for rapid development - excellent overall profile")
        
        if len(risky_targets) == 0:
            insights["opportunities"].append("💎 Rare finding: No significant off-targets detected - highly selective compound")
        
        # Warnings
        if selectivity_score < 0.3:
            insights["warnings"].append("🚨 Critical: Very low selectivity - high risk of off-target effects")
        
        if safety_score < 0.3:
            insights["warnings"].append("🚨 Critical: Very low safety score - high toxicity risk")
        
        if len(risky_targets) > 15:
            insights["warnings"].append("🚨 Critical: Excessive off-target binding - consider alternative scaffolds")
        
        return insights

async def generate_ai_explanation(
    results: Dict[str, Any],
    output_path: str = "ai/explanation_report.json"
) -> Dict[str, Any]:
    """
    Generate AI-powered explanation for synthesis decisions.
    
    Args:
        results: Pipeline results dictionary
        output_path: Path for output file
    
    Returns:
        Comprehensive AI explanation report
    """
    
    logger = structlog.get_logger(__name__)
    logger.info("Generating AI-powered explanation")
    
    generator = ExplanationGenerator()
    
    # Generate explanations
    synthesis_explanation = generator.generate_synthesis_explanation(results)
    ai_insights = generator.generate_ai_insights(results)
    
    # Combine into comprehensive report
    report = {
        "synthesis_explanation": synthesis_explanation,
        "ai_insights": ai_insights,
        "metadata": {
            "generation_timestamp": pd.Timestamp.now().isoformat(),
            "version": "1.0"
        }
    }
    
    # Save report
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(report, f, indent=2)
    
    logger.info(f"✅ AI explanation generated")
    logger.info(f"📁 Report: {output_path}")
    
    return report 