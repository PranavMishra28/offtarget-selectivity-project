#!/usr/bin/env python3
"""
Simple test script for AI Explanations
Demonstrates how to test the AI explanation generator with different scenarios.
"""

import asyncio
import json
import os
from ai.explanation_generator import generate_ai_explanation

async def test_ai_explanations():
    """Test AI explanations with different scenarios"""
    
    print("🧠 Testing AI Explanation System")
    print("=" * 50)
    
    # Test scenarios
    test_scenarios = [
        {
            "name": "Synthesis Recommended",
            "data": {
                "results": {
                    "impact_risk": {
                        "decision_flag": "Synthesize",
                        "selectivity_score": 0.85,
                        "safety_score": 0.78,
                        "risky_offtargets": []
                    }
                }
            }
        },
        {
            "name": "Proceed with Caution",
            "data": {
                "results": {
                    "impact_risk": {
                        "decision_flag": "Watch",
                        "selectivity_score": 0.65,
                        "safety_score": 0.72,
                        "risky_offtargets": [
                            {"uniprot_id": "P08172", "combined_score": 0.45}
                        ]
                    }
                }
            }
        },
        {
            "name": "Modifications Required",
            "data": {
                "results": {
                    "impact_risk": {
                        "decision_flag": "Modify",
                        "selectivity_score": 0.45,
                        "safety_score": 0.38,
                        "risky_offtargets": [
                            {"uniprot_id": "P08172", "combined_score": 0.78},
                            {"uniprot_id": "P08173", "combined_score": 0.65}
                        ]
                    }
                }
            }
        },
        {
            "name": "Synthesis Not Recommended",
            "data": {
                "results": {
                    "impact_risk": {
                        "decision_flag": "Reject",
                        "selectivity_score": 0.25,
                        "safety_score": 0.15,
                        "risky_offtargets": [
                            {"uniprot_id": "P08172", "combined_score": 0.92},
                            {"uniprot_id": "P08173", "combined_score": 0.88},
                            {"uniprot_id": "P08174", "combined_score": 0.85}
                        ]
                    }
                }
            }
        }
    ]
    
    # Create output directory
    os.makedirs("test_outputs/ai", exist_ok=True)
    
    for i, scenario in enumerate(test_scenarios):
        print(f"\n📋 Test {i+1}: {scenario['name']}")
        print("-" * 30)
        
        try:
            # Generate AI explanation
            result = await generate_ai_explanation(
                results=scenario['data'],
                output_path=f"test_outputs/ai/explanation_{i+1}.json"
            )
            
            # Display key information
            if "synthesis_explanation" in result:
                explanation = result["synthesis_explanation"]
                decision = explanation.get("decision", {})
                
                print(f"✅ Decision: {decision.get('title', 'Unknown')}")
                print(f"📝 Summary: {decision.get('summary', 'No summary')}")
                print(f"🎯 Confidence: {decision.get('confidence', 'Unknown')}")
                
                # Show key factors
                key_factors = explanation.get("key_factors", [])
                if key_factors:
                    print("🔍 Key Factors:")
                    for factor in key_factors[:3]:  # Show first 3
                        print(f"   • {factor}")
                
                # Show next steps
                next_steps = explanation.get("next_steps", [])
                if next_steps:
                    print("📋 Next Steps:")
                    for step in next_steps[:3]:  # Show first 3
                        print(f"   • {step}")
                
                print(f"💾 Saved to: test_outputs/ai/explanation_{i+1}.json")
                
            else:
                print("❌ No synthesis explanation found in result")
                
        except Exception as e:
            print(f"❌ Error testing scenario {i+1}: {e}")
    
    print(f"\n🎉 AI Explanation Testing Complete!")
    print(f"📁 Check the 'test_outputs/ai/' directory for detailed JSON reports")

if __name__ == "__main__":
    asyncio.run(test_ai_explanations()) 