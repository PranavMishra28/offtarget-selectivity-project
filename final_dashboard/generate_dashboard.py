"""
Enhanced Final Dashboard - Comprehensive Visualization & Analytics
Implements advanced dashboard generation with interactive visualizations, automated reporting, and comprehensive analytics.
"""

import os
import json
import asyncio
from typing import Dict, Any, List, Optional, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import structlog
from tqdm import tqdm
from pathlib import Path

# Import configuration
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.config_manager import config_manager

class DashboardDataAggregator:
    """Aggregates data from all pipeline components"""
    
    def __init__(self):
        self.logger = structlog.get_logger(__name__)
    
    def aggregate_pipeline_data(self, base_path: str = ".") -> Dict[str, Any]:
        """Aggregate data from all pipeline components"""
        
        aggregated_data = {
            "nebula": self._load_nebula_data(base_path),
            "sparrow": self._load_sparrow_data(base_path),
            "empirical_binding": self._load_empirical_data(base_path),
            "structure_modeling": self._load_structure_data(base_path),
            "impact_risk": self._load_impact_data(base_path),
            "expression_filter": self._load_expression_data(base_path),
            "conflict_resolution": self._load_conflict_data(base_path)
        }
        
        # Add metadata
        aggregated_data["metadata"] = {
            "aggregation_timestamp": pd.Timestamp.now().isoformat(),
            "pipeline_version": "2.0",
            "data_sources": list(aggregated_data.keys())
        }
        
        return aggregated_data
    
    def _load_nebula_data(self, base_path: str) -> Dict[str, Any]:
        """Load NEBULA generation data"""
        try:
            metadata_path = os.path.join(base_path, "nebula/generation_metadata.json")
            if os.path.exists(metadata_path):
                with open(metadata_path, 'r') as f:
                    return json.load(f)
            return {"status": "not_found"}
        except Exception as e:
            self.logger.error(f"Failed to load NEBULA data: {e}")
            return {"status": "error", "error": str(e)}
    
    def _load_sparrow_data(self, base_path: str) -> Dict[str, Any]:
        """Load SPARROW synthesis data"""
        try:
            csv_path = os.path.join(base_path, "sparrow/ranked_candidates.csv")
            metadata_path = os.path.join(base_path, "sparrow/synthesis_analysis.json")
            
            data = {"status": "loaded"}
            
            if os.path.exists(csv_path):
                data["ranked_candidates"] = pd.read_csv(csv_path).to_dict('records')
            
            if os.path.exists(metadata_path):
                with open(metadata_path, 'r') as f:
                    data["synthesis_analysis"] = json.load(f)
            
            return data
        except Exception as e:
            self.logger.error(f"Failed to load SPARROW data: {e}")
            return {"status": "error", "error": str(e)}
    
    def _load_empirical_data(self, base_path: str) -> Dict[str, Any]:
        """Load empirical binding data"""
        try:
            predictions_path = os.path.join(base_path, "empirical_binding/offtarget_predictions.json")
            metadata_path = os.path.join(base_path, "empirical_binding/prediction_metadata.json")
            
            data = {"status": "loaded"}
            
            if os.path.exists(predictions_path):
                with open(predictions_path, 'r') as f:
                    data["predictions"] = json.load(f)
            
            if os.path.exists(metadata_path):
                with open(metadata_path, 'r') as f:
                    data["metadata"] = json.load(f)
            
            return data
        except Exception as e:
            self.logger.error(f"Failed to load empirical data: {e}")
            return {"status": "error", "error": str(e)}
    
    def _load_structure_data(self, base_path: str) -> Dict[str, Any]:
        """Load structure modeling data"""
        try:
            risk_path = os.path.join(base_path, "structure_modeling/binding_risk.json")
            data = {"status": "loaded"}
            
            if os.path.exists(risk_path):
                with open(risk_path, 'r') as f:
                    data["binding_risk"] = json.load(f)
            
            return data
        except Exception as e:
            self.logger.error(f"Failed to load structure data: {e}")
            return {"status": "error", "error": str(e)}
    
    def _load_impact_data(self, base_path: str) -> Dict[str, Any]:
        """Load IMPACT risk data"""
        try:
            summary_path = os.path.join(base_path, "impact_risk/impact_summary.json")
            detailed_path = os.path.join(base_path, "impact_risk/detailed_analysis.json")
            
            data = {"status": "loaded"}
            
            if os.path.exists(summary_path):
                with open(summary_path, 'r') as f:
                    data["summary"] = json.load(f)
            
            if os.path.exists(detailed_path):
                with open(detailed_path, 'r') as f:
                    data["detailed"] = json.load(f)
            
            return data
        except Exception as e:
            self.logger.error(f"Failed to load IMPACT data: {e}")
            return {"status": "error", "error": str(e)}
    
    def _load_expression_data(self, base_path: str) -> Dict[str, Any]:
        """Load expression filter data"""
        try:
            risk_path = os.path.join(base_path, "expression_filter/tissue_weighted_risk.json")
            analysis_path = os.path.join(base_path, "expression_filter/expression_analysis.json")
            
            data = {"status": "loaded"}
            
            if os.path.exists(risk_path):
                with open(risk_path, 'r') as f:
                    data["weighted_risks"] = json.load(f)
            
            if os.path.exists(analysis_path):
                with open(analysis_path, 'r') as f:
                    data["analysis"] = json.load(f)
            
            return data
        except Exception as e:
            self.logger.error(f"Failed to load expression data: {e}")
            return {"status": "error", "error": str(e)}
    
    def _load_conflict_data(self, base_path: str) -> Dict[str, Any]:
        """Load conflict resolution data"""
        try:
            conflict_path = os.path.join(base_path, "conflict_resolution/conflict_summary.json")
            data = {"status": "loaded"}
            
            if os.path.exists(conflict_path):
                with open(conflict_path, 'r') as f:
                    data["conflicts"] = json.load(f)
            
            return data
        except Exception as e:
            self.logger.error(f"Failed to load conflict data: {e}")
            return {"status": "error", "error": str(e)}

class DashboardVisualizer:
    """Advanced dashboard visualization generation"""
    
    def __init__(self):
        self.logger = structlog.get_logger(__name__)
        self.colors = {
            'primary': '#1f77b4',
            'secondary': '#ff7f0e',
            'success': '#2ca02c',
            'warning': '#d62728',
            'info': '#9467bd'
        }
    
    def generate_comprehensive_dashboard(
        self, 
        aggregated_data: Dict[str, Any], 
        output_dir: str = "final_dashboard"
    ) -> Dict[str, str]:
        """Generate comprehensive dashboard with all visualizations"""
        os.makedirs(output_dir, exist_ok=True)
        visualization_paths = {}
        # 1. Pipeline Overview
        self.logger.info("Generating pipeline overview...")
        overview_path = os.path.join(output_dir, "pipeline_overview.html")
        self._generate_pipeline_overview(aggregated_data, overview_path)
        visualization_paths["pipeline_overview"] = overview_path
        # 2. Risk Assessment Dashboard
        self.logger.info("Generating risk assessment dashboard...")
        risk_path = os.path.join(output_dir, "risk_assessment.html")
        self._generate_risk_dashboard(aggregated_data, risk_path)
        visualization_paths["risk_assessment"] = risk_path
        # 3. Selectivity Analysis
        self.logger.info("Generating selectivity analysis...")
        selectivity_path = os.path.join(output_dir, "selectivity_analysis.html")
        self._generate_selectivity_analysis(aggregated_data, selectivity_path)
        visualization_paths["selectivity_analysis"] = selectivity_path
        # 4. Expression Analysis
        self.logger.info("Generating expression analysis...")
        expression_path = os.path.join(output_dir, "expression_analysis.html")
        self._generate_expression_analysis(aggregated_data, expression_path)
        visualization_paths["expression_analysis"] = expression_path
        # 5. Synthesis Feasibility
        self.logger.info("Generating synthesis feasibility analysis...")
        synthesis_path = os.path.join(output_dir, "synthesis_feasibility.html")
        self._generate_synthesis_analysis(aggregated_data, synthesis_path)
        visualization_paths["synthesis_feasibility"] = synthesis_path
        # 6. Static Summary Plot
        self.logger.info("Generating static summary plot...")
        summary_plot_path = os.path.join(output_dir, "summary_dashboard.png")
        self._generate_summary_plot(aggregated_data, summary_plot_path)
        visualization_paths["summary_plot"] = summary_plot_path
        # 7. Dashboard Index
        index_path = self.generate_dashboard_index(output_dir, visualization_paths)
        visualization_paths["dashboard_index"] = index_path
        return visualization_paths
    
    def generate_dashboard_index(self, output_dir: str, visualization_paths: Dict[str, str]) -> str:
        """Generate an index HTML page linking to all dashboard outputs"""
        index_path = os.path.join(output_dir, "dashboard_index.html")
        with open(index_path, 'w') as f:
            f.write("<html><head><title>Off-Target & Selectivity Dashboard</title></head><body>")
            f.write("<h1>Off-Target & Selectivity Pipeline Dashboard</h1>")
            f.write("<ul>")
            for name, path in visualization_paths.items():
                label = name.replace('_', ' ').title()
                rel_path = os.path.basename(path)
                f.write(f'<li><a href="{rel_path}">{label}</a></li>')
            f.write("</ul>")
            f.write("<p>Generated by Off-Target & Selectivity Pipeline v2.0</p>")
            f.write("</body></html>")
        return index_path
    
    def _generate_pipeline_overview(self, data: Dict[str, Any], output_path: str):
        """Generate pipeline overview dashboard"""
        try:
            # Create pipeline status overview
            fig = make_subplots(
                rows=2, cols=2,
                subplot_titles=('Pipeline Status', 'Data Flow', 'Key Metrics', 'Timeline'),
                specs=[[{"type": "indicator"}, {"type": "sankey"}],
                       [{"type": "bar"}, {"type": "scatter"}]]
            )
            
            # Pipeline status indicators
            pipeline_status = {
                "NEBULA": data.get("nebula", {}).get("status", "unknown"),
                "SPARROW": data.get("sparrow", {}).get("status", "unknown"),
                "Empirical": data.get("empirical_binding", {}).get("status", "unknown"),
                "Structure": data.get("structure_modeling", {}).get("status", "unknown"),
                "IMPACT": data.get("impact_risk", {}).get("status", "unknown"),
                "Expression": data.get("expression_filter", {}).get("status", "unknown")
            }
            
            # Status indicator
            fig.add_trace(
                go.Indicator(
                    mode="gauge+number+delta",
                    value=sum(1 for status in pipeline_status.values() if status == "loaded"),
                    domain={'x': [0, 1], 'y': [0, 1]},
                    title={'text': "Pipeline Completion"},
                    delta={'reference': 6},
                    gauge={'axis': {'range': [None, 6]},
                           'bar': {'color': "darkblue"},
                           'steps': [{'range': [0, 3], 'color': "lightgray"},
                                    {'range': [3, 5], 'color': "yellow"},
                                    {'range': [5, 6], 'color': "green"}]}
                ),
                row=1, col=1
            )
            
            # Key metrics
            impact_data = data.get("impact_risk", {}).get("summary", {})
            if impact_data:
                metrics = [
                    impact_data.get("selectivity_score", 0),
                    impact_data.get("safety_score", 0),
                    1 - impact_data.get("ion_channel_risks", {}).get("overall_ion_channel_risk", 0)
                ]
                
                fig.add_trace(
                    go.Bar(
                        x=['Selectivity', 'Safety', 'Ion Channel Safety'],
                        y=metrics,
                        marker_color=['#1f77b4', '#2ca02c', '#ff7f0e']
                    ),
                    row=2, col=1
                )
            
            fig.update_layout(height=800, title_text="Pipeline Overview Dashboard")
            fig.write_html(output_path)
            
        except Exception as e:
            self.logger.error(f"Failed to generate pipeline overview: {e}")
    
    def _generate_risk_dashboard(self, data: Dict[str, Any], output_path: str):
        """Generate risk assessment dashboard"""
        try:
            impact_data = data.get("impact_risk", {}).get("summary", {})
            if not impact_data:
                return
            
            fig = make_subplots(
                rows=2, cols=2,
                subplot_titles=('Risk Distribution', 'Ion Channel Risks', 'Safety Profile', 'Decision Factors'),
                specs=[[{"type": "histogram"}, {"type": "bar"}],
                       [{"type": "scatter"}, {"type": "scatter"}]]
            )
            
            # Risk distribution
            risky_targets = impact_data.get("risky_offtargets", [])
            if risky_targets:
                scores = [target["combined_score"] for target in risky_targets]
                fig.add_trace(
                    go.Histogram(x=scores, nbinsx=20, name="Risk Distribution"),
                    row=1, col=1
                )
            
            # Ion channel risks
            ion_risks = impact_data.get("ion_channel_risks", {})
            if ion_risks:
                risk_types = ['hERG', 'Nav', 'Ca2+']
                risk_values = [
                    ion_risks.get("herg_risk", 0),
                    ion_risks.get("nav_risk", 0),
                    ion_risks.get("ca_risk", 0)
                ]
                
                fig.add_trace(
                    go.Bar(x=risk_types, y=risk_values, marker_color='red'),
                    row=1, col=2
                )
            
            # Safety radar chart (converted to scatter)
            safety_metrics = [
                impact_data.get("safety_score", 0),
                1 - ion_risks.get("overall_ion_channel_risk", 0),
                impact_data.get("selectivity_score", 0),
                0.8  # Mock metabolic safety
            ]
            
            fig.add_trace(
                go.Scatter(
                    x=['Overall Safety', 'Ion Channel Safety', 'Selectivity', 'Metabolic Safety'],
                    y=safety_metrics,
                    mode='lines+markers',
                    fill='tonexty',
                    name='Safety Profile',
                    line=dict(color='blue', width=2),
                    marker=dict(size=8)
                ),
                row=2, col=1
            )
            
            fig.update_layout(height=800, title_text="Risk Assessment Dashboard")
            fig.write_html(output_path)
            
        except Exception as e:
            self.logger.error(f"Failed to generate risk dashboard: {e}")
    
    def _generate_selectivity_analysis(self, data: Dict[str, Any], output_path: str):
        """Generate selectivity analysis dashboard"""
        try:
            impact_data = data.get("impact_risk", {}).get("summary", {})
            if not impact_data:
                return
            
            fig = make_subplots(
                rows=2, cols=2,
                subplot_titles=('Selectivity Score', 'Off-Target Distribution', 'Target Comparison', 'Selectivity Index'),
                specs=[[{"type": "indicator"}, {"type": "pie"}],
                       [{"type": "bar"}, {"type": "scatter"}]]
            )
            
            # Selectivity score indicator
            selectivity_score = impact_data.get("selectivity_score", 0)
            fig.add_trace(
                go.Indicator(
                    mode="gauge+number",
                    value=selectivity_score,
                    domain={'x': [0, 1], 'y': [0, 1]},
                    title={'text': "Selectivity Score"},
                    gauge={'axis': {'range': [None, 1]},
                           'bar': {'color': "darkblue"},
                           'steps': [{'range': [0, 0.3], 'color': "red"},
                                    {'range': [0.3, 0.7], 'color': "yellow"},
                                    {'range': [0.7, 1], 'color': "green"}]}
                ),
                row=1, col=1
            )
            
            # Off-target distribution
            risky_targets = impact_data.get("risky_offtargets", [])
            if risky_targets:
                target_names = [target["uniprot_id"] for target in risky_targets[:10]]
                target_scores = [target["combined_score"] for target in risky_targets[:10]]
                
                fig.add_trace(
                    go.Bar(x=target_names, y=target_scores, marker_color='orange'),
                    row=2, col=1
                )
            
            fig.update_layout(height=800, title_text="Selectivity Analysis Dashboard")
            fig.write_html(output_path)
            
        except Exception as e:
            self.logger.error(f"Failed to generate selectivity analysis: {e}")
    
    def _generate_expression_analysis(self, data: Dict[str, Any], output_path: str):
        """Generate expression analysis dashboard"""
        try:
            expression_data = data.get("expression_filter", {})
            if not expression_data:
                return
            
            fig = make_subplots(
                rows=2, cols=2,
                subplot_titles=('Expression Weighting', 'Tissue Distribution', 'Risk Comparison', 'Expression Heatmap'),
                specs=[[{"type": "scatter"}, {"type": "bar"}],
                       [{"type": "scatter"}, {"type": "heatmap"}]]
            )
            
            # Expression weighting scatter
            weighted_risks = expression_data.get("weighted_risks", {})
            if weighted_risks:
                original_scores = []
                weighted_scores = []
                
                for uniprot_id, score in weighted_risks.items():
                    # Get original score from impact data
                    impact_data = data.get("impact_risk", {}).get("summary", {})
                    risky_targets = impact_data.get("risky_offtargets", [])
                    original_score = next((t["combined_score"] for t in risky_targets if t["uniprot_id"] == uniprot_id), score)
                    
                    original_scores.append(original_score)
                    weighted_scores.append(score)
                
                fig.add_trace(
                    go.Scatter(x=original_scores, y=weighted_scores, mode='markers', name='Risk Scores'),
                    row=1, col=1
                )
            
            fig.update_layout(height=800, title_text="Expression Analysis Dashboard")
            fig.write_html(output_path)
            
        except Exception as e:
            self.logger.error(f"Failed to generate expression analysis: {e}")
    
    def _generate_synthesis_analysis(self, data: Dict[str, Any], output_path: str):
        """Generate synthesis feasibility analysis"""
        try:
            sparrow_data = data.get("sparrow", {})
            if not sparrow_data:
                return
            
            fig = make_subplots(
                rows=2, cols=2,
                subplot_titles=('Synthesis Scores', 'Complexity Distribution', 'Feasibility vs Risk', 'Recommendations'),
                specs=[[{"type": "histogram"}, {"type": "scatter"}],
                       [{"type": "bar"}, {"type": "pie"}]]
            )
            
            ranked_candidates = sparrow_data.get("ranked_candidates", [])
            if ranked_candidates:
                scores = [candidate.get("synthesis_score", 0) for candidate in ranked_candidates]
                fig.add_trace(
                    go.Histogram(x=scores, nbinsx=20, name="Synthesis Scores"),
                    row=1, col=1
                )
            
            fig.update_layout(height=800, title_text="Synthesis Feasibility Dashboard")
            fig.write_html(output_path)
            
        except Exception as e:
            self.logger.error(f"Failed to generate synthesis analysis: {e}")
    
    def _generate_summary_plot(self, data: Dict[str, Any], output_path: str):
        """Generate static summary plot"""
        try:
            # Normalize the output path
            output_path = os.path.normpath(output_path)
            
            fig, axes = plt.subplots(2, 2, figsize=(15, 12))
            
            # Extract key data
            impact_data = data.get("impact_risk", {}).get("summary", {})
            
            # 1. Risk distribution
            risky_targets = impact_data.get("risky_offtargets", [])
            if risky_targets:
                scores = [target["combined_score"] for target in risky_targets[:10]]
                target_ids = [target["uniprot_id"] for target in risky_targets[:10]]
                
                axes[0, 0].barh(target_ids, scores, color='skyblue')
                axes[0, 0].set_title('Top 10 Off-Target Risks')
                axes[0, 0].set_xlabel('Risk Score')
            
            # 2. Selectivity metrics
            selectivity_score = impact_data.get("selectivity_score", 0)
            safety_score = impact_data.get("safety_score", 0)
            
            metrics = ['Selectivity', 'Safety']
            values = [selectivity_score, safety_score]
            colors = ['green' if v > 0.5 else 'red' for v in values]
            
            axes[0, 1].bar(metrics, values, color=colors)
            axes[0, 1].set_title('Key Metrics')
            axes[0, 1].set_ylim(0, 1)
            
            # 3. Decision breakdown
            decision = impact_data.get("decision_flag", "Unknown")
            decision_counts = {"Synthesize": 0, "Watch": 0, "Modify": 0, "Reject": 0}
            decision_counts[decision] = 1
            
            axes[1, 0].pie(decision_counts.values(), labels=decision_counts.keys(), autopct='%1.1f%%')
            axes[1, 0].set_title('Synthesis Decision')
            
            # 4. Ion channel risks
            ion_risks = impact_data.get("ion_channel_risks", {})
            if ion_risks:
                risk_types = ['hERG', 'Nav', 'Ca2+']
                risk_values = [
                    ion_risks.get("herg_risk", 0),
                    ion_risks.get("nav_risk", 0),
                    ion_risks.get("ca_risk", 0)
                ]
                
                axes[1, 1].bar(risk_types, risk_values, color=['red', 'orange', 'yellow'])
                axes[1, 1].set_title('Ion Channel Risks')
                axes[1, 1].set_ylim(0, 1)
            
            plt.tight_layout()
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            self.logger.error(f"Failed to generate summary plot: {e}")

class ReportGenerator:
    """Automated report generation"""
    
    def __init__(self):
        self.logger = structlog.get_logger(__name__)
    
    def generate_comprehensive_report(
        self, 
        aggregated_data: Dict[str, Any], 
        visualization_paths: Dict[str, str],
        output_path: str = "final_dashboard/comprehensive_report.html"
    ) -> str:
        """Generate comprehensive HTML report"""
        
        try:
            # Extract key information
            impact_data = aggregated_data.get("impact_risk", {}).get("summary", {})
            
            # Generate HTML report
            html_content = self._generate_html_report(aggregated_data, visualization_paths, impact_data)
            
            # Save report
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(html_content)
            
            self.logger.info(f"✅ Comprehensive report generated: {output_path}")
            return output_path
            
        except Exception as e:
            self.logger.error(f"Failed to generate comprehensive report: {e}")
            return ""
    
    def _generate_html_report(
        self, 
        aggregated_data: Dict[str, Any], 
        visualization_paths: Dict[str, str],
        impact_data: Dict[str, Any]
    ) -> str:
        """Generate HTML report content"""
        
        decision = impact_data.get("decision_flag", "Unknown")
        selectivity_score = impact_data.get("selectivity_score", 0)
        safety_score = impact_data.get("safety_score", 0)
        
        html = f"""
        <!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>Off-Target Selectivity Analysis Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; }}
                .section {{ margin: 20px 0; padding: 15px; border: 1px solid #ddd; border-radius: 5px; }}
                .metric {{ display: inline-block; margin: 10px; padding: 10px; background-color: #e8f4f8; border-radius: 5px; }}
                .decision {{ font-size: 24px; font-weight: bold; padding: 20px; border-radius: 5px; }}
                .decision.synthesize {{ background-color: #d4edda; color: #155724; }}
                .decision.watch {{ background-color: #fff3cd; color: #856404; }}
                .decision.reject {{ background-color: #f8d7da; color: #721c24; }}
                .visualization {{ margin: 20px 0; text-align: center; }}
                iframe {{ width: 100%; height: 600px; border: none; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>🧪 Off-Target Selectivity Analysis Report</h1>
                <p>Comprehensive analysis of compound selectivity, safety, and synthesis feasibility</p>
                <p><strong>Generated:</strong> {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            </div>
            
            <div class="section">
                <h2>🎯 Executive Summary</h2>
                <div class="decision {decision.lower()}">
                    Synthesis Decision: {decision}
                </div>
                <div class="metric">
                    <strong>Selectivity Score:</strong> {selectivity_score:.3f}
                </div>
                <div class="metric">
                    <strong>Safety Score:</strong> {safety_score:.3f}
                </div>
            </div>
            
            <div class="section">
                <h2>📊 Pipeline Overview</h2>
                <div class="visualization">
                    <iframe src="pipeline_overview.html"></iframe>
                </div>
            </div>
            
            <div class="section">
                <h2>⚠️ Risk Assessment</h2>
                <div class="visualization">
                    <iframe src="risk_assessment.html"></iframe>
                </div>
            </div>
            
            <div class="section">
                <h2>🎯 Selectivity Analysis</h2>
                <div class="visualization">
                    <iframe src="selectivity_analysis.html"></iframe>
                </div>
            </div>
            
            <div class="section">
                <h2>🧬 Expression Analysis</h2>
                <div class="visualization">
                    <iframe src="expression_analysis.html"></iframe>
                </div>
            </div>
            
            <div class="section">
                <h2>⚗️ Synthesis Feasibility</h2>
                <div class="visualization">
                    <iframe src="synthesis_feasibility.html"></iframe>
                </div>
            </div>
            
            <div class="section">
                <h2>📈 Summary Dashboard</h2>
                <div class="visualization">
                    <img src="summary_dashboard.png" alt="Summary Dashboard" style="max-width: 100%;">
                </div>
            </div>
            
            <div class="section">
                <h2>📋 Detailed Results</h2>
                <pre>{json.dumps(aggregated_data, indent=2, default=str)}</pre>
            </div>
        </body>
        </html>
        """
        
        return html

async def generate_final_dashboard(
    impact_path: str = "impact_risk/impact_summary.json",
    expression_path: str = "expression_filter/tissue_weighted_risk.json",
    output_dir: str = "final_dashboard"
) -> Dict[str, Any]:
    """
    Enhanced final dashboard generation with comprehensive analytics and visualization.
    
    Args:
        impact_path: Path to impact risk summary
        expression_path: Path to expression weighted risks
        output_dir: Output directory for dashboard
    
    Returns:
        Dictionary with dashboard paths and metadata
    """
    
    logger = structlog.get_logger(__name__)
    logger.info("Starting enhanced final dashboard generation")
    
    # Initialize components
    aggregator = DashboardDataAggregator()
    visualizer = DashboardVisualizer()
    report_generator = ReportGenerator()
    
    # Aggregate all pipeline data
    logger.info("Aggregating pipeline data...")
    aggregated_data = aggregator.aggregate_pipeline_data()
    
    # Generate comprehensive visualizations
    logger.info("Generating visualizations...")
    visualizations = visualizer.generate_comprehensive_dashboard(aggregated_data, output_dir)
    
    # Generate comprehensive report
    logger.info("Generating comprehensive report...")
    report_path = report_generator.generate_comprehensive_report(
        aggregated_data, visualizations, 
        os.path.join(output_dir, "comprehensive_report.html")
    )
    
    # Create legacy summary files for compatibility
    logger.info("Creating legacy summary files...")
    
    # Legacy JSON summary
    legacy_summary = {
        "compound_on_target": aggregated_data.get("impact_risk", {}).get("summary", {}).get("on_target_id", "Unknown"),
        "on_target_score": aggregated_data.get("impact_risk", {}).get("summary", {}).get("on_target_score", 0),
        "avg_off_target_score": aggregated_data.get("impact_risk", {}).get("summary", {}).get("avg_off_target_score", 0),
        "selectivity_score": aggregated_data.get("impact_risk", {}).get("summary", {}).get("selectivity_score", 0),
        "decision_flag": aggregated_data.get("impact_risk", {}).get("summary", {}).get("decision_flag", "Unknown"),
        "top_risky_offtargets": aggregated_data.get("impact_risk", {}).get("summary", {}).get("risky_offtargets", []),
        "expression_weighted_risks": aggregated_data.get("expression_filter", {}).get("weighted_risks", {}),
        "metadata": {
            "dashboard_version": "2.0",
            "generation_timestamp": pd.Timestamp.now().isoformat(),
            "visualization_paths": visualizations
        }
    }
    
    legacy_json_path = os.path.join(output_dir, "compound_summary.json")
    with open(legacy_json_path, 'w') as f:
        json.dump(legacy_summary, f, indent=2)
    
    # Legacy CSV summary
    if legacy_summary["expression_weighted_risks"]:
        flat_data = []
        for uid, score in legacy_summary["expression_weighted_risks"].items():
            flat_data.append({
                "uniprot_id": uid,
                "expression_weighted_score": score
            })
        
        df = pd.DataFrame(flat_data)
        legacy_csv_path = os.path.join(output_dir, "compound_summary.csv")
        df.to_csv(legacy_csv_path, index=False)
    
    # Compile results
    results = {
        "dashboard_dir": output_dir,
        "comprehensive_report": report_path,
        "visualizations": visualizations,
        "legacy_files": {
            "json_summary": legacy_json_path,
            "csv_summary": legacy_csv_path if 'legacy_csv_path' in locals() else None
        },
        "metadata": {
            "generation_timestamp": pd.Timestamp.now().isoformat(),
            "dashboard_version": "2.0",
            "components_analyzed": len(aggregated_data)
        }
    }
    
    logger.info(f"✅ Enhanced final dashboard generation complete")
    logger.info(f"📁 Dashboard directory: {output_dir}")
    logger.info(f"📄 Comprehensive report: {report_path}")
    logger.info(f"🖼 Visualizations: {len(visualizations)} generated")
    
    return results
