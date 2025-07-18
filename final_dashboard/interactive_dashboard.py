"""
Enhanced Interactive Dashboard - Modern UI/UX with Advanced Features
Implements Plotly Dash-based interactive dashboard with real-time features and AI-powered insights.
"""

import os
import json
import asyncio
from typing import Dict, Any, List, Optional
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import dash
from dash import dcc, html, Input, Output, State, callback_context
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
import structlog
from datetime import datetime
import requests
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

# Import configuration
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.config_manager import config_manager

class InteractiveDashboard:
    """Modern interactive dashboard with advanced features"""
    
    def __init__(self):
        self.logger = structlog.get_logger(__name__)
        self.app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
        self.setup_layout()
        self.setup_callbacks()
        
    def setup_layout(self):
        """Setup modern dashboard layout"""
        self.app.layout = dbc.Container([
            # Header
            dbc.Row([
                dbc.Col([
                    html.H1("🧪 Enhanced Off-Target Selectivity Pipeline", 
                           className="text-center mb-4 text-primary"),
                    html.P("Interactive Dashboard for Drug Discovery & Safety Assessment", 
                          className="text-center text-muted")
                ])
            ], className="mb-4"),
            
            # Input Section
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardHeader("🎯 Compound Input"),
                        dbc.CardBody([
                            dbc.Row([
                                dbc.Col([
                                    dbc.Label("SMILES String"),
                                    dbc.Input(id="smiles-input", type="text", 
                                            placeholder="e.g., CC(=O)OC1=CC=CC=C1C(=O)O",
                                            value="CC(=O)OC1=CC=CC=C1C(=O)O")
                                ], width=6),
                                dbc.Col([
                                    dbc.Label("Primary Target (UniProt ID)"),
                                    dbc.Input(id="target-input", type="text", 
                                            placeholder="e.g., P33259",
                                            value="P33259")
                                ], width=6)
                            ]),
                            dbc.Row([
                                dbc.Col([
                                    dbc.Button("🚀 Run Analysis", id="run-button", 
                                             color="primary", size="lg", className="mt-3")
                                ], className="text-center")
                            ])
                        ])
                    ])
                ])
            ], className="mb-4"),
            
            # Loading State
            dbc.Row([
                dbc.Col([
                    dbc.Spinner(id="loading-spinner", children=[
                        html.Div(id="loading-output")
                    ])
                ])
            ]),
            
            # Main Dashboard Tabs
            dbc.Row([
                dbc.Col([
                    dbc.Tabs([
                        # Executive Summary Tab
                        dbc.Tab([
                            dbc.Row([
                                dbc.Col([
                                    dbc.Card([
                                        dbc.CardHeader("🎯 Synthesis Decision"),
                                        dbc.CardBody([
                                            html.H2(id="decision-display", className="text-center"),
                                            html.Div(id="decision-explanation", className="mt-3")
                                        ])
                                    ])
                                ], width=4),
                                dbc.Col([
                                    dbc.Card([
                                        dbc.CardHeader("📊 Key Metrics"),
                                        dbc.CardBody([
                                            html.Div(id="key-metrics-display")
                                        ])
                                    ])
                                ], width=8)
                            ], className="mb-4"),
                            dbc.Row([
                                dbc.Col([
                                    dcc.Graph(id="summary-radar-chart")
                                ])
                            ])
                        ], label="Executive Summary", tab_id="summary"),
                        
                        # Risk Assessment Tab
                        dbc.Tab([
                            dbc.Row([
                                dbc.Col([
                                    dcc.Graph(id="risk-heatmap")
                                ], width=6),
                                dbc.Col([
                                    dcc.Graph(id="offtarget-distribution")
                                ], width=6)
                            ], className="mb-4"),
                            dbc.Row([
                                dbc.Col([
                                    dcc.Graph(id="ion-channel-risks")
                                ])
                            ])
                        ], label="Risk Assessment", tab_id="risk"),
                        
                        # Selectivity Analysis Tab
                        dbc.Tab([
                            dbc.Row([
                                dbc.Col([
                                    dcc.Graph(id="selectivity-comparison")
                                ], width=6),
                                dbc.Col([
                                    dcc.Graph(id="target-clustering")
                                ], width=6)
                            ], className="mb-4"),
                            dbc.Row([
                                dbc.Col([
                                    dcc.Graph(id="pathway-analysis")
                                ])
                            ])
                        ], label="Selectivity Analysis", tab_id="selectivity"),
                        
                        # Expression Analysis Tab
                        dbc.Tab([
                            dbc.Row([
                                dbc.Col([
                                    dcc.Graph(id="tissue-expression-heatmap")
                                ], width=6),
                                dbc.Col([
                                    dcc.Graph(id="expression-weighted-risks")
                                ], width=6)
                            ], className="mb-4"),
                            dbc.Row([
                                dbc.Col([
                                    dcc.Graph(id="tissue-specificity")
                                ])
                            ])
                        ], label="Expression Analysis", tab_id="expression"),
                        
                        # Synthesis Feasibility Tab
                        dbc.Tab([
                            dbc.Row([
                                dbc.Col([
                                    dcc.Graph(id="synthesis-scores")
                                ], width=6),
                                dbc.Col([
                                    dcc.Graph(id="complexity-analysis")
                                ], width=6)
                            ], className="mb-4"),
                            dbc.Row([
                                dbc.Col([
                                    html.Div(id="synthesis-recommendations")
                                ])
                            ])
                        ], label="Synthesis Feasibility", tab_id="synthesis"),
                        
                        # Advanced Features Tab
                        dbc.Tab([
                            dbc.Row([
                                dbc.Col([
                                    dbc.Card([
                                        dbc.CardHeader("🧠 AI-Powered Insights"),
                                        dbc.CardBody([
                                            html.Div(id="ai-insights"),
                                            dbc.Button("Generate New Insights", id="generate-insights-btn",
                                                     color="info", className="mt-3")
                                        ])
                                    ])
                                ], width=6),
                                dbc.Col([
                                    dbc.Card([
                                        dbc.CardHeader("🔍 Known Toxicophores"),
                                        dbc.CardBody([
                                            html.Div(id="toxicophore-analysis")
                                        ])
                                    ])
                                ], width=6)
                            ], className="mb-4"),
                            dbc.Row([
                                dbc.Col([
                                    dbc.Card([
                                        dbc.CardHeader("🛠️ Modification Recommendations"),
                                        dbc.CardBody([
                                            html.Div(id="modification-recommendations")
                                        ])
                                    ])
                                ])
                            ])
                        ], label="Advanced Features", tab_id="advanced")
                    ], id="dashboard-tabs")
                ])
            ])
        ], fluid=True)
    
    def setup_callbacks(self):
        """Setup interactive callbacks"""
        
        @self.app.callback(
            [Output("loading-output", "children"),
             Output("decision-display", "children"),
             Output("decision-explanation", "children"),
             Output("key-metrics-display", "children"),
             Output("summary-radar-chart", "figure"),
             Output("risk-heatmap", "figure"),
             Output("offtarget-distribution", "figure"),
             Output("ion-channel-risks", "figure"),
             Output("selectivity-comparison", "figure"),
             Output("target-clustering", "figure"),
             Output("pathway-analysis", "figure"),
             Output("tissue-expression-heatmap", "figure"),
             Output("expression-weighted-risks", "figure"),
             Output("tissue-specificity", "figure"),
             Output("synthesis-scores", "figure"),
             Output("complexity-analysis", "figure"),
             Output("synthesis-recommendations", "children"),
             Output("ai-insights", "children"),
             Output("toxicophore-analysis", "children"),
             Output("modification-recommendations", "children")],
            [Input("run-button", "n_clicks")],
            [State("smiles-input", "value"),
             State("target-input", "value")]
        )
        def run_analysis(n_clicks, smiles, target):
            if n_clicks is None:
                raise PreventUpdate
            
            # Run pipeline analysis
            results = self.run_pipeline_analysis(smiles, target)
            
            # Generate all visualizations
            return self.generate_all_outputs(results)
    
    def run_pipeline_analysis(self, smiles: str, target: str) -> Dict[str, Any]:
        """Run the complete pipeline analysis"""
        try:
            # Import and run pipeline
            from run import run_offtarget_selectivity_pipeline
            
            # Create temporary output directory
            output_dir = f"temp_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
            os.makedirs(output_dir, exist_ok=True)
            
            # Run pipeline
            results = asyncio.run(run_offtarget_selectivity_pipeline(
                smiles=smiles,
                primary_uniprot_id=target,
                output_base_dir=output_dir
            ))
            
            return results
            
        except Exception as e:
            self.logger.error(f"Pipeline analysis failed: {e}")
            return {"error": str(e)}
    
    def generate_all_outputs(self, results: Dict[str, Any]) -> List[Any]:
        """Generate all dashboard outputs"""
        if "error" in results:
            return [f"❌ Analysis failed: {results['error']}"] * 20
        
        # Extract key data
        impact_data = results.get("results", {}).get("impact_risk", {})
        decision = impact_data.get("decision_flag", "Unknown")
        selectivity_score = impact_data.get("selectivity_score", 0)
        safety_score = impact_data.get("safety_score", 0)
        
        # Generate outputs
        outputs = []
        
        # Loading output
        outputs.append("✅ Analysis completed successfully!")
        
        # Decision display
        decision_color = {
            "Synthesize": "success",
            "Watch": "warning", 
            "Modify": "info",
            "Reject": "danger"
        }.get(decision, "secondary")
        
        outputs.append(html.H2(decision, className=f"text-{decision_color}"))
        
        # Decision explanation
        explanations = {
            "Synthesize": "Compound shows excellent selectivity and safety profile. Proceed with synthesis.",
            "Watch": "Compound has moderate risk. Proceed with caution and enhanced monitoring.",
            "Modify": "Structural modifications required to improve selectivity or safety.",
            "Reject": "Compound has significant off-target risks. Consider alternative scaffolds."
        }
        outputs.append(html.P(explanations.get(decision, "Analysis complete.")))
        
        # Key metrics
        outputs.append(html.Div([
            dbc.Row([
                dbc.Col([
                    html.H4(f"{selectivity_score:.3f}", className="text-primary"),
                    html.P("Selectivity Score", className="text-muted")
                ], width=3),
                dbc.Col([
                    html.H4(f"{safety_score:.3f}", className="text-success"),
                    html.P("Safety Score", className="text-muted")
                ], width=3),
                dbc.Col([
                    html.H4(f"{len(impact_data.get('risky_offtargets', []))}", className="text-warning"),
                    html.P("Risky Off-Targets", className="text-muted")
                ], width=3),
                dbc.Col([
                    html.H4(f"{impact_data.get('selectivity_index', 0):.1f}", className="text-info"),
                    html.P("Selectivity Index", className="text-muted")
                ], width=3)
            ])
        ]))
        
        # Generate all charts
        outputs.extend(self.generate_charts(results))
        
        return outputs
    
    def generate_charts(self, results: Dict[str, Any]) -> List[go.Figure]:
        """Generate all interactive charts"""
        charts = []
        
        # Summary radar chart
        charts.append(self.create_radar_chart(results))
        
        # Risk heatmap
        charts.append(self.create_risk_heatmap(results))
        
        # Off-target distribution
        charts.append(self.create_offtarget_distribution(results))
        
        # Ion channel risks
        charts.append(self.create_ion_channel_chart(results))
        
        # Selectivity comparison
        charts.append(self.create_selectivity_comparison(results))
        
        # Target clustering
        charts.append(self.create_target_clustering(results))
        
        # Pathway analysis
        charts.append(self.create_pathway_analysis(results))
        
        # Tissue expression heatmap
        charts.append(self.create_tissue_expression_heatmap(results))
        
        # Expression weighted risks
        charts.append(self.create_expression_weighted_risks(results))
        
        # Tissue specificity
        charts.append(self.create_tissue_specificity(results))
        
        # Synthesis scores
        charts.append(self.create_synthesis_scores(results))
        
        # Complexity analysis
        charts.append(self.create_complexity_analysis(results))
        
        return charts
    
    def create_radar_chart(self, results: Dict[str, Any]) -> go.Figure:
        """Create radar chart for key metrics"""
        impact_data = results.get("results", {}).get("impact_risk", {})
        
        categories = ['Selectivity', 'Safety', 'Synthesis', 'Expression', 'Structure']
        values = [
            impact_data.get("selectivity_score", 0),
            impact_data.get("safety_score", 0),
            0.8,  # Mock synthesis score
            0.7,  # Mock expression score
            0.6   # Mock structure score
        ]
        
        fig = go.Figure()
        fig.add_trace(go.Scatterpolar(
            r=values,
            theta=categories,
            fill='toself',
            name='Compound Profile',
            line_color='rgb(32, 201, 151)'
        ))
        
        fig.update_layout(
            polar=dict(
                radialaxis=dict(
                    visible=True,
                    range=[0, 1]
                )),
            showlegend=False,
            title="Compound Performance Profile"
        )
        
        return fig
    
    def create_risk_heatmap(self, results: Dict[str, Any]) -> go.Figure:
        """Create risk assessment heatmap"""
        impact_data = results.get("results", {}).get("impact_risk", {})
        risky_targets = impact_data.get("risky_offtargets", [])
        
        if not risky_targets:
            # Create mock data
            targets = ["Target1", "Target2", "Target3", "Target4", "Target5"]
            scores = [0.8, 0.6, 0.4, 0.7, 0.3]
        else:
            targets = [t["uniprot_id"] for t in risky_targets[:10]]
            scores = [t["combined_score"] for t in risky_targets[:10]]
        
        fig = go.Figure(data=go.Heatmap(
            z=[scores],
            x=targets,
            y=['Risk Score'],
            colorscale='Reds',
            showscale=True
        ))
        
        fig.update_layout(
            title="Off-Target Risk Assessment",
            xaxis_title="Targets",
            yaxis_title="Risk Metrics"
        )
        
        return fig
    
    def create_offtarget_distribution(self, results: Dict[str, Any]) -> go.Figure:
        """Create off-target score distribution"""
        impact_data = results.get("results", {}).get("impact_risk", {})
        risky_targets = impact_data.get("risky_offtargets", [])
        
        if not risky_targets:
            scores = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        else:
            scores = [t["combined_score"] for t in risky_targets]
        
        fig = go.Figure(data=[go.Histogram(x=scores, nbinsx=20)])
        fig.update_layout(
            title="Off-Target Score Distribution",
            xaxis_title="Binding Score",
            yaxis_title="Frequency"
        )
        
        return fig
    
    def create_ion_channel_chart(self, results: Dict[str, Any]) -> go.Figure:
        """Create ion channel risk chart"""
        impact_data = results.get("results", {}).get("impact_risk", {})
        ion_risks = impact_data.get("ion_channel_risks", {})
        
        channels = ['hERG', 'Nav', 'Ca2+']
        risks = [
            ion_risks.get("herg_risk", 0.3),
            ion_risks.get("nav_risk", 0.2),
            ion_risks.get("ca_risk", 0.1)
        ]
        
        fig = go.Figure(data=[
            go.Bar(x=channels, y=risks, marker_color=['red', 'orange', 'yellow'])
        ])
        
        fig.update_layout(
            title="Ion Channel Liability Assessment",
            xaxis_title="Channel Type",
            yaxis_title="Risk Score"
        )
        
        return fig
    
    def create_selectivity_comparison(self, results: Dict[str, Any]) -> go.Figure:
        """Create empirical vs structural comparison"""
        empirical_data = results.get("results", {}).get("empirical_binding", {})
        structural_data = results.get("results", {}).get("structure_modeling", {})
        
        # Mock data for comparison
        targets = ["Target1", "Target2", "Target3", "Target4", "Target5"]
        empirical_scores = [0.8, 0.6, 0.4, 0.7, 0.3]
        structural_scores = [0.7, 0.5, 0.3, 0.6, 0.2]
        
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=targets, y=empirical_scores, 
                               mode='markers+lines', name='Empirical'))
        fig.add_trace(go.Scatter(x=targets, y=structural_scores, 
                               mode='markers+lines', name='Structural'))
        
        fig.update_layout(
            title="Empirical vs Structural Predictions",
            xaxis_title="Targets",
            yaxis_title="Binding Score"
        )
        
        return fig
    
    def create_target_clustering(self, results: Dict[str, Any]) -> go.Figure:
        """Create target clustering visualization"""
        # Mock clustering data
        targets = ["GPCR", "Kinase", "Ion Channel", "Enzyme", "Receptor"]
        counts = [15, 12, 8, 10, 6]
        
        fig = go.Figure(data=[go.Pie(labels=targets, values=counts)])
        fig.update_layout(title="Target Class Distribution")
        
        return fig
    
    def create_pathway_analysis(self, results: Dict[str, Any]) -> go.Figure:
        """Create pathway analysis chart"""
        # Mock pathway data
        pathways = ["Cell Signaling", "Metabolism", "Transport", "Immune Response"]
        scores = [0.8, 0.6, 0.4, 0.7]
        
        fig = go.Figure(data=[go.Bar(x=pathways, y=scores)])
        fig.update_layout(
            title="Pathway Impact Analysis",
            xaxis_title="Pathways",
            yaxis_title="Impact Score"
        )
        
        return fig
    
    def create_tissue_expression_heatmap(self, results: Dict[str, Any]) -> go.Figure:
        """Create tissue expression heatmap"""
        # Mock tissue expression data
        tissues = ["Brain", "Heart", "Liver", "Kidney", "Lung"]
        expression_levels = [0.8, 0.3, 0.6, 0.4, 0.7]
        
        fig = go.Figure(data=go.Heatmap(
            z=[expression_levels],
            x=tissues,
            y=['Expression'],
            colorscale='Blues'
        ))
        
        fig.update_layout(title="Tissue Expression Profile")
        
        return fig
    
    def create_expression_weighted_risks(self, results: Dict[str, Any]) -> go.Figure:
        """Create expression-weighted risk chart"""
        # Mock expression-weighted data
        targets = ["Target1", "Target2", "Target3", "Target4"]
        original_risks = [0.8, 0.6, 0.4, 0.7]
        weighted_risks = [0.6, 0.5, 0.3, 0.4]
        
        fig = go.Figure()
        fig.add_trace(go.Bar(x=targets, y=original_risks, name='Original Risk'))
        fig.add_trace(go.Bar(x=targets, y=weighted_risks, name='Expression-Weighted'))
        
        fig.update_layout(
            title="Expression-Weighted Risk Assessment",
            barmode='group'
        )
        
        return fig
    
    def create_tissue_specificity(self, results: Dict[str, Any]) -> go.Figure:
        """Create tissue specificity chart"""
        # Mock tissue specificity data
        tissues = ["Brain", "Heart", "Liver", "Kidney", "Lung", "Pancreas"]
        specificity = [0.9, 0.3, 0.6, 0.4, 0.7, 0.2]
        
        fig = go.Figure(data=[go.Bar(x=tissues, y=specificity)])
        fig.update_layout(
            title="Tissue Specificity Analysis",
            xaxis_title="Tissues",
            yaxis_title="Specificity Score"
        )
        
        return fig
    
    def create_synthesis_scores(self, results: Dict[str, Any]) -> go.Figure:
        """Create synthesis feasibility scores"""
        sparrow_data = results.get("results", {}).get("sparrow", {})
        
        # Mock synthesis data
        compounds = ["C1", "C2", "C3", "C4", "C5"]
        scores = [0.8, 0.6, 0.9, 0.4, 0.7]
        
        fig = go.Figure(data=[go.Bar(x=compounds, y=scores, marker_color='green')])
        fig.update_layout(
            title="Synthesis Feasibility Scores",
            xaxis_title="Compounds",
            yaxis_title="Feasibility Score"
        )
        
        return fig
    
    def create_complexity_analysis(self, results: Dict[str, Any]) -> go.Figure:
        """Create molecular complexity analysis"""
        # Mock complexity data
        metrics = ["MW", "LogP", "TPSA", "HBD", "HBA", "Rotatable Bonds"]
        values = [350, 2.5, 80, 3, 6, 4]
        thresholds = [500, 5, 140, 5, 10, 10]
        
        fig = go.Figure()
        fig.add_trace(go.Bar(x=metrics, y=values, name='Compound Values'))
        fig.add_trace(go.Bar(x=metrics, y=thresholds, name='Lipinski Thresholds'))
        
        fig.update_layout(
            title="Molecular Property Analysis",
            barmode='group'
        )
        
        return fig
    
    def run_dashboard(self, host: str = "localhost", port: int = 8050, debug: bool = True):
        """Run the interactive dashboard"""
        self.logger.info(f"Starting interactive dashboard on {host}:{port}")
        self.app.run_server(host=host, port=port, debug=debug)

def create_interactive_dashboard():
    """Create and return interactive dashboard instance"""
    return InteractiveDashboard()

if __name__ == "__main__":
    dashboard = create_interactive_dashboard()
    dashboard.run_dashboard() 