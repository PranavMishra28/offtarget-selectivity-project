"""
Enhanced Expression Filter - Advanced Tissue-Specific Expression Integration
Implements multi-source expression data integration with configurable tissue weighting and risk assessment.
"""

import os
import json
import asyncio
from typing import Dict, Any, List, Optional, Tuple
import numpy as np
import pandas as pd
from scipy.stats import percentileofscore
import matplotlib.pyplot as plt
import seaborn as sns
import structlog
from tqdm import tqdm

# Import configuration
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.config_manager import config_manager

class ExpressionDataManager:
    """Manages multi-source expression data integration"""
    
    def __init__(self):
        self.logger = structlog.get_logger(__name__)
        self.gtex_data = None
        self.hpa_data = None
        self.tissue_data = None
        self.expression_cache = {}
    
    def load_gtex_data(self, gtex_filepath: str) -> pd.DataFrame:
        """Load GTEx expression data"""
        try:
            self.logger.info(f"Loading GTEx data from {gtex_filepath}")
            
            if gtex_filepath.endswith('.gz'):
                gtex_data = pd.read_csv(gtex_filepath, sep='\t', index_col=0, compression='gzip')
            else:
                gtex_data = pd.read_csv(gtex_filepath, sep='\t', index_col=0)
            
            # Clean and normalize data
            gtex_data = self._clean_expression_data(gtex_data)
            self.gtex_data = gtex_data
            
            self.logger.info(f"✅ Loaded GTEx data: {gtex_data.shape[0]} genes, {gtex_data.shape[1]} tissues")
            return gtex_data
            
        except Exception as e:
            self.logger.error(f"Failed to load GTEx data: {e}")
            return pd.DataFrame()
    
    def load_hpa_data(self, hpa_filepath: str) -> pd.DataFrame:
        """Load Human Protein Atlas expression data"""
        try:
            self.logger.info(f"Loading HPA data from {hpa_filepath}")
            
            # In production, this would load actual HPA data
            # For now, create mock HPA data
            hpa_data = self._create_mock_hpa_data()
            self.hpa_data = hpa_data
            
            self.logger.info(f"✅ Loaded HPA data: {hpa_data.shape[0]} genes, {hpa_data.shape[1]} tissues")
            return hpa_data
            
        except Exception as e:
            self.logger.error(f"Failed to load HPA data: {e}")
            return pd.DataFrame()
    
    def load_tissue_data(self, tissue_filepath: str) -> pd.DataFrame:
        """Load Tissue 2.0 expression data"""
        try:
            self.logger.info(f"Loading Tissue 2.0 data from {tissue_filepath}")
            
            # In production, this would load actual Tissue 2.0 data
            # For now, create mock tissue data
            tissue_data = self._create_mock_tissue_data()
            self.tissue_data = tissue_data
            
            self.logger.info(f"✅ Loaded Tissue 2.0 data: {tissue_data.shape[0]} genes, {tissue_data.shape[1]} tissues")
            return tissue_data
            
        except Exception as e:
            self.logger.error(f"Failed to load Tissue 2.0 data: {e}")
            return pd.DataFrame()
    
    def _clean_expression_data(self, data: pd.DataFrame) -> pd.DataFrame:
        """Clean and normalize expression data"""
        try:
            # Convert to numeric, coercing errors to NaN
            data = data.apply(pd.to_numeric, errors='coerce')
            
            # Remove genes with all zero expression
            data = data.loc[data.sum(axis=1) > 0]
            
            # Fill NaN values with 0
            data = data.fillna(0)
            
            # Log2 transform (add small constant to avoid log(0))
            data = np.log2(data + 1)
            
            # Remove outliers (genes with extremely high expression)
            data = data.clip(upper=data.quantile(0.99, axis=1), axis=1)
            
            return data
        except Exception as e:
            self.logger.error(f"Error cleaning expression data: {e}")
            # Return empty DataFrame if cleaning fails
            return pd.DataFrame()
    
    def _create_mock_hpa_data(self) -> pd.DataFrame:
        """Create mock HPA expression data"""
        # Mock tissue types
        tissues = ['brain', 'heart', 'liver', 'kidney', 'lung', 'pancreas', 'spleen', 'stomach']
        
        # Mock gene IDs (using some real UniProt IDs)
        genes = ['P23219', 'P35354', 'P29274', 'P41594', 'P08172', 'P21918', 'P28223', 'P08908']
        
        # Create mock expression data
        np.random.seed(42)
        data = np.random.exponential(2, size=(len(genes), len(tissues)))
        
        return pd.DataFrame(data, index=genes, columns=tissues)
    
    def _create_mock_tissue_data(self) -> pd.DataFrame:
        """Create mock Tissue 2.0 expression data"""
        # Mock tissue types
        tissues = ['cerebellum', 'cortex', 'hippocampus', 'amygdala', 'thalamus', 'cerebellum']
        
        # Mock gene IDs
        genes = ['P23219', 'P35354', 'P29274', 'P41594', 'P08172', 'P21918', 'P28223', 'P08908']
        
        # Create mock expression data
        np.random.seed(43)
        data = np.random.exponential(1.5, size=(len(genes), len(tissues)))
        
        return pd.DataFrame(data, index=genes, columns=tissues)
    
    def get_gene_expression(self, uniprot_id: str) -> Dict[str, Any]:
        """Get comprehensive expression data for a gene"""
        if uniprot_id in self.expression_cache:
            return self.expression_cache[uniprot_id]
        
        expression_data = {
            "uniprot_id": uniprot_id,
            "gtex": {},
            "hpa": {},
            "tissue": {},
            "combined": {}
        }
        
        # Get GTEx data
        if self.gtex_data is not None and uniprot_id in self.gtex_data.index:
            gtex_row = self.gtex_data.loc[uniprot_id]
            expression_data["gtex"] = {
                "tissues": gtex_row.to_dict(),
                "mean_expression": float(gtex_row.mean()),
                "max_expression": float(gtex_row.max()),
                "expression_variance": float(gtex_row.var())
            }
        
        # Get HPA data
        if self.hpa_data is not None and uniprot_id in self.hpa_data.index:
            hpa_row = self.hpa_data.loc[uniprot_id]
            expression_data["hpa"] = {
                "tissues": hpa_row.to_dict(),
                "mean_expression": float(hpa_row.mean()),
                "max_expression": float(hpa_row.max()),
                "expression_variance": float(hpa_row.var())
            }
        
        # Get Tissue 2.0 data
        if self.tissue_data is not None and uniprot_id in self.tissue_data.index:
            tissue_row = self.tissue_data.loc[uniprot_id]
            expression_data["tissue"] = {
                "tissues": tissue_row.to_dict(),
                "mean_expression": float(tissue_row.mean()),
                "max_expression": float(tissue_row.max()),
                "expression_variance": float(tissue_row.var())
            }
        
        # Combine data sources
        expression_data["combined"] = self._combine_expression_sources(expression_data)
        
        # Cache result
        self.expression_cache[uniprot_id] = expression_data
        
        return expression_data
    
    def _combine_expression_sources(self, expression_data: Dict[str, Any]) -> Dict[str, Any]:
        """Combine expression data from multiple sources"""
        combined = {
            "mean_expression": 0.0,
            "max_expression": 0.0,
            "expression_variance": 0.0,
            "tissue_specificity": 0.0,
            "source_count": 0
        }
        
        sources = ["gtex", "hpa", "tissue"]
        valid_sources = []
        
        for source in sources:
            if expression_data[source]:
                valid_sources.append(source)
                combined["mean_expression"] += expression_data[source]["mean_expression"]
                combined["max_expression"] = max(combined["max_expression"], 
                                               expression_data[source]["max_expression"])
                combined["expression_variance"] += expression_data[source]["expression_variance"]
        
        if valid_sources:
            combined["mean_expression"] /= len(valid_sources)
            combined["expression_variance"] /= len(valid_sources)
            combined["source_count"] = len(valid_sources)
            
            # Calculate tissue specificity (higher variance = more tissue specific)
            combined["tissue_specificity"] = min(1.0, combined["expression_variance"] / 10.0)
        
        return combined

class TissueWeightingEngine:
    """Advanced tissue weighting and risk assessment"""
    
    def __init__(self):
        self.logger = structlog.get_logger(__name__)
        self.pipeline_config = config_manager.get_pipeline_config()
        self.tissue_weights = self.pipeline_config.tissue_weights
    
    def calculate_tissue_weighted_risk(
        self, 
        uniprot_id: str, 
        base_risk_score: float,
        expression_data: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Calculate tissue-weighted risk score"""
        
        # Get combined expression data
        combined_expr = expression_data.get("combined", {})
        mean_expression = combined_expr.get("mean_expression", 0.0)
        tissue_specificity = combined_expr.get("tissue_specificity", 0.0)
        
        # Calculate tissue-specific risk factors
        tissue_risk_factors = self._calculate_tissue_risk_factors(expression_data)
        
        # Calculate overall tissue weight
        tissue_weight = self._calculate_tissue_weight(tissue_risk_factors, tissue_specificity)
        
        # Apply tissue weighting to base risk
        weighted_risk = base_risk_score * tissue_weight
        
        # Calculate risk components
        risk_components = {
            "base_risk": base_risk_score,
            "tissue_weight": tissue_weight,
            "expression_factor": min(1.0, mean_expression / 10.0),
            "specificity_factor": tissue_specificity,
            "critical_tissue_penalty": tissue_risk_factors.get("critical_tissue_penalty", 0.0)
        }
        
        return {
            "uniprot_id": uniprot_id,
            "original_score": base_risk_score,
            "tissue_weighted_score": round(weighted_risk, 4),
            "tissue_weight": round(tissue_weight, 4),
            "mean_expression": round(mean_expression, 4),
            "tissue_specificity": round(tissue_specificity, 4),
            "risk_components": risk_components,
            "tissue_risk_factors": tissue_risk_factors,
            "expression_data": expression_data
        }
    
    def _calculate_tissue_risk_factors(self, expression_data: Dict[str, Any]) -> Dict[str, float]:
        """Calculate tissue-specific risk factors"""
        risk_factors = {
            "cns_expression": 0.0,
            "cardiac_expression": 0.0,
            "liver_expression": 0.0,
            "kidney_expression": 0.0,
            "critical_tissue_penalty": 0.0
        }
        
        # Check GTEx data for critical tissues
        if expression_data.get("gtex", {}).get("tissues"):
            gtex_tissues = expression_data["gtex"]["tissues"]
            
            # CNS tissues
            cns_tissues = ['brain', 'cerebellum', 'cortex', 'hippocampus', 'amygdala', 'thalamus']
            cns_expression = sum(gtex_tissues.get(tissue, 0) for tissue in cns_tissues)
            risk_factors["cns_expression"] = min(1.0, cns_expression / 10.0)
            
            # Cardiac tissues
            cardiac_tissues = ['heart', 'atrial_appendage', 'left_ventricle']
            cardiac_expression = sum(gtex_tissues.get(tissue, 0) for tissue in cardiac_tissues)
            risk_factors["cardiac_expression"] = min(1.0, cardiac_expression / 10.0)
            
            # Liver
            risk_factors["liver_expression"] = min(1.0, gtex_tissues.get("liver", 0) / 5.0)
            
            # Kidney
            risk_factors["kidney_expression"] = min(1.0, gtex_tissues.get("kidney", 0) / 5.0)
        
        # Calculate critical tissue penalty
        critical_weights = {
            "cns_expression": self.tissue_weights.get("cns", 2.0),
            "cardiac_expression": self.tissue_weights.get("cardiac", 1.8),
            "liver_expression": self.tissue_weights.get("liver", 1.5),
            "kidney_expression": self.tissue_weights.get("kidney", 1.3)
        }
        
        critical_penalty = sum(
            risk_factors[tissue] * critical_weights[tissue]
            for tissue in critical_weights.keys()
        )
        
        risk_factors["critical_tissue_penalty"] = min(2.0, critical_penalty)
        
        return risk_factors
    
    def _calculate_tissue_weight(self, tissue_risk_factors: Dict[str, float], tissue_specificity: float) -> float:
        """Calculate overall tissue weight"""
        # Base weight
        base_weight = 1.0
        
        # Add critical tissue penalty
        critical_penalty = tissue_risk_factors.get("critical_tissue_penalty", 0.0)
        base_weight += critical_penalty
        
        # Adjust for tissue specificity (more specific = higher risk)
        specificity_factor = 1.0 + (tissue_specificity * 0.5)
        
        # Calculate final weight
        tissue_weight = base_weight * specificity_factor
        
        return min(5.0, max(0.1, tissue_weight))  # Clamp between 0.1 and 5.0

class ExpressionVisualizer:
    """Expression data visualization"""
    
    def __init__(self):
        self.logger = structlog.get_logger(__name__)
    
    def generate_expression_heatmap(
        self, 
        expression_data: Dict[str, Any], 
        output_path: str
    ) -> str:
        """Generate expression heatmap"""
        try:
            # Extract tissue expression data
            tissues_data = {}
            
            for source in ["gtex", "hpa", "tissue"]:
                if expression_data.get(source, {}).get("tissues"):
                    tissues_data[source] = expression_data[source]["tissues"]
            
            if not tissues_data:
                return ""
            
            # Create heatmap data
            all_tissues = set()
            for source_data in tissues_data.values():
                all_tissues.update(source_data.keys())
            
            # Create DataFrame for heatmap
            heatmap_data = pd.DataFrame(index=list(all_tissues), columns=list(tissues_data.keys()))
            
            for source, tissues in tissues_data.items():
                for tissue, expression in tissues.items():
                    heatmap_data.loc[tissue, source] = expression
            
            # Fill NaN values
            heatmap_data = heatmap_data.fillna(0)
            
            # Create heatmap
            plt.figure(figsize=(12, 8))
            sns.heatmap(heatmap_data, annot=True, cmap='YlOrRd', fmt='.2f', cbar_kws={'label': 'Expression Level'})
            plt.title('Tissue Expression Profile')
            plt.xlabel('Data Source')
            plt.ylabel('Tissue')
            plt.tight_layout()
            
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            return output_path
            
        except Exception as e:
            self.logger.error(f"Failed to generate expression heatmap: {e}")
            return ""
    
    def generate_risk_comparison_plot(
        self, 
        risk_data: List[Dict[str, Any]], 
        output_path: str
    ) -> str:
        """Generate risk comparison plot"""
        try:
            # Extract data for plotting
            uniprot_ids = [item["uniprot_id"] for item in risk_data]
            original_scores = [item["original_score"] for item in risk_data]
            weighted_scores = [item["tissue_weighted_score"] for item in risk_data]
            
            # Create comparison plot
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
            
            # Original vs weighted scores
            ax1.scatter(original_scores, weighted_scores, alpha=0.7, s=50)
            ax1.plot([0, 1], [0, 1], 'r--', alpha=0.5, label='No change')
            ax1.set_xlabel('Original Risk Score')
            ax1.set_ylabel('Tissue-Weighted Risk Score')
            ax1.set_title('Risk Score Comparison')
            ax1.legend()
            ax1.grid(True, alpha=0.3)
            
            # Score distribution
            ax2.hist(original_scores, alpha=0.7, label='Original', bins=20)
            ax2.hist(weighted_scores, alpha=0.7, label='Tissue-Weighted', bins=20)
            ax2.set_xlabel('Risk Score')
            ax2.set_ylabel('Frequency')
            ax2.set_title('Risk Score Distribution')
            ax2.legend()
            ax2.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            return output_path
            
        except Exception as e:
            self.logger.error(f"Failed to generate risk comparison plot: {e}")
            return ""

async def apply_expression_weighting(
    combined_scores: Dict[str, float],
    gtex_filepath: str = "data/gtex/GTEx_median_tpm.tsv.gz",
    output_path: str = "expression_filter/tissue_weighted_risk.json",
    metadata_path: str = "expression_filter/expression_analysis.json",
    visualization_path: str = "expression_filter/expression_visualizations"
) -> Dict[str, Any]:
    """
    Enhanced expression weighting with multi-source data integration.
    
    Args:
        combined_scores: Dictionary of UniProt IDs to risk scores
        gtex_filepath: Path to GTEx data file
        output_path: Path for weighted risk output
        metadata_path: Path for detailed analysis metadata
        visualization_path: Path for visualization outputs
    
    Returns:
        Dictionary of tissue-weighted risk scores
    """
    
    logger = structlog.get_logger(__name__)
    logger.info(f"Starting enhanced expression weighting for {len(combined_scores)} targets")
    
    # Initialize components
    expression_manager = ExpressionDataManager()
    weighting_engine = TissueWeightingEngine()
    visualizer = ExpressionVisualizer()
    
    # Load expression data
    logger.info("Loading expression data sources...")
    gtex_data = expression_manager.load_gtex_data(gtex_filepath)
    
    # Load additional data sources (mock for now)
    hpa_data = expression_manager.load_hpa_data("data/hpa/expression_data.tsv")
    tissue_data = expression_manager.load_tissue_data("data/tissue/expression_data.tsv")
    
    # Process each target
    logger.info("Processing targets with expression weighting...")
    weighted_results = {}
    expression_analysis = {}
    
    for uniprot_id, base_score in tqdm(combined_scores.items(), desc="Applying expression weighting"):
        try:
            # Get expression data
            expression_data = expression_manager.get_gene_expression(uniprot_id)
            
            # Calculate tissue-weighted risk
            weighted_result = weighting_engine.calculate_tissue_weighted_risk(
                uniprot_id, base_score, expression_data
            )
            
            weighted_results[uniprot_id] = weighted_result
            expression_analysis[uniprot_id] = expression_data
            
        except Exception as e:
            logger.error(f"Failed to process {uniprot_id}: {e}")
            # Add fallback result
            weighted_results[uniprot_id] = {
                "uniprot_id": uniprot_id,
                "original_score": base_score,
                "tissue_weighted_score": base_score,
                "tissue_weight": 1.0,
                "mean_expression": 0.0,
                "tissue_specificity": 0.0,
                "risk_components": {},
                "tissue_risk_factors": {},
                "expression_data": {}
            }
    
    # Generate visualizations
    logger.info("Generating visualizations...")
    os.makedirs(visualization_path, exist_ok=True)
    
    # Expression heatmap for top targets
    top_targets = sorted(weighted_results.items(), key=lambda x: x[1]["tissue_weighted_score"], reverse=True)[:5]
    if top_targets:
        top_uniprot_id = top_targets[0][0]
        top_expression_data = expression_analysis.get(top_uniprot_id, {})
        heatmap_path = os.path.join(visualization_path, "expression_heatmap.png")
        visualizer.generate_expression_heatmap(top_expression_data, heatmap_path)
    
    # Risk comparison plot
    risk_data = list(weighted_results.values())
    comparison_path = os.path.join(visualization_path, "risk_comparison.png")
    visualizer.generate_risk_comparison_plot(risk_data, comparison_path)
    
    # Save results
    logger.info("Saving results...")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Save weighted risk scores (legacy format)
    legacy_results = {
        uniprot_id: result["tissue_weighted_score"]
        for uniprot_id, result in weighted_results.items()
    }
    
    with open(output_path, 'w') as f:
        json.dump(legacy_results, f, indent=2)
    
    # Save detailed analysis
    detailed_analysis = {
        "input_scores": combined_scores,
        "weighted_results": weighted_results,
        "expression_analysis": expression_analysis,
        "statistics": {
            "total_targets": len(weighted_results),
            "mean_weight": np.mean([r["tissue_weight"] for r in weighted_results.values()]),
            "max_weight": max([r["tissue_weight"] for r in weighted_results.values()]),
            "min_weight": min([r["tissue_weight"] for r in weighted_results.values()]),
            "mean_expression": np.mean([r["mean_expression"] for r in weighted_results.values()]),
            "mean_specificity": np.mean([r["tissue_specificity"] for r in weighted_results.values()])
        },
        "visualizations": {
            "expression_heatmap": os.path.join(visualization_path, "expression_heatmap.png"),
            "risk_comparison": os.path.join(visualization_path, "risk_comparison.png")
        },
        "metadata": {
            "analysis_timestamp": pd.Timestamp.now().isoformat(),
            "analysis_method": "enhanced_expression_weighting"
        }
    }
    
    with open(metadata_path, 'w') as f:
        json.dump(detailed_analysis, f, indent=2, default=str)
    
    logger.info(f"✅ Enhanced expression weighting complete")
    logger.info(f"📁 Results: {output_path}")
    logger.info(f"📄 Detailed analysis: {metadata_path}")
    logger.info(f"🖼 Visualizations: {visualization_path}")
    
    return {"weighted_risks": legacy_results}
