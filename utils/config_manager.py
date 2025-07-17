"""
Configuration Management System for Off-Target & Selectivity Pipeline
Handles loading, validation, and access to all pipeline settings.
"""

import os
import yaml
from pathlib import Path
from typing import Dict, Any, Optional
from dataclasses import dataclass
from dotenv import load_dotenv
import structlog

# Load environment variables
load_dotenv()

@dataclass
class APIConfig:
    """API configuration settings"""
    base_url: str
    rate_limit: int
    timeout: int
    api_key: Optional[str] = None

@dataclass
class ModelConfig:
    """Model configuration settings"""
    library: str
    num_variants: int
    novelty_threshold: float
    validity_check: bool
    pains_filter: bool

@dataclass
class PipelineConfig:
    """Pipeline configuration settings"""
    scoring_weights: Dict[str, float]
    decision_thresholds: Dict[str, float]
    tissue_weights: Dict[str, float]

class ConfigManager:
    """Centralized configuration management"""
    
    def __init__(self, config_path: str = "config.yaml"):
        self.config_path = Path(config_path)
        self.logger = structlog.get_logger()
        self.config = self._load_config()
        
    def _load_config(self) -> Dict[str, Any]:
        """Load configuration from YAML file"""
        if not self.config_path.exists():
            self.logger.warning(f"Config file {self.config_path} not found, using defaults")
            return self._get_default_config()
            
        try:
            with open(self.config_path, 'r') as f:
                config = yaml.safe_load(f)
            self.logger.info(f"Loaded configuration from {self.config_path}")
            return config
        except Exception as e:
            self.logger.error(f"Failed to load config: {e}")
            return self._get_default_config()
    
    def _get_default_config(self) -> Dict[str, Any]:
        """Default configuration settings"""
        return {
            "apis": {
                "swiss_target_prediction": {
                    "base_url": "http://www.swisstargetprediction.ch",
                    "rate_limit": 10,
                    "timeout": 30
                },
                "chembl": {
                    "base_url": "https://www.ebi.ac.uk/chembl/api",
                    "rate_limit": 60,
                    "timeout": 15
                }
            },
            "models": {
                "generative": {
                    "library": "moses",
                    "num_variants": 50,
                    "novelty_threshold": 0.7,
                    "validity_check": True,
                    "pains_filter": True
                }
            },
            "pipeline": {
                "scoring_weights": {
                    "empirical_confidence": 0.4,
                    "structural_binding_score": 0.3,
                    "tissue_expression_weight": 0.2,
                    "synthesis_feasibility": 0.1
                },
                "decision_thresholds": {
                    "synthesize": 2.0,
                    "watch": 1.0,
                    "reject": 0.5
                }
            },
            "performance": {
                "max_concurrent_requests": 10,
                "cache_ttl": 3600,
                "batch_size": 50
            }
        }
    
    def get_api_config(self, api_name: str) -> APIConfig:
        """Get API configuration for specified service"""
        api_config = self.config.get("apis", {}).get(api_name, {})
        
        # Check for environment variable overrides
        api_key = os.getenv(f"{api_name.upper()}_API_KEY")
        
        return APIConfig(
            base_url=api_config.get("base_url", ""),
            rate_limit=api_config.get("rate_limit", 10),
            timeout=api_config.get("timeout", 30),
            api_key=api_key
        )
    
    def get_model_config(self, model_type: str) -> ModelConfig:
        """Get model configuration for specified type"""
        model_config = self.config.get("models", {}).get(model_type, {})
        
        return ModelConfig(
            library=model_config.get("library", "moses"),
            num_variants=model_config.get("num_variants", 50),
            novelty_threshold=model_config.get("novelty_threshold", 0.7),
            validity_check=model_config.get("validity_check", True),
            pains_filter=model_config.get("pains_filter", True)
        )
    
    def get_pipeline_config(self) -> PipelineConfig:
        """Get pipeline configuration"""
        pipeline_config = self.config.get("pipeline", {})
        
        return PipelineConfig(
            scoring_weights=pipeline_config.get("scoring_weights", {}),
            decision_thresholds=pipeline_config.get("decision_thresholds", {}),
            tissue_weights=pipeline_config.get("tissue_weights", {})
        )
    
    def get_performance_config(self) -> Dict[str, Any]:
        """Get performance configuration"""
        return self.config.get("performance", {})
    
    def get_logging_config(self) -> Dict[str, Any]:
        """Get logging configuration"""
        return self.config.get("logging", {})

# Global configuration instance
config_manager = ConfigManager() 