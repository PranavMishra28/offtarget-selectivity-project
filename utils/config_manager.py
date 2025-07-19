"""
Configuration Management System
Handles all configuration loading, validation, and access for the pipeline.
"""

import os
import yaml
import json
from typing import Dict, Any, Optional
from pathlib import Path
import structlog

class ConfigManager:
    """Centralized configuration management"""
    
    def __init__(self, config_path: str = "config.yaml"):
        self.config_path = config_path
        self.logger = structlog.get_logger(__name__)
        self._config = None
        self._load_config()
    
    def _load_config(self):
        """Load configuration from file"""
        try:
            if os.path.exists(self.config_path):
                with open(self.config_path, 'r') as f:
                    self._config = yaml.safe_load(f)
                self.logger.info(f"Loaded configuration from {self.config_path}")
            else:
                self._config = self._get_default_config()
                self._save_default_config()
                self.logger.info(f"Created default configuration at {self.config_path}")
        except Exception as e:
            self.logger.error(f"Failed to load configuration: {e}")
            self._config = self._get_default_config()
    
    def _get_default_config(self) -> Dict[str, Any]:
        """Get default configuration"""
        return {
            "apis": {
                "swiss_target_prediction": {
                    "base_url": "http://www.swisstargetprediction.ch",
                    "rate_limit": 10,
                    "timeout": 30
                },
                "askcos": {
                    "base_url": "https://askcos.mit.edu/api",
                    "rate_limit": 5,
                    "timeout": 60
                },
                "ibm_rxn": {
                    "base_url": "https://rxn.res.ibm.com/rxn/api",
                    "rate_limit": 10,
                    "timeout": 45
                },
                "alphafold": {
                    "base_url": "https://alphafold.ebi.ac.uk/api",
                    "rate_limit": 2,
                    "timeout": 300
                }
            },
            "models": {
                "moses": {
                    "enabled": True,
                    "model_path": "models/moses",
                    "max_molecules": 100
                },
                "chemformer": {
                    "enabled": True,
                    "model_path": "models/chemformer",
                    "max_molecules": 50
                },
                "rdkit": {
                    "enabled": True,
                    "max_molecules": 200
                }
            },
            "pipeline": {
                "scoring_weights": {
                    "empirical": 0.3,
                    "structural": 0.25,
                    "synthesis": 0.2,
                    "expression": 0.25
                },
                "thresholds": {
                    "min_confidence": 0.3,
                    "max_risk_score": 0.7,
                    "min_synthesis_score": 0.4
                }
            },
            "output": {
                "base_directory": "outputs",
                "save_intermediates": True,
                "generate_visualizations": True
            }
        }
    
    def _save_default_config(self):
        """Save default configuration to file"""
        try:
            with open(self.config_path, 'w') as f:
                yaml.dump(self._config, f, default_flow_style=False, indent=2)
        except Exception as e:
            self.logger.error(f"Failed to save default configuration: {e}")
    
    def get_config(self) -> Dict[str, Any]:
        """Get full configuration"""
        return self._config
    
    def get_api_config(self, api_name: Optional[str] = None) -> Dict[str, str]:
        """Get API endpoints configuration"""
        api_config = self._config.get("apis", {})
        if api_name:
            return api_config.get(api_name, {})
        return api_config
    
    def get_model_config(self, model_type: Optional[str] = None) -> Dict[str, Any]:
        """Get model configuration"""
        if model_type:
            return self._config.get("models", {}).get(model_type, {})
        return self._config.get("models", {})
    
    def get_pipeline_config(self) -> Dict[str, Any]:
        """Get pipeline configuration"""
        return self._config.get("pipeline", {})
    
    def get_output_config(self) -> Dict[str, Any]:
        """Get output configuration"""
        return self._config.get("output", {})
    
    def get_scoring_weights(self) -> Dict[str, float]:
        """Get scoring weights"""
        return self._config.get("pipeline", {}).get("scoring_weights", {})
    
    def get_thresholds(self) -> Dict[str, float]:
        """Get threshold values"""
        return self._config.get("pipeline", {}).get("thresholds", {})

# Global config manager instance
config_manager = ConfigManager() 