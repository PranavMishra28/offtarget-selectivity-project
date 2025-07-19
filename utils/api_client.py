"""
API Client Management System
Handles all external API interactions with proper error handling and fallback mechanisms.
"""

import asyncio
import aiohttp
import json
from typing import Dict, Any, Optional, List
from dataclasses import dataclass
import structlog
from utils.config_manager import config_manager

@dataclass
class APIResponse:
    """Standardized API response"""
    success: bool
    data: Dict[str, Any]
    error: Optional[str] = None

class APIClient:
    """Base API client with common functionality"""
    
    def __init__(self, base_url: str, timeout: int = 30):
        self.base_url = base_url
        self.timeout = timeout
        self.logger = structlog.get_logger(__name__)
        self.session = None
    
    async def _get_session(self) -> aiohttp.ClientSession:
        """Get or create HTTP session"""
        if self.session is None or self.session.closed:
            timeout = aiohttp.ClientTimeout(total=self.timeout)
            self.session = aiohttp.ClientSession(timeout=timeout)
        return self.session
    
    async def _make_request(self, method: str, endpoint: str, **kwargs) -> APIResponse:
        """Make HTTP request with error handling"""
        try:
            session = await self._get_session()
            url = f"{self.base_url}{endpoint}"
            
            async with session.request(method, url, **kwargs) as response:
                if response.status == 200:
                    data = await response.json()
                    return APIResponse(success=True, data=data)
                else:
                    error_msg = f"HTTP {response.status}: {response.reason}"
                    return APIResponse(success=False, data={}, error=error_msg)
                    
        except Exception as e:
            return APIResponse(success=False, data={}, error=str(e))
    
    async def close(self):
        """Close HTTP session"""
        if self.session and not self.session.closed:
            await self.session.close()

class SwissTargetClient(APIClient):
    """SwissTargetPrediction API client"""
    
    async def predict_targets(self, smiles: str, max_results: int = 10) -> APIResponse:
        """Predict targets using SwissTargetPrediction"""
        endpoint = "/predict"
        payload = {
            "smiles": smiles,
            "max_results": max_results
        }
        return await self._make_request("POST", endpoint, json=payload)

class SEAClient(APIClient):
    """SEA API client"""
    
    async def predict_targets(self, smiles: str, max_results: int = 10) -> APIResponse:
        """Predict targets using SEA"""
        endpoint = "/predict"
        payload = {
            "smiles": smiles,
            "max_results": max_results
        }
        return await self._make_request("POST", endpoint, json=payload)

class ChemProtClient(APIClient):
    """ChemProt API client"""
    
    async def predict_targets(self, smiles: str, max_results: int = 10) -> APIResponse:
        """Predict targets using ChemProt"""
        endpoint = "/predict"
        payload = {
            "smiles": smiles,
            "max_results": max_results
        }
        return await self._make_request("POST", endpoint, json=payload)

class ASKCOSClient(APIClient):
    """ASKCOS API client"""
    
    async def get_synthesis_routes(self, smiles: str, max_routes: int = 5) -> APIResponse:
        """Get synthesis routes from ASKCOS"""
        endpoint = "/retrosynthesis"
        payload = {
            "smiles": smiles,
            "max_routes": max_routes
        }
        return await self._make_request("POST", endpoint, json=payload)

class IBMRXNClient(APIClient):
    """IBM RXN API client"""
    
    async def get_retrosynthesis(self, smiles: str, max_paths: int = 3) -> APIResponse:
        """Get retrosynthesis paths from IBM RXN"""
        endpoint = "/retrosynthesis"
        payload = {
            "smiles": smiles,
            "max_paths": max_paths
        }
        return await self._make_request("POST", endpoint, json=payload)

class AlphaFoldClient(APIClient):
    """AlphaFold API client"""
    
    async def predict_structure(self, uniprot_id: str, sequence: str = "") -> APIResponse:
        """Predict protein structure using AlphaFold"""
        endpoint = "/predict"
        payload = {
            "uniprot_id": uniprot_id,
            "sequence": sequence
        }
        return await self._make_request("POST", endpoint, json=payload)

class APIManager:
    """Centralized API management"""
    
    def __init__(self):
        self.logger = structlog.get_logger(__name__)
        self.config = config_manager.get_api_config()
        self.clients = {}
        self._initialize_clients()
    
    def _initialize_clients(self):
        """Initialize API clients"""
        try:
            # Initialize SwissTargetPrediction
            if "swiss_target" in self.config:
                self.clients["swiss_target"] = SwissTargetClient(
                    self.config["swiss_target"]
                )
            
            # Initialize SEA
            if "sea" in self.config:
                self.clients["sea"] = SEAClient(
                    self.config["sea"]
                )
            
            # Initialize ChemProt
            if "chemprot" in self.config:
                self.clients["chemprot"] = ChemProtClient(
                    self.config["chemprot"]
                )
            
            # Initialize ASKCOS
            if "askcos" in self.config:
                self.clients["askcos"] = ASKCOSClient(
                    self.config["askcos"]
                )
            
            # Initialize IBM RXN
            if "ibm_rxn" in self.config:
                self.clients["ibm_rxn"] = IBMRXNClient(
                    self.config["ibm_rxn"]
                )
            
            # Initialize AlphaFold
            if "alphafold" in self.config:
                self.clients["alphafold"] = AlphaFoldClient(
                    self.config["alphafold"]
                )
            
            self.logger.info(f"Initialized {len(self.clients)} API clients")
            
        except Exception as e:
            self.logger.error(f"Failed to initialize API clients: {e}")
    
    def get_client(self, client_name: str) -> Optional[APIClient]:
        """Get API client by name"""
        return self.clients.get(client_name)
    
    async def close_all(self):
        """Close all API client sessions"""
        for client in self.clients.values():
            await client.close()

# Global API manager instance
api_manager = APIManager() 