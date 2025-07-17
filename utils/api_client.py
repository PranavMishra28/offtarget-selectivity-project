"""
Advanced API Client System for Off-Target & Selectivity Pipeline
Provides robust, rate-limited access to external APIs with comprehensive error handling.
"""

import asyncio
import aiohttp
import time
import json
from typing import Dict, Any, Optional, List
from dataclasses import dataclass
import structlog
from .config_manager import config_manager, APIConfig

@dataclass
class APIResponse:
    """Standardized API response format"""
    success: bool
    data: Optional[Dict[str, Any]] = None
    error: Optional[str] = None
    status_code: Optional[int] = None
    metadata: Optional[Dict[str, Any]] = None

class RateLimiter:
    """Simple rate limiter using asyncio"""
    
    def __init__(self, rate_limit: int, period: int = 60):
        self.rate_limit = rate_limit
        self.period = period
        self.requests = []
    
    async def acquire(self):
        """Acquire permission to make a request"""
        now = time.time()
        
        # Remove old requests outside the window
        self.requests = [req_time for req_time in self.requests if now - req_time < self.period]
        
        if len(self.requests) >= self.rate_limit:
            # Wait until we can make another request
            sleep_time = self.period - (now - self.requests[0])
            if sleep_time > 0:
                await asyncio.sleep(sleep_time)
        
        self.requests.append(time.time())

class CircuitBreaker:
    """Circuit breaker pattern for API resilience"""
    
    def __init__(self, failure_threshold: int = 5, recovery_timeout: int = 60):
        self.failure_threshold = failure_threshold
        self.recovery_timeout = recovery_timeout
        self.failure_count = 0
        self.last_failure_time = 0
        self.state = "CLOSED"  # CLOSED, OPEN, HALF_OPEN
        
    def can_execute(self) -> bool:
        """Check if the circuit breaker allows execution"""
        if self.state == "CLOSED":
            return True
        elif self.state == "OPEN":
            if time.time() - self.last_failure_time > self.recovery_timeout:
                self.state = "HALF_OPEN"
                return True
            return False
        else:  # HALF_OPEN
            return True
    
    def on_success(self):
        """Handle successful execution"""
        self.failure_count = 0
        self.state = "CLOSED"
        
    def on_failure(self):
        """Handle failed execution"""
        self.failure_count += 1
        self.last_failure_time = time.time()
        
        if self.failure_count >= self.failure_threshold:
            self.state = "OPEN"

class APIClient:
    """Base API client with common functionality"""
    
    def __init__(self, api_name: str):
        self.api_name = api_name
        self.config = config_manager.get_api_config(api_name)
        self.logger = structlog.get_logger(__name__)
        self.circuit_breaker = CircuitBreaker()
        self.rate_limiter = RateLimiter(rate_limit=self.config.rate_limit, period=60)
        
    async def _make_request(
        self, 
        method: str, 
        endpoint: str, 
        params: Optional[Dict[str, Any]] = None,
        data: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None
    ) -> APIResponse:
        """Make a rate-limited, circuit-breaker protected API request"""
        
        if not self.circuit_breaker.can_execute():
            return APIResponse(
                success=False,
                error=f"Circuit breaker is OPEN for {self.api_name}"
            )
        
        url = f"{self.config.base_url}{endpoint}"
        
        # Add API key if available
        if self.config.api_key:
            if headers is None:
                headers = {}
            headers["Authorization"] = f"Bearer {self.config.api_key}"
        
        try:
            # Apply rate limiting
            await self.rate_limiter.acquire()
            
            async with aiohttp.ClientSession() as session:
                async with session.request(
                    method=method,
                    url=url,
                    params=params,
                    json=data,
                    headers=headers,
                    timeout=aiohttp.ClientTimeout(total=self.config.timeout)
                ) as response:
                    
                    if response.status == 200:
                        data = await response.json()
                        self.circuit_breaker.on_success()
                        return APIResponse(
                            success=True,
                            data=data,
                            status_code=response.status,
                            metadata={"api": self.api_name, "endpoint": endpoint}
                        )
                    else:
                        error_text = await response.text()
                        self.circuit_breaker.on_failure()
                        return APIResponse(
                            success=False,
                            error=f"HTTP {response.status}: {error_text}",
                            status_code=response.status
                        )
                        
        except asyncio.TimeoutError:
            self.circuit_breaker.on_failure()
            return APIResponse(
                success=False,
                error=f"Timeout after {self.config.timeout}s"
            )
        except Exception as e:
            self.circuit_breaker.on_failure()
            return APIResponse(
                success=False,
                error=f"Request failed: {str(e)}"
            )

class SwissTargetPredictionClient(APIClient):
    """Client for SwissTargetPrediction API"""
    
    def __init__(self):
        super().__init__("swiss_target_prediction")
    
    async def predict_targets(self, smiles: str) -> APIResponse:
        """Predict protein targets for a SMILES string"""
        endpoint = "/predict"
        data = {
            "smiles": smiles,
            "species": "Homo sapiens"
        }
        return await self._make_request("POST", endpoint, data=data)

class ChEMBLClient(APIClient):
    """Client for ChEMBL API"""
    
    def __init__(self):
        super().__init__("chembl")
    
    async def get_similar_molecules(self, smiles: str, similarity: float = 0.7) -> APIResponse:
        """Get similar molecules from ChEMBL"""
        endpoint = f"/data/similarity/{smiles}/{int(similarity * 100)}"
        params = {"format": "json", "limit": 50}
        return await self._make_request("GET", endpoint, params=params)
    
    async def get_molecule_activities(self, chembl_id: str) -> APIResponse:
        """Get activities for a ChEMBL molecule"""
        endpoint = "/data/activity.json"
        params = {"molecule_chembl_id": chembl_id, "limit": 100}
        return await self._make_request("GET", endpoint, params=params)
    
    async def get_target_info(self, chembl_id: str) -> APIResponse:
        """Get target information from ChEMBL"""
        endpoint = f"/data/target/{chembl_id}.json"
        return await self._make_request("GET", endpoint)

class ASKCOSClient(APIClient):
    """Client for ASKCOS retrosynthesis API"""
    
    def __init__(self):
        super().__init__("askcos")
    
    async def get_retrosynthesis_routes(self, smiles: str) -> APIResponse:
        """Get retrosynthesis routes for a molecule"""
        endpoint = "/api/v2/retro/"
        data = {
            "smiles": smiles,
            "max_num_trees": 10,
            "max_branching": 25,
            "expansion_time": 60,
            "max_ppg": 10
        }
        return await self._make_request("POST", endpoint, data=data)

class IBMRXNClient(APIClient):
    """Client for IBM RXN retrosynthesis API"""
    
    def __init__(self):
        super().__init__("ibm_rxn")
    
    async def get_retrosynthesis_routes(self, smiles: str) -> APIResponse:
        """Get retrosynthesis routes using IBM RXN"""
        endpoint = "/api/rxn/v1/retrosynthesis"
        data = {
            "smiles": smiles,
            "max_steps": 5
        }
        return await self._make_request("POST", endpoint, data=data)

class AlphaFoldClient(APIClient):
    """Client for AlphaFold API"""
    
    def __init__(self):
        super().__init__("alphafold")
    
    async def predict_structure(self, sequence: str) -> APIResponse:
        """Predict protein structure using AlphaFold"""
        endpoint = "/api/predict"
        data = {
            "sequence": sequence,
            "model_type": "alphafold2"
        }
        return await self._make_request("POST", endpoint, data=data)

class APIManager:
    """Centralized API management"""
    
    def __init__(self):
        self.clients = {
            "swiss_target_prediction": SwissTargetPredictionClient(),
            "chembl": ChEMBLClient(),
            "askcos": ASKCOSClient(),
            "ibm_rxn": IBMRXNClient(),
            "alphafold": AlphaFoldClient()
        }
        self.logger = structlog.get_logger(__name__)
    
    def get_client(self, api_name: str) -> Optional[APIClient]:
        """Get API client by name"""
        return self.clients.get(api_name)
    
    async def batch_request(
        self, 
        api_name: str, 
        method: str, 
        requests: List[Dict[str, Any]]
    ) -> List[APIResponse]:
        """Execute batch requests with rate limiting"""
        client = self.get_client(api_name)
        if not client:
            return [APIResponse(success=False, error=f"Unknown API: {api_name}")] * len(requests)
        
        results = []
        for request in requests:
            result = await client._make_request(
                method=method,
                endpoint=request.get("endpoint", ""),
                params=request.get("params"),
                data=request.get("data"),
                headers=request.get("headers")
            )
            results.append(result)
        
        return results

# Global API manager instance
api_manager = APIManager() 