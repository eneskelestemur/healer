'''
    FastAPI routes for HEALER web application.
    
    Supports two modes:
    - Local mode (default): Jobs run synchronously, no Redis/Celery needed
    - Server mode: Jobs run via Celery workers with Redis backend
    
    Set HEALER_SERVER_MODE=true to enable server mode.
'''
import os
import uuid
import json
import io
from pathlib import Path
from typing import List, Optional, Dict, Any

import pandas as pd
from fastapi import APIRouter, HTTPException
from fastapi.responses import StreamingResponse
from pydantic import BaseModel
from rdkit import Chem
from rdkit.Chem import rdDepictor, Descriptors, QED

from healer.web.models import MoleculeRequest, SiteRequest, JobSubmitResponse, JobStatusResponse
from healer.web.interface import (
    SERVER_MODE,
    run_molecule_enumeration,
    run_site_enumeration,
    format_enumeration_results,
    get_server_limits,
    discover_building_blocks,
)
from healer.utils import utils

router = APIRouter(prefix="/api")

# ============================================================================
# Mode Detection
# ============================================================================

USE_CELERY = SERVER_MODE

# Try to import Celery components only if needed
celery_app = None
task_enumerate_molecule = None
task_enumerate_site = None
AsyncResult = None

if USE_CELERY:
    try:
        from celery.result import AsyncResult as _AsyncResult
        from healer.web.celery_worker import (
            celery_app as _celery_app,
            task_enumerate_molecule as _task_enumerate_molecule,
            task_enumerate_site as _task_enumerate_site,
        )
        celery_app = _celery_app
        task_enumerate_molecule = _task_enumerate_molecule
        task_enumerate_site = _task_enumerate_site
        AsyncResult = _AsyncResult
        print("HEALER Web: Running in SERVER mode (Celery/Redis)")
    except ImportError as e:
        print(f"Warning: Celery import failed ({e}), falling back to local mode")
        USE_CELERY = False

if not USE_CELERY:
    print("HEALER Web: Running in LOCAL mode (synchronous)")

# ============================================================================
# In-Memory Job Store (for local mode)
# ============================================================================

_local_jobs: Dict[str, Dict[str, Any]] = {}


def _run_job_sync(job_id: str, job_type: str, params: dict) -> None:
    """Run a job synchronously and store results."""
    _local_jobs[job_id] = {"status": "STARTED", "result": None, "error": None}
    
    try:
        if job_type == "molecule":
            raw_results = run_molecule_enumeration(**params)
            display_res, complete_res = format_enumeration_results(raw_results, 'molecule')
        else:  # site
            raw_results = run_site_enumeration(**params)
            display_res, complete_res = format_enumeration_results(raw_results, 'site')
        
        _local_jobs[job_id] = {
            "status": "SUCCESS",
            "result": {"display": display_res, "complete": complete_res},
            "error": None,
        }
    except Exception as e:
        _local_jobs[job_id] = {
            "status": "FAILURE",
            "result": None,
            "error": str(e),
        }


# ============================================================================
# Enumeration Endpoints
# ============================================================================

@router.post("/enumerate/molecule", response_model=JobSubmitResponse)
async def submit_molecule_enumeration(request: MoleculeRequest):
    params = request.model_dump() if hasattr(request, "model_dump") else request.dict()
    
    if USE_CELERY:
        task = task_enumerate_molecule.delay(params)
        return JobSubmitResponse(job_id=task.id, status="submitted")
    else:
        # Local synchronous mode
        job_id = str(uuid.uuid4())
        _run_job_sync(job_id, "molecule", params)
        return JobSubmitResponse(job_id=job_id, status="submitted")


@router.post("/enumerate/site", response_model=JobSubmitResponse)
async def submit_site_enumeration(request: SiteRequest):
    params = request.model_dump() if hasattr(request, "model_dump") else request.dict()
    
    if USE_CELERY:
        task = task_enumerate_site.delay(params)
        return JobSubmitResponse(job_id=task.id, status="submitted")
    else:
        # Local synchronous mode
        job_id = str(uuid.uuid4())
        _run_job_sync(job_id, "site", params)
        return JobSubmitResponse(job_id=job_id, status="submitted")


@router.get("/jobs/{job_id}", response_model=JobStatusResponse)
async def get_job_status(job_id: str):
    if USE_CELERY:
        task_result = AsyncResult(job_id, app=celery_app)
        response = JobStatusResponse(job_id=job_id, status=task_result.status)
        
        if task_result.status == 'SUCCESS':
            response.result = task_result.result
        elif task_result.status == 'FAILURE':
            response.error = str(task_result.result)
        
        return response
    else:
        # Local mode
        if job_id not in _local_jobs:
            raise HTTPException(status_code=404, detail="Job not found")
        
        job = _local_jobs[job_id]
        response = JobStatusResponse(job_id=job_id, status=job["status"])
        
        if job["status"] == "SUCCESS":
            response.result = job["result"]
        elif job["status"] == "FAILURE":
            response.error = job["error"]
        
        return response


@router.post("/jobs/{job_id}/cancel")
async def cancel_job(job_id: str):
    """Cancel a running or pending job."""
    if USE_CELERY:
        try:
            # Revoke the task - terminate=True kills running tasks
            celery_app.control.revoke(job_id, terminate=True)
            return {"job_id": job_id, "status": "cancelled"}
        except Exception as e:
            raise HTTPException(status_code=500, detail=f"Failed to cancel job: {str(e)}")
    else:
        # Local mode - can't really cancel synchronous jobs
        # But we can mark it as cancelled if it exists
        if job_id in _local_jobs:
            _local_jobs[job_id]["status"] = "CANCELLED"
            return {"job_id": job_id, "status": "cancelled", "note": "Local mode - job may have already completed"}
        raise HTTPException(status_code=400, detail="Cannot cancel jobs in local mode")


@router.get("/info/mode")
async def get_server_mode():
    """Return the current server mode (for UI to know if cancel is available)."""
    return {"mode": "celery" if USE_CELERY else "local"}


@router.get("/info/limits")
async def get_server_limits_endpoint():
    """Return the server parameter limits for UI validation."""
    return {
        "server_mode": SERVER_MODE,
        "limits": get_server_limits()
    }


@router.get("/info/building-blocks")
async def get_available_building_blocks():
    """Return list of available building block libraries."""
    return {"building_blocks": discover_building_blocks()}


@router.get("/jobs/{job_id}/download")
async def download_job_results(job_id: str):
    if USE_CELERY:
        task_result = AsyncResult(job_id, app=celery_app)
        if task_result.status != 'SUCCESS':
            raise HTTPException(status_code=400, detail="Job not completed or failed")
        results = task_result.result.get('complete', [])
    else:
        if job_id not in _local_jobs:
            raise HTTPException(status_code=404, detail="Job not found")
        job = _local_jobs[job_id]
        if job["status"] != "SUCCESS":
            raise HTTPException(status_code=400, detail="Job not completed or failed")
        results = job["result"].get('complete', [])
    
    if not results:
        raise HTTPException(status_code=404, detail="No results found")
    
    df = pd.DataFrame(results)
    stream = io.StringIO()
    df.to_csv(stream, index=False)
    
    response = StreamingResponse(iter([stream.getvalue()]), media_type="text/csv")
    response.headers["Content-Disposition"] = f"attachment; filename=healer_results_{job_id}.csv"
    return response


# ============================================================================
# Utility Endpoints
# ============================================================================

class SmilesRequest(BaseModel):
    smiles: str

class RenderRequest(BaseModel):
    smiles: str
    bbs: Optional[List[str]] = None
    alpha: float = 0.4
    bgColor: str = 'rgba(255, 255, 255, 1.0)'


@router.get("/utils/reaction-tags")
async def get_reaction_tags():
    """Return a list of available reaction tags."""
    try:
        # Reaction tags always come from package data
        healer_pkg = Path(__file__).parent.parent
        reaction_tags_path = healer_pkg / 'data' / 'reactions' / 'reaction_tags.txt'
        
        if not reaction_tags_path.exists():
            return ["amide coupling", "amide", "C-N bond formation", "C-N",
                    "alkylation", "N-arylation", "azole", "amination"]

        with open(reaction_tags_path, 'r') as f:
            tags = [line.strip() for line in f if line.strip()]
        
        # Exclude "all" in server mode
        if SERVER_MODE:
            tags = [tag for tag in tags if tag.lower() != 'all']
            
        return tags
            
    except Exception as e:
        print(f"Error loading reaction tags: {e}")
        return []


@router.post("/utils/smiles-to-mol")
async def smiles_to_molfile(request: SmilesRequest):
    try:
        smiles = request.smiles.strip()
        if not smiles:
            raise HTTPException(status_code=400, detail="No SMILES provided")
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise HTTPException(status_code=400, detail="Invalid SMILES format")

        rdDepictor.Compute2DCoords(mol)
        molblock = Chem.MolToMolBlock(mol)
        
        return {"molblock": molblock}
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error converting SMILES: {str(e)}")


@router.post("/utils/render-mol-with-indices")
async def render_mol_with_indices(request: SmilesRequest):
    """Return a base64 SVG of the molecule with atom indices labeled, plus properties."""
    try:
        smiles = request.smiles.strip()
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise HTTPException(status_code=400, detail="Invalid SMILES")

        props = {
            "MW": round(Descriptors.MolWt(mol), 2),
            "LogP": round(Descriptors.MolLogP(mol), 2),
            "HBA": Descriptors.NumHAcceptors(mol),
            "HBD": Descriptors.NumHDonors(mol),
            "TPSA": round(Descriptors.TPSA(mol), 2),
            "QED": round(QED.qed(mol), 3)
        }

        svg_data_uri = utils.get_svg_mol(mol, legend="", show_idx=True, width=250, height=125)
        return {"svg": svg_data_uri, "properties": props}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error rendering molecule: {str(e)}")


@router.post("/utils/render-result")
async def render_result(request: RenderRequest):
    """Return a base64 SVG of the result molecule, highlighting BBs if provided."""
    try:
        smiles = request.smiles.strip()
        alpha = float(request.alpha)
        # get tuple from bgColor: 'rgba(235, 64, 52, 0.06)' -> (235, 64, 52)
        bg_color = request.bgColor.replace(' ', '').replace('rgba(', '').replace(')', '').split(',')
        bg_color = tuple(int(c)/255.0 for c in bg_color[:3])  # ignore alpha for bg color
        if request.bbs:
            valid_bbs = [bb for bb in request.bbs if bb and bb.strip()]
            if valid_bbs:
                try:
                    svg_data_uri = utils.get_svg_mol_with_bbs(
                        smiles, valid_bbs, legend="", alpha=alpha, bg_color_for_transparency=bg_color
                    )
                    return {"svg": svg_data_uri}
                except Exception as e:
                    print(f"Error highlighting BBs: {e}")
        
        svg_data_uri = utils.get_svg_mol(smiles, legend="")
        return {"svg": svg_data_uri}
    except Exception as e:
        try:
            svg_data_uri = utils.get_svg_mol(request.smiles, legend="")
            return {"svg": svg_data_uri}
        except:
            raise HTTPException(status_code=500, detail=f"Error rendering result: {str(e)}")

