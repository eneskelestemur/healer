from fastapi import APIRouter, HTTPException, Body, Response
from fastapi.responses import StreamingResponse
from celery.result import AsyncResult
from rdkit import Chem
from rdkit.Chem import rdDepictor, Descriptors, QED
from pydantic import BaseModel
from typing import List, Optional, Dict
import pandas as pd
import io
import os
import json
from pathlib import Path

from healer.web.models import MoleculeRequest, SiteRequest, JobSubmitResponse, JobStatusResponse
from healer.utils import utils
from healer.web.celery_worker import celery_app, task_enumerate_molecule, task_enumerate_site

router = APIRouter(prefix="/api")

# ... (Submit endpoints remain unchanged) ...
@router.post("/enumerate/molecule", response_model=JobSubmitResponse)
async def submit_molecule_enumeration(request: MoleculeRequest):
    params = request.model_dump() if hasattr(request, "model_dump") else request.dict()
    task = task_enumerate_molecule.delay(params)
    return JobSubmitResponse(job_id=task.id, status="submitted")

@router.post("/enumerate/site", response_model=JobSubmitResponse)
async def submit_site_enumeration(request: SiteRequest):
    params = request.model_dump() if hasattr(request, "model_dump") else request.dict()
    task = task_enumerate_site.delay(params)
    return JobSubmitResponse(job_id=task.id, status="submitted")

@router.get("/jobs/{job_id}", response_model=JobStatusResponse)
async def get_job_status(job_id: str):
    task_result = AsyncResult(job_id, app=celery_app)
    
    response = JobStatusResponse(
        job_id=job_id,
        status=task_result.status,
    )
    
    if task_result.status == 'SUCCESS':
        response.result = task_result.result
    elif task_result.status == 'FAILURE':
        response.error = str(task_result.result)
    
    return response

@router.get("/jobs/{job_id}/download")
async def download_job_results(job_id: str):
    task_result = AsyncResult(job_id, app=celery_app)
    
    if task_result.status != 'SUCCESS':
        raise HTTPException(status_code=400, detail="Job not completed or failed")
    
    # Get complete results
    results = task_result.result.get('complete', [])
    if not results:
        raise HTTPException(status_code=404, detail="No results found")
        
    df = pd.DataFrame(results)
    
    # Convert to CSV
    stream = io.StringIO()
    df.to_csv(stream, index=False)
    
    response = StreamingResponse(iter([stream.getvalue()]), media_type="text/csv")
    response.headers["Content-Disposition"] = f"attachment; filename=healer_results_{job_id}.csv"
    return response

class SmilesRequest(BaseModel):
    smiles: str

class RenderRequest(BaseModel):
    smiles: str
    bbs: Optional[List[str]] = None

@router.get("/utils/reaction-tags")
async def get_reaction_tags():
    """Return a list of available reaction tags."""
    try:
        # Locate data directory similar to interface.py
        _env_data_dir = os.environ.get('HEALER_DATA_DIR')
        if _env_data_dir:
            reactions_path = Path(_env_data_dir) / 'reactions' / 'reactions.json'
        else:
            # Assuming we are in healer/web/routes.py, root is ../../
            healer_root = Path(__file__).parent.parent.parent
            reactions_path = healer_root / 'data' / 'reactions' / 'reactions.json'
        
        if not reactions_path.exists():
             # Fallback to hardcoded list if file not found
             return ["amide coupling", "amide", "C-N bond formation", "C-N",
                     "alkylation", "N-arylation", "azole", "amination"]

        with open(reactions_path, 'r') as f:
            data = json.load(f)
            # Return keys as tags, sorted
            return sorted(list(data.keys()))
            
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

        # Calculate properties
        props = {
            "MW": round(Descriptors.MolWt(mol), 2),
            "LogP": round(Descriptors.MolLogP(mol), 2),
            "HBA": Descriptors.NumHAcceptors(mol),
            "HBD": Descriptors.NumHDonors(mol),
            "TPSA": round(Descriptors.TPSA(mol), 2),
            "QED": round(QED.qed(mol), 3)
        }

        # Use utils.get_svg_mol with show_idx=True and slightly smaller size
        # This returns a data URI string "data:image/svg+xml;base64,..."
        svg_data_uri = utils.get_svg_mol(mol, legend="", show_idx=True, width=250, height=125)
        return {"svg": svg_data_uri, "properties": props}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error rendering molecule: {str(e)}")

@router.post("/utils/render-result")
async def render_result(request: RenderRequest):
    """Return a base64 SVG of the result molecule, highlighting BBs if provided."""
    try:
        smiles = request.smiles.strip()
        if request.bbs:
            # Filter out empty strings
            valid_bbs = [bb for bb in request.bbs if bb and bb.strip()]
            if valid_bbs:
                try:
                    # Removed legend per user request
                    svg_data_uri = utils.get_svg_mol_with_bbs(smiles, valid_bbs, legend="")
                    return {"svg": svg_data_uri}
                except Exception as e:
                    print(f"Error highlighting BBs: {e}")
                    # Fallback to normal rendering if highlighting fails
        
        # Fallback to normal rendering
        svg_data_uri = utils.get_svg_mol(smiles, legend="")
        return {"svg": svg_data_uri}
    except Exception as e:
        # Return a placeholder or simple render on error
        try:
             svg_data_uri = utils.get_svg_mol(request.smiles, legend="")
             return {"svg": svg_data_uri}
        except:
             raise HTTPException(status_code=500, detail=f"Error rendering result: {str(e)}")
