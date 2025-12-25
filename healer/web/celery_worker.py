'''
    Celery Task definitions to run HEALER enumerations in the background.
'''
import healer.utils.rdkit_monkey_patch  # noqa: F401

import os

from celery import Celery

from healer.web.interface import run_molecule_enumeration, run_site_enumeration, format_enumeration_results

# Get Redis URL from environment or default to localhost
REDIS_URL = os.environ.get("HEALER_REDIS_URL", "redis://localhost:6379/0")

celery_app = Celery(
    "healer_worker",
    broker=REDIS_URL,
    backend=REDIS_URL
)

celery_app.conf.update(
    task_serializer="json",
    accept_content=["json"],
    result_serializer="json",
    timezone="UTC",
    enable_utc=True,
    task_track_started=True,
)

@celery_app.task(bind=True, name="healer.web.celery_worker.task_enumerate_molecule")
def task_enumerate_molecule(self, params: dict):
    """
    Celery task to run molecule enumeration.
    params: Dictionary matching MoleculeRequest model
    """
    try:
        # Run the enumeration
        # Note: JSON serialization converts tuples to lists, but healer code 
        # generally handles list-of-lists for sites fine.
        raw_results = run_molecule_enumeration(**params)
        
        # Format results for the frontend
        display_res, complete_res = format_enumeration_results(raw_results, 'molecule')
        
        return {
            "display": display_res,
            "complete": complete_res
        }
    except Exception as e:
        # Log the error and re-raise so Celery marks task as FAILED
        raise e

@celery_app.task(bind=True, name="healer.web.celery_worker.task_enumerate_site")
def task_enumerate_site(self, params: dict):
    """
    Celery task to run site enumeration.
    params: Dictionary matching SiteRequest model
    """
    try:
        raw_results = run_site_enumeration(**params)
        display_res, complete_res = format_enumeration_results(raw_results, 'site')
        
        return {
            "display": display_res,
            "complete": complete_res
        }
    except Exception as e:
        raise e
