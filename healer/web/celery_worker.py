"""
Celery Task definitions to run HEALER enumerations in the background.
"""

import os
from typing import Callable

from celery import Celery

import healer.utils.rdkit_monkey_patch  # noqa: F401
from healer.web.interface import run_enumeration_job

# Get Redis URL from environment or default to localhost
REDIS_URL = os.environ.get("HEALER_REDIS_URL", "redis://localhost:6379/0")

celery_app = Celery("healer_worker", broker=REDIS_URL, backend=REDIS_URL)

celery_app.conf.update(
    task_serializer="json",
    accept_content=["json"],
    result_serializer="json",
    timezone="UTC",
    enable_utc=True,
    task_track_started=True,
    result_expires=7200,  # 2 hours — explicit TTL to prevent Redis bloat
)


def _stage_reporter(task) -> Callable[[str], None]:
    """Build a callback that publishes the current stage as task state."""

    def report(stage: str) -> None:
        task.update_state(state="PROGRESS", meta={"stage": stage})

    return report


@celery_app.task(bind=True, name="healer.web.celery_worker.task_enumerate_molecule")
def task_enumerate_molecule(self, params: dict):
    """
    Celery task to run molecule enumeration.
    params: Dictionary matching MoleculeRequest model
    """
    return run_enumeration_job("molecule", params, on_stage=_stage_reporter(self))


@celery_app.task(bind=True, name="healer.web.celery_worker.task_enumerate_site")
def task_enumerate_site(self, params: dict):
    """
    Celery task to run site enumeration.
    params: Dictionary matching SiteRequest model
    """
    return run_enumeration_job("site", params, on_stage=_stage_reporter(self))
