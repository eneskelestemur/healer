"""Tests for the Celery task definitions."""

import pytest

celery = pytest.importorskip("celery", reason="celery is an optional dependency")

from healer.web import celery_worker  # noqa: E402


class TestConfiguration:
    def test_the_app_is_named_and_wired_to_a_broker(self):
        assert celery_worker.celery_app.main == "healer_worker"
        assert celery_worker.celery_app.conf.broker_url

    def test_json_serialization_is_enforced(self):
        conf = celery_worker.celery_app.conf
        assert conf.task_serializer == "json"
        assert conf.result_serializer == "json"
        assert conf.accept_content == ["json"]

    def test_results_expire_so_redis_does_not_grow_without_bound(self):
        assert celery_worker.celery_app.conf.result_expires > 0

    def test_the_broker_url_follows_the_environment(self, monkeypatch):
        import importlib

        monkeypatch.setenv("HEALER_REDIS_URL", "redis://example.test:6379/1")
        reloaded = importlib.reload(celery_worker)
        assert reloaded.REDIS_URL == "redis://example.test:6379/1"

        monkeypatch.delenv("HEALER_REDIS_URL")
        importlib.reload(celery_worker)


class TestTaskRegistration:
    @pytest.mark.parametrize(
        "name",
        [
            "healer.web.celery_worker.task_enumerate_molecule",
            "healer.web.celery_worker.task_enumerate_site",
        ],
    )
    def test_tasks_are_registered_under_their_documented_names(self, name):
        assert name in celery_worker.celery_app.tasks

    def test_callers_pass_only_the_request_payload(self):
        """bind=True absorbs `self`, so callers supply just the params dict."""
        import inspect

        task = celery_worker.celery_app.tasks[
            "healer.web.celery_worker.task_enumerate_molecule"
        ]
        assert list(inspect.signature(task.run).parameters) == ["params"]
