'''
    FastAPI server entry point for the HEALER web application.
'''
import os
from dotenv import load_dotenv

# Load environment variables from the root directory
root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
env_path = os.path.join(root_dir, ".env")
if os.path.exists(env_path):
    load_dotenv(env_path)

# Check if web dependencies are available
try:
    from fastapi import FastAPI
    from fastapi.staticfiles import StaticFiles
    from fastapi.middleware.cors import CORSMiddleware
    _WEB_AVAILABLE = True
except ImportError:
    _WEB_AVAILABLE = False


def _create_app():
    """Create and configure the FastAPI application."""
    import healer.utils.rdkit_monkey_patch  # noqa: F401
    from healer.web.routes import router

    app = FastAPI(title="HEALER Web API")

    # Enable CORS for development
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["http://localhost:3000", "http://127.0.0.1:3000", "http://localhost:5173"],
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    app.include_router(router)

    @app.get("/api/health")
    async def health_check():
        return {"status": "ok", "message": "HEALER API is running"}

    # Mount static files if they exist
    static_dir = os.path.join(os.path.dirname(__file__), "static")
    if os.path.exists(static_dir) and os.listdir(static_dir):
        app.mount("/", StaticFiles(directory=static_dir, html=True), name="static")

    return app


# Create app only if dependencies are available (for uvicorn import)
app = _create_app() if _WEB_AVAILABLE else None


def start():
    """Entry point for the 'healer-ui' command."""
    if not _WEB_AVAILABLE:
        print("Error: Web dependencies not installed.")
        print("Install with: pip install mol-healer[web]")
        raise SystemExit(1)
    
    import uvicorn
    uvicorn.run("healer.web.app:app", host="0.0.0.0", port=8000, reload=False)


if __name__ == "__main__":
    start()
