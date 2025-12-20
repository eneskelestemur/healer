'''
    FastAPI server entry point for the HEALER web application.
'''
import healer.utils.rdkit_monkey_patch  # noqa: F401

import os

from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware

from healer.web.routes import router

# Create the FastAPI app
app = FastAPI(title="HEALER Web API")

# Enable CORS for development (allows frontend on port 3000 to talk to backend on 8000)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000", "http://127.0.0.1:3000"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include the API routes
app.include_router(router)

# API Routes
@app.get("/api/health")
async def health_check():
    return {"status": "ok", "message": "HEALER API is running"}

# Mount static files (The compiled React app)
# We only mount this if the directory exists (i.e., in production/installed mode)
static_dir = os.path.join(os.path.dirname(__file__), "static")
if os.path.exists(static_dir) and os.listdir(static_dir):
    app.mount("/", StaticFiles(directory=static_dir, html=True), name="static")

def start():
    """Entry point for the 'healer ui' command"""
    import uvicorn
    # In development, we want reload=True, but for the installed CLI tool, defaults are fine
    uvicorn.run("healer.web.app:app", host="0.0.0.0", port=8000, reload=False)

if __name__ == "__main__":
    start()
