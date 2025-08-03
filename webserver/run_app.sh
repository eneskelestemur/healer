#!/bin/bash
# Startup script for the HEALER web application

echo "Starting HEALER Dashboard..."
echo "Access the application at: http://localhost:8053"
echo "Press Ctrl+C to stop the server"
echo ""

cd "$(dirname "$0")"
python app.py
