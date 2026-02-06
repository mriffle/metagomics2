#!/bin/bash
set -e

echo "Starting Metagomics 2 v${METAGOMICS_VERSION:-unknown}"

# Start the background worker
echo "Starting worker..."
python -m metagomics2.worker.worker &
WORKER_PID=$!

# Trap signals to forward to both processes
cleanup() {
    echo "Shutting down..."
    kill -TERM "$WORKER_PID" 2>/dev/null || true
    wait "$WORKER_PID" 2>/dev/null || true
    exit 0
}
trap cleanup SIGTERM SIGINT

# Start the web server (foreground)
echo "Starting web server on port 8000..."
exec uvicorn metagomics2.server.app:app --host 0.0.0.0 --port 8000
