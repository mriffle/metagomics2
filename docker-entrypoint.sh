#!/bin/bash
set -e

echo "Starting Metagomics 2 v${METAGOMICS_VERSION:-unknown}"

 if [ "$#" -gt 0 ]; then
     case "$1" in
         metagomics2)
             shift
             exec metagomics2 "$@"
             ;;
         run|version|-*)
             exec metagomics2 "$@"
             ;;
         *)
             exec "$@"
             ;;
     esac
 fi

# Run the worker and the web server as sibling processes and supervise both.
# If either one dies, stop the other and exit non-zero so Docker's restart
# policy restarts the container, instead of silently running without a worker.

echo "Starting worker..."
python -m metagomics2.worker.worker &
WORKER_PID=$!

echo "Starting web server on port 8000..."
uvicorn metagomics2.server.app:app --host 0.0.0.0 --port 8000 &
SERVER_PID=$!

shutdown() {
    echo "Shutting down..."
    kill -TERM "$WORKER_PID" "$SERVER_PID" 2>/dev/null || true
    wait "$WORKER_PID" "$SERVER_PID" 2>/dev/null || true
    exit 0
}
trap shutdown SIGTERM SIGINT

# Wait for whichever process exits first.
set +e
wait -n "$WORKER_PID" "$SERVER_PID"
EXIT_CODE=$?
set -e

if kill -0 "$WORKER_PID" 2>/dev/null; then
    echo "FATAL: web server exited with code ${EXIT_CODE}; stopping worker" >&2
    kill -TERM "$WORKER_PID" 2>/dev/null || true
else
    echo "FATAL: worker exited with code ${EXIT_CODE}; stopping web server so the container restarts" >&2
    echo "       (an exit code of 137 or a missing 'Worker stopped' line usually means it was killed for memory)" >&2
    kill -TERM "$SERVER_PID" 2>/dev/null || true
fi
wait "$WORKER_PID" "$SERVER_PID" 2>/dev/null || true
exit 1
