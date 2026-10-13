#!/bin/bash
# Serve the Chemical LR Resource site on the local network.
#   ./serve.sh [port]      (default 8080)
# Bound to 0.0.0.0 so anyone on the same network can open it.
PORT="${1:-8080}"
cd "$(dirname "$0")"
IP=$(hostname -I | awk '{print $1}')
echo "Chemical LR Resource site:"
echo "  this machine : http://localhost:$PORT"
echo "  shareable    : http://$IP:$PORT"
echo "Ctrl-C to stop."
exec python3 -m http.server "$PORT" --bind 0.0.0.0
