#!/usr/bin/env bash
set -euo pipefail

echo "=== gnomAD MT Browser - Setup ==="

# Check Python 3
if ! command -v python3 &>/dev/null; then
    echo "Error: python3 is required but not found."
    exit 1
fi

PYTHON_VERSION=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
echo "Found Python $PYTHON_VERSION"

# Create virtual environment if it doesn't exist
VENV_DIR="$(cd "$(dirname "$0")" && pwd)/.venv"
if [ ! -d "$VENV_DIR" ]; then
    echo "Creating virtual environment in $VENV_DIR..."
    python3 -m venv "$VENV_DIR"
else
    echo "Virtual environment already exists at $VENV_DIR"
fi

# Activate and install
echo "Installing dependencies..."
"$VENV_DIR/bin/pip" install --upgrade pip -q
"$VENV_DIR/bin/pip" install -r "$(dirname "$0")/requirements.txt" -q

echo ""
echo "=== Setup complete ==="
echo ""
echo "To run the browser:"
echo "  source $VENV_DIR/bin/activate"
echo "  python3 app.py"
echo ""
echo "Then open http://127.0.0.1:5001 in your browser."
