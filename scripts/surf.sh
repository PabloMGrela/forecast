#!/usr/bin/env bash
# Uses VENV_PYTHON env var if set, otherwise falls back to python3
PYTHON="${VENV_PYTHON:-python3}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Forward all arguments: [tomorrow] [HH]
$PYTHON "$SCRIPT_DIR/oceanographic_engine.py" "$@"
