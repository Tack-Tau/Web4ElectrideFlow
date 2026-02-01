#!/bin/bash
# 
# Inorganic Electride Database Web Application
# Based on the electride_flow workflow: https://github.com/Tack-Tau/electride_flow
#

set -e  # Exit immediately on error

# 1. Setup JSmol assets (similar to lego-crystal)
TARGET_LINK="ase_root/ase/db/static/jsmol"
mkdir -p "$(dirname "$TARGET_LINK")"

# If a bundled jsmol folder is present, link it
if [ -d "jsmol" ]; then
    if [ ! -e "$TARGET_LINK" ] && [ ! -L "$TARGET_LINK" ]; then
        echo "Linking JSmol to $TARGET_LINK"
        ln -s "$PWD/jsmol" "$TARGET_LINK"
    elif [ -L "$TARGET_LINK" ]; then
        echo "JSmol link already exists"
    fi
else
    echo "WARNING: JSmol directory not found. Structure viewer will be unavailable."
fi

# 2. Determine port (CLI arg > $PORT > default 5000)
if [ -n "$1" ]; then
    PORT="$1"
elif [ -z "$PORT" ]; then
    PORT=5000
fi
export PORT

# 3. Run the Flask app
echo "Starting Electride Database web server on port $PORT..."
export PYTHONPATH="$PWD/ase_root${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONWARNINGS="ignore::FutureWarning"
python3 app.py
