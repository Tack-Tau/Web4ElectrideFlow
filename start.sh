#!/bin/bash
# 
# Inorganic Electride Database Web Application
# Based on the electride_flow workflow: https://github.com/Tack-Tau/electride_flow
#

set -e  # Exit immediately on error

# 1. Setup JSmol assets (copy instead of symlink for Render compatibility)
TARGET_DIR="ase_root/ase/db/static/jsmol"
mkdir -p "$(dirname "$TARGET_DIR")"

# If a bundled jsmol folder is present, copy it
if [ -d "jsmol" ]; then
    # Remove any existing symlink or directory
    if [ -e "$TARGET_DIR" ] || [ -L "$TARGET_DIR" ]; then
        echo "Removing existing JSmol symlink/directory..."
        rm -rf "$TARGET_DIR"
    fi
    
    echo "Copying JSmol files to $TARGET_DIR..."
    cp -r jsmol "$TARGET_DIR"
    echo "JSmol files copied successfully"
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
