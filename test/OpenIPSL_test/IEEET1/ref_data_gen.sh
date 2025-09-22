#!/bin/bash

# Get script directory
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Set temporary directory
TMPDIR="/tmp/ieeet1_sim"

# Create temporary directory
mkdir -p "$TMPDIR" && cd "$TMPDIR"

# Clone OpenIPSL at specific version for reproducibility
git clone --depth 1 --branch v3.0.1 https://github.com/OpenIPSL/OpenIPSL.git

# Set OpenIPSL library path
export OPENIPSL_PATH="$TMPDIR/OpenIPSL/OpenIPSL"

# Run OpenModelica simulation
omc "$SCRIPTDIR/ieeet1_simulation.mos"

# Copy results back, compress, and cleanup
# Copy minimal data (for git repo)
cp modelica_results.csv "$SCRIPTDIR/" || { echo "Error: Failed to copy simulation results"; exit 1; }

# Remove existing compressed files if they exist to avoid conflicts
rm -f "$SCRIPTDIR/modelica_results.csv.gz"

# Compress the file
gzip "$SCRIPTDIR/modelica_results.csv"

rm -rf "$TMPDIR"