#!/bin/bash

# Get script directory
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Set temporary directory
TMPDIR="/tmp/ieeest_sim"

# Create temporary directory
mkdir -p "$TMPDIR" && cd "$TMPDIR"

# Remove existing OpenIPSL directory if it exists
rm -rf OpenIPSL

# Clone OpenIPSL at specific commit for reproducibility
git clone --depth 1 https://github.com/OpenIPSL/OpenIPSL.git
cd OpenIPSL
git fetch --depth 1 origin fe8aa5c
git checkout fe8aa5c

# Set OpenIPSL library path
export OPENIPSL_PATH="$TMPDIR/OpenIPSL/OpenIPSL"

# Run OpenModelica simulation
omc "$SCRIPTDIR/ieeest_simulation.mos"

# Copy results back, compress, and cleanup
# Copy minimal data (for git repo)
cp modelica_results.csv "$SCRIPTDIR/" || { echo "Error: Failed to copy simulation results"; exit 1; }

# Remove existing compressed files if they exist to avoid conflicts
rm -f "$SCRIPTDIR/modelica_results.csv.gz"

# Compress the file
gzip "$SCRIPTDIR/modelica_results.csv"

rm -rf "$TMPDIR"