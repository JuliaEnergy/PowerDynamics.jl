#!/bin/bash

# Get script directory
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Set temporary directory
TMPDIR="/tmp/genrou_sim"

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
omc "$SCRIPTDIR/genrou_simulation.mos"

# Copy results back, compress, and cleanup
# Copy minimal generator data (for git repo)
cp modelica_results.csv "$SCRIPTDIR/" || { echo "Error: Failed to copy minimal simulation results"; exit 1; }
# Copy extended data with all bus variables (for debugging)
cp modelica_results_extended.csv "$SCRIPTDIR/" || { echo "Error: Failed to copy extended simulation results"; exit 1; }

# Remove existing compressed files if they exist to avoid conflicts
rm -f "$SCRIPTDIR/modelica_results.csv.gz"
rm -f "$SCRIPTDIR/modelica_results_extended.csv.gz"

# Compress both files
gzip "$SCRIPTDIR/modelica_results.csv"
gzip "$SCRIPTDIR/modelica_results_extended.csv"

rm -rf "$TMPDIR"
