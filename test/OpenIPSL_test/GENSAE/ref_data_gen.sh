#!/bin/bash

# Get script directory
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Set temporary directory
TMPDIR="/tmp/gensae_sim"

# Create temporary directory
mkdir -p "$TMPDIR" && cd "$TMPDIR"

# Clone OpenIPSL at specific version for reproducibility
git clone --depth 1 --branch v3.0.1 https://github.com/OpenIPSL/OpenIPSL.git

# Set OpenIPSL library path
export OPENIPSL_PATH="$TMPDIR/OpenIPSL/OpenIPSL"

# Run OpenModelica simulation
omc "$SCRIPTDIR/gensae_simulation.mos"

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