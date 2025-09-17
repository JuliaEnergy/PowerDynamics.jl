#!/bin/bash

# Get script directory
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Set temporary directory
TMPDIR="/tmp/gencls_sim"

# Create temporary directory
mkdir -p "$TMPDIR" && cd "$TMPDIR"

# Clone OpenIPSL at specific version for reproducibility
git clone --depth 1 --branch v3.0.1 https://github.com/OpenIPSL/OpenIPSL.git

# Set OpenIPSL library path
export OPENIPSL_PATH="$TMPDIR/OpenIPSL/OpenIPSL"

# Run OpenModelica simulation
omc "$SCRIPTDIR/gencls_simulation.mos"

# Copy results back, compress, and cleanup
cp modelica_results.csv "$SCRIPTDIR/"
gzip "$SCRIPTDIR/modelica_results.csv"
rm -rf "$TMPDIR"