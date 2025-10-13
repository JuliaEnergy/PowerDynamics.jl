#!/bin/bash

# Get script directory
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Set temporary directory
TMPDIR="/tmp/ggov1_sim"

# Create temporary directory
mkdir -p "$TMPDIR" && cd "$TMPDIR"

# Remove existing OpenIPSL directory if it exists
rm -rf OpenIPSL

# Clone OpenIPSL at specific commit for reproducibility
git clone --depth 1 https://github.com/OpenIPSL/OpenIPSL.git
cd OpenIPSL
git fetch --depth 1 origin fe8aa5c
git checkout fe8aa5c
cd ..

# Set OpenIPSL library path
export OPENIPSL_PATH="$TMPDIR/OpenIPSL/OpenIPSL"

# Apply patch to bypass delay completely by directly connecting gain1.y to s4.u
cd "$TMPDIR/OpenIPSL"
cat << 'EOF' | git apply --ignore-whitespace
diff --git a/OpenIPSL/Electrical/Controls/PSSE/TG/BaseClasses/GGOV1/Turbine.mo b/OpenIPSL/Electrical/Controls/PSSE/TG/BaseClasses/GGOV1/Turbine.mo
index 1234567..abcdefg 100644
--- a/OpenIPSL/Electrical/Controls/PSSE/TG/BaseClasses/GGOV1/Turbine.mo
+++ b/OpenIPSL/Electrical/Controls/PSSE/TG/BaseClasses/GGOV1/Turbine.mo
@@ -141,9 +141,7 @@ equation
   connect(SPEED, Dw2w.u1) annotation (
     Line(points = {{-220, 80}, {-130, 80}, {-130, 26}, {-122, 26}}, color = {0, 0, 127}));
   connect(FSR, add8.u1) annotation (
     Line(points = {{-220, 0}, {-180, 0}, {-180, -34}, {-172, -34}, {-172, -34}}, color = {0, 0, 127}));
-  connect(fixedDelay.y, s4.u) annotation (Line(points={{120.6,-26},{130,-26},{130,-40},{138,-40}}, color={0,0,127}));
-  connect(padeDelay.y, s4.u) annotation (Line(points={{120.6,-54},{130,-54},{130,-40},{138,-40}}, color={0,0,127},
-      pattern=LinePattern.Dash));
+  connect(gain1.y, s4.u) annotation (Line(points={{89,-40},{138,-40}}, color={0,0,127}));
   connect(gain1.y, fixedDelay.u) annotation (Line(points={{89,-40},{98,-40},{98,-26},{106.8,-26}}, color={0,0,127}));
   connect(gain1.y, padeDelay.u) annotation (Line(points={{89,-40},{98,-40},{98,-54},{106.8,-54}}, color={0,0,127},
       pattern=LinePattern.Dash));
EOF
cd ..

# Run OpenModelica simulation
omc "$SCRIPTDIR/model_simulation.mos"

# Copy results back, compress, and cleanup
# Copy minimal data (for git repo)
cp modelica_results_modified.csv "$SCRIPTDIR/" || { echo "Error: Failed to copy simulation results"; exit 1; }

# Remove existing compressed files if they exist to avoid conflicts
rm -f "$SCRIPTDIR/modelica_results_modified.csv.gz"

# Compress the file
gzip "$SCRIPTDIR/modelica_results_modified.csv"

rm -rf "$TMPDIR"