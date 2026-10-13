#!/bin/bash
set -euo pipefail

echo "=========================================================="
echo "  Abaqus: ODB -> VTK export & ZIP                         "
echo "=========================================================="

JOB_NAME="${1:-Gear_Shaft_assembly_simulation}"
ODB_FILE="${JOB_NAME%.odb}.odb"

# 1. Abaqus launcher (override with the ABAQUS_CMD environment variable)
LAUNCHER="${ABAQUS_CMD:-abaqus}"
if ! command -v "${LAUNCHER}" >/dev/null 2>&1; then
    echo "[ERROR] Abaqus launcher '${LAUNCHER}' not found. Set ABAQUS_CMD to the launcher of your installation."
    exit 1
fi

# 2. Run the ODB -> VTK conversion script with the Abaqus Python interpreter
echo "[INFO] Running: ${LAUNCHER} python odb_to_vtk.py ${ODB_FILE}"
"${LAUNCHER}" python odb_to_vtk.py "${ODB_FILE}"

# 3. Pack every VTK file (frames and fatigue cycle summary) into gear_helical_vtk.zip
echo "[INFO] Packing the VTK files into 'gear_helical_vtk.zip'..."
if command -v zip >/dev/null 2>&1; then
    zip -q -9 gear_helical_vtk.zip gear_helical_*.vtk
else
    "${LAUNCHER}" python -c "import zipfile, glob; z = zipfile.ZipFile('gear_helical_vtk.zip', 'w', zipfile.ZIP_DEFLATED); [z.write(f) for f in sorted(glob.glob('gear_helical_*.vtk'))]; z.close()"
fi

echo "=========================================================="
echo "  Post-processing completed: gear_helical_vtk.zip         "
echo "=========================================================="
