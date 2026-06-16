#!/bin/sh
set -e

cd "${0%/*}" || exit 1

RECON_START="${RECON_START:-0.01}"
RECON_END="${RECON_END:-0.5}"
RECON_DT="${RECON_DT:-0.01}"
N_MODES="${N_MODES:-5}"

command -v python3 >/dev/null 2>&1 || {
    echo "python3 is required" >&2
    exit 1
}

python3 -c 'import numpy' >/dev/null 2>&1 || {
    echo "Python package numpy is required" >&2
    exit 1
}

command -v mpirun >/dev/null 2>&1 || {
    echo "mpirun is required" >&2
    exit 1
}

echo "Building foamToNumpy and numpyToFoam..."
(cd .. && ./Allwmake)

echo "Cleaning generated tutorial outputs..."
./reset_tutorial.sh

echo "Running cavity simulation..."
(cd of_cavity && ./Allclean && ./Allrun)

NPROC="$(find of_cavity -maxdepth 1 -type d -name 'processor*' | wc -l)"
NPROC="$(printf '%s' "$NPROC" | tr -d '[:space:]')"
if [ "$NPROC" -eq 0 ]; then
    echo "No processor directories found in tutorials/of_cavity" >&2
    exit 1
fi
echo "Detected $NPROC processor directories."

echo "Exporting pressure and velocity snapshots with foamToNumpy..."
(cd of_cavity && mpirun -np "$NPROC" foamToNumpy -parallel)

echo "Computing POD reconstruction and modes..."
python3 compute_pod.py \
    --input of_cavity/exported_data \
    --reconstruction-output reconstruction_data \
    --modes-output mode_data \
    --rank 15 \
    --n-modes "$N_MODES"

echo "Writing reconstructed pressure and velocity snapshots..."
cp -r of_cavity of_cavity_reconstructed
(cd of_cavity_reconstructed && ./Clean_proc_data)
cp -r reconstruction_data of_cavity_reconstructed/reconstruction_data

cat > of_cavity_reconstructed/system/numpyToFoamDict <<EOF
FoamFile
{
    format      ascii;
    class       dictionary;
    object      numpyToFoamDict;
}

dataDir       reconstruction_data;
fields        (p U);

time
{
    startTime   ${RECON_START};
    endTime     ${RECON_END};
    deltaT      ${RECON_DT};
}
EOF

(cd of_cavity_reconstructed && mpirun -np "$NPROC" numpyToFoam -parallel)
touch of_cavity_reconstructed/of_cavity_reconstructed.foam

echo "Writing POD modes as OpenFOAM time directories..."
cp -r of_cavity of_cavity_modes
(cd of_cavity_modes && ./Clean_proc_data)
cp -r mode_data of_cavity_modes/mode_data

cat > of_cavity_modes/system/numpyToFoamDict <<EOF
FoamFile
{
    format      ascii;
    class       dictionary;
    object      numpyToFoamDict;
}

dataDir       mode_data;
fields        (p_modes U_modes);

time
{
    startTime   1;
    endTime     ${N_MODES};
    deltaT      1;
}
EOF

(cd of_cavity_modes && mpirun -np "$NPROC" numpyToFoam -parallel)
touch of_cavity_modes/of_cavity_modes.foam

echo "Tutorial complete."
echo "Reconstructed case: tutorials/of_cavity_reconstructed"
echo "POD modes case:     tutorials/of_cavity_modes"
