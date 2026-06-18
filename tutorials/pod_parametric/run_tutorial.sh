#!/bin/sh
set -e

cd "${0%/*}" || exit 1

N_CASES="${N_CASES:-8}"
POD_SPLINE_SMOOTHING="${POD_SPLINE_SMOOTHING:-0}"
POD_MODES="${POD_MODES:-0}"

command -v python3 >/dev/null 2>&1 || {
    echo "python3 is required" >&2
    exit 1
}

python3 -c 'import numpy' >/dev/null 2>&1 || {
    echo "Python package numpy is required" >&2
    exit 1
}

python3 -c 'import matplotlib' >/dev/null 2>&1 || {
    echo "Python package matplotlib is required for pressure-distribution plots" >&2
    exit 1
}

command -v mpirun >/dev/null 2>&1 || {
    echo "mpirun is required" >&2
    exit 1
}

echo "Building foamToNumpy and numpyToFoam..."
(cd ../.. && ./Allwmake)

echo "Cleaning generated tutorial outputs..."
./reset_tutorial.sh

python3 parametric_airfoil.py \
    --cases "$N_CASES" \
    --spline-smoothing "$POD_SPLINE_SMOOTHING" \
    --modes "$POD_MODES"
