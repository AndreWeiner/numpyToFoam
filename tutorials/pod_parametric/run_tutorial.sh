#!/bin/sh
set -e

cd "${0%/*}" || exit 1

(cd ../.. && ./Allwmake)
./reset_tutorial.sh
python3 parametric_airfoil.py
