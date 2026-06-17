#!/bin/sh
set -e

cd "${0%/*}" || exit 1

rm -rf results prediction_data of_airfoil_podi
rm -rf of_airfoil/exported_data of_airfoil/postProcessing

if [ -x of_airfoil/Allclean ]; then
    (cd of_airfoil && ./Allclean)
fi
