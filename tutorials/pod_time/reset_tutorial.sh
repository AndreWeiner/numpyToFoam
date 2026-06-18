#!/bin/sh
set -e

cd "${0%/*}" || exit 1

rm -rf reconstruction_data mode_data \
    of_cavity_reconstructed of_cavity_modes \
    of_cavity/exported_data

if [ -x of_cavity/Allclean ]; then
    (cd of_cavity && ./Allclean)
fi
