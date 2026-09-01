#!/usr/bin/env bash

POLY="${1:-POLY_SerialDualDesc}"

"$POLY" \
    input/input_G3_1_3_path_G313-P025_k7_d8.nml \
    2>&1 | tee run.log

STATUS=${PIPESTATUS[0]}

echo "exit_status=$STATUS"

if grep -q \
    'Normal termination of POLY_SerialDualDesc' \
    run.log
then
    echo "normal_termination=yes"
else
    echo "normal_termination=no"
fi

echo "FINISHED — WSL REMAINS OPEN"
