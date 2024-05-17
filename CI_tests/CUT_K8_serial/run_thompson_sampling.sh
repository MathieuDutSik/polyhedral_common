#!/bin/bash

rm -rf Saving_Bank && mkdir Saving_Bank
rm -rf Saving_Polyhedral && mkdir Saving_Polyhedral

time (../../src_dualdesc/POLY_SerialDualDesc CUT_K8_ts.nml 2>&1 | tee err2_polydualdesc_ts)
