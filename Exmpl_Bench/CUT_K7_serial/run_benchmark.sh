#!/bin/bash

rm -rf Saving_Bank && mkdir Saving_Bank
rm -rf Saving_Polyhedral && mkdir Saving_Polyhedral

time (../../src_dualdesc/POLY_SerialDualDesc CUT_K7.nml 2>&1 | tee err2_polydualdesc)
