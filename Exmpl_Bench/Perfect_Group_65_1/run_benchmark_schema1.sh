#!/bin/bash

rm -rf Saving_Bank && mkdir Saving_Bank
rm -rf Saving_Polyhedral && mkdir Saving_Polyhedral

time (../../src_dualdesc/POLY_SerialDualDesc Perfect_schema1_stepA.nml 2>&1 | tee err2_benchmark_schema1_stepA)


time (../../src_dualdesc/POLY_SerialDualDesc Perfect_canonic_initial_triv.nml 2>&1 | tee err2_benchmark_schema1_stepB)
