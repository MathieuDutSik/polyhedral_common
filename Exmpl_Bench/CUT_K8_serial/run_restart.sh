#!/bin/bash

rm -rf Saving_Bank && mkdir Saving_Bank
rm -rf Saving_Polyhedral && mkdir Saving_Polyhedral

# Normally two blocks should suffice. max_runtime is 1100 while
# usually, it runs in less than 2000 seconds.

time (../../src_dualdesc/POLY_SerialDualDesc CUT_K8_limit.nml 2>&1 | tee err2_part1)
time (../../src_dualdesc/POLY_SerialDualDesc CUT_K8_limit.nml 2>&1 | tee err2_part2)
