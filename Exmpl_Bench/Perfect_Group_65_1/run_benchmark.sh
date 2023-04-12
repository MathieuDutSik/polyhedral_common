#!/bin/bash

rm -rf Saving_Bank && mkdir Saving_Bank
rm -rf Saving_Polyhedral && mkdir Saving_Polyhedral

time (./POLY_SerialDualDesc Perfect.nml 2>&1 | tee err2_benchmark)
