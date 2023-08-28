#!/bin/bash

rm -rf Saving_Bank && mkdir Saving_Bank
rm -rf Saving_Polyhedral && mkdir Saving_Polyhedral

../../src_dualdesc/POLY_SerialDualDesc serial.nml 2> RUN1
../../src_dualdesc/POLY_SerialDualDesc serial.nml 2> RUN2
