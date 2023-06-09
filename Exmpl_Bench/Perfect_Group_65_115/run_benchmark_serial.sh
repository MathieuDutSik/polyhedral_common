#!/bin/bash

#rm -rf Saving_Bank && mkdir Saving_Bank
#rm -rf Saving_Polyhedral && mkdir Saving_Polyhedral

/usr/bin/time --verbose ../../src_dualdesc/POLY_SerialDualDesc mpi_serial.nml 2>&1 | tee err2_polydualdesc
