#!/bin/bash

time (../../src_dualdesc/POLY_SerialDualDesc CUT_K8.nml 2>&1 | tee err2_restart)
