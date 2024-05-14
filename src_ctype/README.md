Iso-edge domains
================

This is the set of functionality for dealing with iso-edge
domains.

Enumeration of iso-edge domains
-------------------------------

The program **CTYP_MPI_Enumeration_c** is used for enumerating
the iso-edge domains in dimension six. The result was *55083358*
types. It was done in 20 processors in one week of computation.
The result is written in netCDF files.

Another enumeration program
---------------------------

Another program that should do the same enumeration but with
a more general program is **CTYP_MPI_AdjScheme**.

Finding iso-edge domains without free vectors
---------------------------------------------

The program **CTYP_LookForNoFreeVector** allows to find iso-edge
domains in dimension **n >= 6** by using random walks.
