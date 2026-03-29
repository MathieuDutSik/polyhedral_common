Developer Guide
===============

Tests
-----

The directory `Exmpl_Bench` contains a bunch of test and development.
They are typically work in progress and not necessarily in a state of being
finished.

The directory `CI_tests` contains some tests that are run in CI on GitHub.
Their runtime should be short, from 5 minutes to 1 hour. Together they
should cover as much as possible of the functionality of the code. If the
test is working like normally, then it should be scheduled in the cron to
run once per month. We do not want to overflow the credit that we have.


Compilation options related to debug
------------------------------------

There are several environment variables that can be used during runtime.
* `DEBUG` for making some print statements and making some checks. Those checks have to be fast.
* `KEY_VALUE` for printing some `KEY(....) VALUE=(....)` that can be used for postprocessing of the options and heuristic optimization.
* `TIMINGS` for printing some runtime information.
* `SANITY_CHECK` for doing some checks and stopping if incoherence are detected.
* `SANITY_CHECK_EXTENSIVE` for doing some long running checks and stopping if incoherence are detected. Since those are long running, a print statement must mark their beginning and another one their end.
* `TRACK_INFO` for printing stuff for further work down the line.
* `METHOD_COMPARISON` for considering two (or more) different methods and comparing their performance (mainly runtime but not only).

The options `TIMINGS` and `DEBUG` enable all the timings and debugging statements.
For a more granular debugging, stuff like `DEBUG_LINEAR_PROGRAM` can be used. See
the top of the header files.

A printout to `std::cerr` should occur if an error has been identified and the program
will terminate with a call to `TerminalException`. Other print statement should be
encapsulated in `std::ostream & os` that should be passed by reference from the initial
call. So typically for serial output we pass the `std::cerr` while for the parallel runs
we pass a stream to the output of that process. That way we avoid mixing between
different sources.
