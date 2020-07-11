cmake_minimum_required(VERSION 3.9)

project(bliss VERSION 0.73 DESCRIPTION "bliss library")


add_library(bliss SHARED
    src/bliss.cc
    src/heap.cc
    src/uintseqhash.cc
    src/bliss_C.cc
    src/orbit.cc
    src/utils.cc
    src/defs.cc
    src/partition.cc
    src/graph.cc
    src/timer.cc
)

set_target_properties(bliss PROPERTIES VERSION ${PROJECT_VERSION})

set_target_properties(bliss PROPERTIES PUBLIC_HEADER include/mylib.h)

