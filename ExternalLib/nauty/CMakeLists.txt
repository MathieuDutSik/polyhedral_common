cmake_minimum_required(VERSION 3.9)

project(nauty VERSION 1.0.0 DESCRIPTION "libshort library")

add_library(nauty SHARED
    gtools.c
    naurng.c
    nautinv.c
    nautil.c
    nauty.c
    naututil.c
    nausparse.c
    naugraph.c
    schreier.c
    traces.c
)

set_target_properties(nauty PROPERTIES VERSION ${PROJECT_VERSION})
#set_target_properties(nauty PROPERTIES PUBLIC_HEADER la_support.h lll_basis.h lrsgmp.h lrslib.h panic.h shvec.h)


# Install library
install(TARGETS nauty DESTINATION lib)

# Install library headers
file(GLOB HEADERS *.h)
install(FILES ${HEADERS} DESTINATION include)


