cmake_minimum_required(VERSION 3.9)

project(libshort VERSION 0.73 DESCRIPTION "libshort library")

add_library(libshort SHARED
    shvec.c
    panic.c
    lrslib.c
    lrsgmp.c
    lll_basis.c
    la_support.c
)

set_target_properties(libshort PROPERTIES VERSION ${PROJECT_VERSION})

#set_target_properties(libshort PROPERTIES PUBLIC_HEADER la_support.h lll_basis.h lrsgmp.h lrslib.h panic.h shvec.h)


# Install library
install(TARGETS libshort DESTINATION lib)

# Install library headers
file(GLOB HEADERS *.h)
install(FILES ${HEADERS} DESTINATION include)


