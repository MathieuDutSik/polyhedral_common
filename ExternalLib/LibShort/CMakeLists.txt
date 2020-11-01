cmake_minimum_required(VERSION 3.9)

project(libshort VERSION 0.73 DESCRIPTION "libshort library")

add_library(libshort SHARED
    src/shvec.c
    src/panic.c
    src/lrslib.c
    src/lrsgmp.c
    src/lll_basis.c
    src/la_support.c
)

set_target_properties(libshort PROPERTIES VERSION ${PROJECT_VERSION})

#set_target_properties(libshort PROPERTIES PUBLIC_HEADER src/la_support.h src/lll_basis.h src/lrsgmp.h src/lrslib.h src/panic.h src/shvec.h)


# Install library
install(TARGETS libshort DESTINATION lib)

# Install library headers
file(GLOB HEADERS src/*.h)
install(FILES ${HEADERS} DESTINATION include)


