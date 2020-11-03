cmake_minimum_required(VERSION 3.9)

project(short VERSION 1.0.0 DESCRIPTION "libshort library")

add_library(short SHARED
    shvec.c
    panic.c
    lrslib.c
    lrsgmp.c
    lll_basis.c
    la_support.c
)

set_target_properties(short PROPERTIES VERSION ${PROJECT_VERSION})
target_include_directories(short PUBLIC /usr/local/Cellar/gmp/6.2.0/include)
target_link_libraries(short PUBLIC -L/usr/local/Cellar/gmp/6.2.0/lib -lgmp)
#set_target_properties(short PROPERTIES PUBLIC_HEADER la_support.h lll_basis.h lrsgmp.h lrslib.h panic.h shvec.h)


# Install library
install(TARGETS short DESTINATION lib)

# Install library headers
file(GLOB HEADERS *.h)
install(FILES ${HEADERS} DESTINATION include)


