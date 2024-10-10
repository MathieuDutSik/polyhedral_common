cmake_minimum_required(VERSION 3.30)
project(MyProject)
find_package(PkgConfig REQUIRED)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen3 (make sure you have Eigen installed and CMake can find it)
find_package(Eigen3 3.3 REQUIRED)

# Find Boost (adjust components if you need specific ones like filesystem, system, etc.)
find_package(Boost REQUIRED COMPONENTS serialization)

# Find the GMP library (add REQUIRED to stop if it's not found)
pkg_check_modules(GMP REQUIRED gmpxx)
find_library(GMP_LIBRARY NAMES gmp HINTS ${GMP_LIBRARY_DIRS})
find_library(GMPXX_LIBRARY NAMES gmpxx HINTS ${GMP_LIBRARY_DIRS})
message(STATUS "GMPXX library names: ${GMPXX_LIBRARY}")

# Add include directories (adjust paths to your actual include locations)
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/basic_common_cpp/src_basic)
include_directories(${CMAKE_SOURCE_DIR}/basic_common_cpp/src_number)
include_directories(${CMAKE_SOURCE_DIR}/basic_common_cpp/src_matrix)
include_directories(${CMAKE_SOURCE_DIR}/basic_common_cpp/src_comb)
include_directories(${CMAKE_SOURCE_DIR}/basic_common_cpp/src_graph)
include_directories(${CMAKE_SOURCE_DIR}/basic_common_cpp/sparse-map/include/tsl)
include_directories(${CMAKE_SOURCE_DIR}/basic_common_cpp/robin-map/include/tsl)
include_directories(${CMAKE_SOURCE_DIR}/basic_common_cpp/hopscotch-map/include/tsl)
include_directories(${CMAKE_SOURCE_DIR}/src_copos)
include_directories(${CMAKE_SOURCE_DIR}/src_dualdesc)
include_directories(${CMAKE_SOURCE_DIR}/src_group)
include_directories(${CMAKE_SOURCE_DIR}/src_indefinite)
include_directories(${CMAKE_SOURCE_DIR}/src_isotropy)
include_directories(${CMAKE_SOURCE_DIR}/src_perfect)
include_directories(${CMAKE_SOURCE_DIR}/src_short)
include_directories(${CMAKE_SOURCE_DIR}/src_latt)
include_directories(${CMAKE_SOURCE_DIR}/src_lorentzian)
include_directories(${CMAKE_SOURCE_DIR}/src_poly)
include_directories(${CMAKE_SOURCE_DIR}/permutalib/src)
include_directories(${CMAKE_SOURCE_DIR}/src)

include(ExternalProject)

ExternalProject_Add(
    NautyDirectPackage
    PREFIX "${CMAKE_BINARY_DIR}/external"
    GIT_REPOSITORY https://github.com/MathieuDutSik/nauty
    GIT_TAG main # or the specific commit, tag, or branch you want to use
    UPDATE_DISCONNECTED TRUE
    CONFIGURE_COMMAND ./configure --prefix=${CMAKE_BINARY_DIR}/external/nauty/install
    BUILD_IN_SOURCE 1
)

# Add the executable
add_executable(POLY_SerialDualDesc ${CMAKE_SOURCE_DIR}/src_dualdesc/POLY_SerialDualDesc.cpp)
add_executable(CP_TestCopositivity ${CMAKE_SOURCE_DIR}/src_copos/CP_TestCopositivity.cpp)
add_executable(CP_TestCompletePositivity ${CMAKE_SOURCE_DIR}/src_copos/CP_TestCompletePositivity.cpp)

# Grouped the libraries
set(MY_LIBRARIES_GRP ${CMAKE_BINARY_DIR}/external/nauty/install/lib/libnauty.a ${Boost_LIBRARIES} ${GMPXX_LIBRARY} ${GMP_LIBRARY})
set(MY_LIBRARIES ${Boost_LIBRARIES} ${GMPXX_LIBRARY} ${GMP_LIBRARY})


# Add dependency to your target
add_dependencies(POLY_SerialDualDesc NautyDirectPackage)
target_link_libraries(POLY_SerialDualDesc ${MY_LIBRARIES_GRP})
target_link_libraries(CP_TestCopositivity ${MY_LIBRARIES})
target_link_libraries(CP_TestCompletePositivity ${MY_LIBRARIES})

target_include_directories(POLY_SerialDualDesc PUBLIC ${GMP_INCLUDE_DIRS})
