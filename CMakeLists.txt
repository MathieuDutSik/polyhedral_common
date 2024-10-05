cmake_minimum_required(VERSION 3.30)
project(MyProject)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen3 (make sure you have Eigen installed and CMake can find it)
find_package(Eigen3 3.3 REQUIRED)

# Find Boost (adjust components if you need specific ones like filesystem, system, etc.)
find_package(Boost REQUIRED COMPONENTS system filesystem)

# Add include directories (adjust paths to your actual include locations)
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/basic_common_cpp/src_basic)
include_directories(${CMAKE_SOURCE_DIR}/basic_common_cpp/src_number)
include_directories(${CMAKE_SOURCE_DIR}/basic_common_cpp/src_matrix)
include_directories(${CMAKE_SOURCE_DIR}/basic_common_cpp/src_comb)
include_directories(${CMAKE_SOURCE_DIR}/basic_common_cpp/src_graph)
include_directories(${CMAKE_SOURCE_DIR}/src_copos)
include_directories(${CMAKE_SOURCE_DIR}/src_dualdesc)
include_directories(${CMAKE_SOURCE_DIR}/src_group)
include_directories(${CMAKE_SOURCE_DIR}/src_indefinite)
include_directories(${CMAKE_SOURCE_DIR}/src_isotropy)
include_directories(${CMAKE_SOURCE_DIR}/src_latt)
include_directories(${CMAKE_SOURCE_DIR}/src_lorentzian)
include_directories(${CMAKE_SOURCE_DIR}/src_poly)
include_directories(${CMAKE_SOURCE_DIR}/permutalib/src)
include_directories(${CMAKE_SOURCE_DIR}/src)

include(ExternalProject)

ExternalProject_Add(
    NautyGitPackage
    PREFIX "${CMAKE_BINARY_DIR}/external"
    URL "https://pallini.di.uniroma1.it/nauty2_8_9.tar.gz"
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/external
    BUILD_IN_SOURCE 1
)

# Add the source codes.
file(GLOB_RECURSE SOURCES ${CMAKE_SOURCE_DIR}/src_copos/CP_TestCopositivity.cpp)
file(GLOB_RECURSE SOURCES ${CMAKE_SOURCE_DIR}/src_copos/CP_TestCompletePositivity.cpp)
file(GLOB_RECURSE SOURCES ${CMAKE_SOURCE_DIR}/src_dualdesc/POLY_SerialDualDesc.cpp)

# Add the executable
add_executable(POLY_SerialDualDesc ${SOURCES})
add_executable(CP_TestCopositivity ${SOURCES})
add_executable(CP_TestCompletePositivity ${SOURCES})

# Add dependency to your target
add_dependencies(POLY_SerialDualDesc NautyGitPackage)
