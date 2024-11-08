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
    GIT_TAG 2.8.9
    UPDATE_DISCONNECTED TRUE
    CONFIGURE_COMMAND ./configure --prefix=${CMAKE_BINARY_DIR}/external/nauty/install
    BUILD_IN_SOURCE 1
)

# Grouped the libraries
set(MY_LIBRARIES_GRP ${CMAKE_BINARY_DIR}/external/nauty/install/lib/libnauty.a ${Boost_LIBRARIES} ${GMPXX_LIBRARY} ${GMP_LIBRARY})
set(MY_LIBRARIES ${Boost_LIBRARIES} ${GMPXX_LIBRARY} ${GMP_LIBRARY})

# Add the executable
add_executable(POLY_DirectSerialDualDesc ${CMAKE_SOURCE_DIR}/src_dualdesc/POLY_DirectSerialDualDesc.cpp)
add_executable(CP_TestCopositivity ${CMAKE_SOURCE_DIR}/src_copos/CP_TestCopositivity.cpp)
add_executable(CP_TestCompletePositivity ${CMAKE_SOURCE_DIR}/src_copos/CP_TestCompletePositivity.cpp)
add_executable(LORENTZ_ReflectiveEdgewalk ${CMAKE_SOURCE_DIR}/src_lorentzian/LORENTZ_ReflectiveEdgewalk.cpp)
add_executable(POLY_DirectFaceLattice ${CMAKE_SOURCE_DIR}/src_poly/POLY_DirectFaceLattice.cpp)
add_executable(INDEF_FORM_AutomorphismGroup ${CMAKE_SOURCE_DIR}/src_indefinite/INDEF_FORM_AutomorphismGroup.cpp)
add_executable(INDEF_FORM_TestEquivalence ${CMAKE_SOURCE_DIR}/src_indefinite/INDEF_FORM_TestEquivalence.cpp)
add_executable(INDEF_FORM_GetOrbitRepresentative ${CMAKE_SOURCE_DIR}/src_indefinite/INDEF_FORM_GetOrbitRepresentative.cpp)
add_executable(INDEF_FORM_GetOrbit_IsotropicKplane ${CMAKE_SOURCE_DIR}/src_indefinite/INDEF_FORM_GetOrbit_IsotropicKplane.cpp)
add_executable(LATT_Canonicalize ${CMAKE_SOURCE_DIR}/src_latt/LATT_Canonicalize.cpp)
add_executable(LATT_FindIsotropic ${CMAKE_SOURCE_DIR}/src_isotropy/LATT_FindIsotropic.cpp)
add_executable(LATT_SerialComputeDelaunay ${CMAKE_SOURCE_DIR}/src_latt/LATT_SerialComputeDelaunay.cpp)
add_executable(LATT_SerialLattice_IsoDelaunayDomain ${CMAKE_SOURCE_DIR}/src_latt/LATT_SerialLattice_IsoDelaunayDomain.cpp)

# ---Delaunay polytopes enumeration (Serial)
# ---Iso-Delaunay spacse (Serial)

# Add dependency to your target
add_dependencies(POLY_SerialDualDesc NautyDirectPackage)
target_link_libraries(POLY_SerialDualDesc ${MY_LIBRARIES_GRP})
target_link_libraries(CP_TestCopositivity ${MY_LIBRARIES})
target_link_libraries(CP_TestCompletePositivity ${MY_LIBRARIES})
target_link_libraries(LORENTZ_FundDomain_AllcockEdgewalk ${MY_LIBRARIES_GRP})
target_link_libraries(POLY_FaceLatticeGen ${MY_LIBRARIES_GRP})
target_link_libraries(INDEF_FORM_AutomorphismGroup ${MY_LIBRARIES_GRP})
target_link_libraries(INDEF_FORM_TestEquivalence ${MY_LIBRARIES_GRP})
target_link_libraries(INDEF_FORM_GetOrbitRepresentative ${MY_LIBRARIES_GRP})
target_link_libraries(INDEF_FORM_GetOrbit_IsotropicKplane ${MY_LIBRARIES_GRP})
target_link_libraries(LATT_Canonicalize ${MY_LIBRARIES_GRP})
target_link_libraries(LATT_FindIsotropic ${MY_LIBRARIES})

target_include_directories(POLY_SerialDualDesc PUBLIC ${GMP_INCLUDE_DIRS})
