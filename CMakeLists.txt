cmake_minimum_required(VERSION 3.30)

project(polyhedral_common LANGUAGES CXX)

include(ExternalProject)
include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

find_package(PkgConfig REQUIRED)
find_package(Threads REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(Boost REQUIRED COMPONENTS mpi serialization)
find_package(MPI REQUIRED COMPONENTS CXX)

pkg_check_modules(GMPXX REQUIRED IMPORTED_TARGET gmpxx)
pkg_check_modules(CDDGMP REQUIRED cddgmp)
pkg_check_modules(CDDLIB REQUIRED cddlib)

find_library(GMP_LIBRARY NAMES gmp REQUIRED)
find_path(GLPK_INCLUDE_DIR NAMES glpk.h REQUIRED)
find_library(GLPK_LIBRARY NAMES glpk REQUIRED)
find_path(
  CDD_HEADER_DIR
  NAMES setoper.h cdd.h
  HINTS
    ${CDDGMP_INCLUDE_DIRS}
    ${CDDLIB_INCLUDE_DIRS}
  PATH_SUFFIXES
    cddlib
  REQUIRED
)

find_path(
  NETCDF_INCLUDE_DIR
  NAMES netcdf
  HINTS
    /opt/homebrew/include
    /usr/local/include
)
find_library(
  NETCDF_LIBRARY
  NAMES netcdf
  HINTS
    /opt/homebrew/lib
    /usr/local/lib
)
find_library(
  NETCDF_CXX_LIBRARY
  NAMES netcdf-cxx4
  HINTS
    /opt/homebrew/lib
    /usr/local/lib
)
if(NOT NETCDF_INCLUDE_DIR OR NOT NETCDF_LIBRARY OR NOT NETCDF_CXX_LIBRARY)
  message(
    FATAL_ERROR
      "netCDF C/C++ libraries were not found. Set NETCDF_INCLUDE_DIR, "
      "NETCDF_LIBRARY and NETCDF_CXX_LIBRARY."
  )
endif()

set(NAUTY_PREFIX "${CMAKE_BINARY_DIR}/external/nauty")
set(NAUTY_INSTALL_DIR "${NAUTY_PREFIX}/install")
set(NAUTY_INCLUDE_DIR "${NAUTY_INSTALL_DIR}/include")
set(NAUTY_LIBRARY "${NAUTY_INSTALL_DIR}/lib/libnauty.a")

ExternalProject_Add(
  NautyDirectPackage
  PREFIX "${CMAKE_BINARY_DIR}/external"
  GIT_REPOSITORY https://github.com/MathieuDutSik/nauty
  GIT_TAG 2.8.9
  UPDATE_DISCONNECTED TRUE
  CONFIGURE_COMMAND ./configure --prefix=${NAUTY_INSTALL_DIR}
  BUILD_IN_SOURCE 1
  BUILD_BYPRODUCTS "${NAUTY_LIBRARY}"
)

add_library(polyhedral_nauty STATIC IMPORTED GLOBAL)
set_target_properties(
  polyhedral_nauty
  PROPERTIES
    IMPORTED_LOCATION "${NAUTY_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${NAUTY_INCLUDE_DIR}"
)
add_dependencies(polyhedral_nauty NautyDirectPackage)

add_library(polyhedral_common INTERFACE)
target_include_directories(
  polyhedral_common
  INTERFACE
    "${CMAKE_SOURCE_DIR}"
    "${CMAKE_SOURCE_DIR}/basic_common_cpp"
    "${CMAKE_SOURCE_DIR}/basic_common_cpp/src_basic"
    "${CMAKE_SOURCE_DIR}/basic_common_cpp/src_number"
    "${CMAKE_SOURCE_DIR}/basic_common_cpp/src_matrix"
    "${CMAKE_SOURCE_DIR}/basic_common_cpp/src_comb"
    "${CMAKE_SOURCE_DIR}/basic_common_cpp/src_graph"
    "${CMAKE_SOURCE_DIR}/basic_common_cpp/sparse-map/include/tsl"
    "${CMAKE_SOURCE_DIR}/basic_common_cpp/robin-map/include/tsl"
    "${CMAKE_SOURCE_DIR}/basic_common_cpp/hopscotch-map/include/tsl"
    "${CMAKE_SOURCE_DIR}/permutalib/src"
    "${CMAKE_SOURCE_DIR}/src_basic"
    "${CMAKE_SOURCE_DIR}/src_number"
    "${CMAKE_SOURCE_DIR}/src_matrix"
    "${CMAKE_SOURCE_DIR}/src_comb"
    "${CMAKE_SOURCE_DIR}/src_graph"
    "${CMAKE_SOURCE_DIR}/src_copos"
    "${CMAKE_SOURCE_DIR}/src_ctype"
    "${CMAKE_SOURCE_DIR}/src_delaunay"
    "${CMAKE_SOURCE_DIR}/src_dualdesc"
    "${CMAKE_SOURCE_DIR}/src_enum_schemes"
    "${CMAKE_SOURCE_DIR}/src_group"
    "${CMAKE_SOURCE_DIR}/src_indefinite"
    "${CMAKE_SOURCE_DIR}/src_isotropy"
    "${CMAKE_SOURCE_DIR}/src_latt"
    "${CMAKE_SOURCE_DIR}/src_lorentzian"
    "${CMAKE_SOURCE_DIR}/src_perfect"
    "${CMAKE_SOURCE_DIR}/src_poincare_polyhedron"
    "${CMAKE_SOURCE_DIR}/src_poly"
    "${CMAKE_SOURCE_DIR}/src_polydecomp"
    "${CMAKE_SOURCE_DIR}/src_polygen"
    "${CMAKE_SOURCE_DIR}/src_rankin"
    "${CMAKE_SOURCE_DIR}/src_robust_covering"
    "${CMAKE_SOURCE_DIR}/src_short"
    "${CMAKE_SOURCE_DIR}/src_sparse_solver"
)
target_link_libraries(
  polyhedral_common
  INTERFACE
    Eigen3::Eigen
    Boost::serialization
    Threads::Threads
    PkgConfig::GMPXX
    "${GMP_LIBRARY}"
)

add_library(polyhedral_group_support INTERFACE)
target_link_libraries(
  polyhedral_group_support
  INTERFACE
    polyhedral_common
    polyhedral_nauty
)

add_library(polyhedral_indefinite_support INTERFACE)
target_link_libraries(
  polyhedral_indefinite_support
  INTERFACE
    polyhedral_common
    polyhedral_nauty
)

add_library(polyhedral_isotropy_support INTERFACE)
target_link_libraries(
  polyhedral_isotropy_support
  INTERFACE
    polyhedral_common
    polyhedral_nauty
)

add_library(polyhedral_mpi_support INTERFACE)
target_link_libraries(
  polyhedral_mpi_support
  INTERFACE
    Boost::mpi
    MPI::MPI_CXX
)

add_library(polyhedral_glpk_support INTERFACE)
target_include_directories(polyhedral_glpk_support INTERFACE "${GLPK_INCLUDE_DIR}")
target_link_libraries(polyhedral_glpk_support INTERFACE "${GLPK_LIBRARY}")

add_library(polyhedral_cddgmp_support INTERFACE)
target_include_directories(
  polyhedral_cddgmp_support
  INTERFACE
    ${CDDGMP_INCLUDE_DIRS}
    "${CDD_HEADER_DIR}"
)
target_link_directories(
  polyhedral_cddgmp_support
  INTERFACE
    ${CDDGMP_LIBRARY_DIRS}
)
target_link_libraries(
  polyhedral_cddgmp_support
  INTERFACE
    ${CDDGMP_LINK_LIBRARIES}
)

add_library(polyhedral_cdd_double_support INTERFACE)
target_include_directories(
  polyhedral_cdd_double_support
  INTERFACE
    ${CDDLIB_INCLUDE_DIRS}
    "${CDD_HEADER_DIR}"
)
target_link_directories(
  polyhedral_cdd_double_support
  INTERFACE
    ${CDDLIB_LIBRARY_DIRS}
)
target_link_libraries(
  polyhedral_cdd_double_support
  INTERFACE
    ${CDDLIB_LINK_LIBRARIES}
)


add_library(polyhedral_netcdf_support INTERFACE)
target_include_directories(polyhedral_netcdf_support INTERFACE "${NETCDF_INCLUDE_DIR}")
target_link_libraries(
  polyhedral_netcdf_support
  INTERFACE
    "${NETCDF_CXX_LIBRARY}"
    "${NETCDF_LIBRARY}"
)

add_library(polyhedral_copos_support INTERFACE)
target_link_libraries(
  polyhedral_copos_support
  INTERFACE
    polyhedral_group_support
    polyhedral_glpk_support
)

add_library(polyhedral_dualdesc_support INTERFACE)
target_link_libraries(
  polyhedral_dualdesc_support
  INTERFACE
    polyhedral_group_support
    polyhedral_glpk_support
    polyhedral_cddgmp_support
    polyhedral_mpi_support
)

add_library(polyhedral_perfect_support INTERFACE)
target_link_libraries(
  polyhedral_perfect_support
  INTERFACE
    polyhedral_group_support
    polyhedral_glpk_support
    polyhedral_cddgmp_support
    polyhedral_mpi_support
)

add_library(polyhedral_short_support INTERFACE)
target_link_libraries(
  polyhedral_short_support
  INTERFACE
    polyhedral_group_support
    polyhedral_glpk_support
)

add_library(polyhedral_latt_support INTERFACE)
target_link_libraries(
  polyhedral_latt_support
  INTERFACE
    polyhedral_group_support
    polyhedral_glpk_support
    polyhedral_mpi_support
)

add_library(polyhedral_delaunay_support INTERFACE)
target_link_libraries(
  polyhedral_delaunay_support
  INTERFACE
    polyhedral_group_support
    polyhedral_glpk_support
    polyhedral_mpi_support
)

add_library(polyhedral_lorentzian_support INTERFACE)
target_link_libraries(
  polyhedral_lorentzian_support
  INTERFACE
    polyhedral_common
    polyhedral_nauty
    polyhedral_mpi_support
)

add_library(polyhedral_poly_support INTERFACE)
target_link_libraries(
  polyhedral_poly_support
  INTERFACE
    polyhedral_group_support
    polyhedral_glpk_support
    polyhedral_cddgmp_support
)

add_library(polyhedral_ctype_support INTERFACE)
target_link_libraries(
  polyhedral_ctype_support
  INTERFACE
    polyhedral_group_support
    polyhedral_cdd_double_support
    polyhedral_mpi_support
    polyhedral_netcdf_support
)

add_library(polyhedral_sparse_solver_support INTERFACE)
target_link_libraries(polyhedral_sparse_solver_support INTERFACE polyhedral_common)

add_library(polyhedral_polydecomp_support INTERFACE)
target_link_libraries(
  polyhedral_polydecomp_support
  INTERFACE
    polyhedral_group_support
    polyhedral_glpk_support
    polyhedral_cddgmp_support
)

add_library(polyhedral_polygen_support INTERFACE)
target_link_libraries(
  polyhedral_polygen_support
  INTERFACE
    polyhedral_group_support
    polyhedral_glpk_support
    polyhedral_mpi_support
)

add_library(polyhedral_rankin_support INTERFACE)
target_link_libraries(
  polyhedral_rankin_support
  INTERFACE
    polyhedral_common
    polyhedral_nauty
)

add_library(polyhedral_poincare_support INTERFACE)
target_link_libraries(
  polyhedral_poincare_support
  INTERFACE
    polyhedral_common
    polyhedral_nauty
)

add_library(polyhedral_robust_support INTERFACE)
target_link_libraries(
  polyhedral_robust_support
  INTERFACE
    polyhedral_group_support
    polyhedral_glpk_support
    polyhedral_mpi_support
)

function(add_polyhedral_binary target source support_target)
  add_executable(${target} ${source} ${ARGN})
  target_link_libraries(${target} PRIVATE ${support_target})
endfunction()

set(POLY_C_CDDLIB_MPQ_SOURCE "${CMAKE_SOURCE_DIR}/src_poly/POLY_c_cddlib_mpq.cpp")

set(
  COPOS_PROGRAMS
  CP_CopositiveMin
  CP_CopositiveListCone
  CP_TestCompletePositivity
  CP_CopositiveKernelMin
  CP_CopositiveMaxNorm
  CP_TestCopositivity
  CP_TestStrictCopositivity
)
foreach(target IN LISTS COPOS_PROGRAMS)
  add_polyhedral_binary(
    ${target}
    "${CMAKE_SOURCE_DIR}/src_copos/${target}.cpp"
    polyhedral_copos_support
  )
endforeach()

set(
  DUALDESC_PROGRAMS
  POLY_CheckCanonicalGAP
  POLY_SerialDualDesc
  POLY_MPI_DualDesc
  POLY_RunTheBank
  POLY_DatabaseRestructuration
  POLY_FindNeededBalinski
  POLY_EvaluateDualDesc
  POLY_ReadPartialEnum
  Read_MPI_sizes
  POLY_UpgradeDatabase
  POLY_EvaluateBalinski
  POLY_sampling_facets
  POLY_DirectSerialDualDesc
)
foreach(target IN LISTS DUALDESC_PROGRAMS)
  add_polyhedral_binary(
    ${target}
    "${CMAKE_SOURCE_DIR}/src_dualdesc/${target}.cpp"
    polyhedral_dualdesc_support
    "${POLY_C_CDDLIB_MPQ_SOURCE}"
  )
endforeach()

set(
  GROUP_PROGRAMS
  TEST_WeightMatrixLimited
  GRP_OrbitSplitting
  GRP_VectorSplitting
  GRP_LinearSpace_Equivalence
  GRP_LinearSpace_Stabilizer
  GRP_LinearSpace_Stabilizer_RightCoset
  GRP_ListMat_Vdiag_EXT_Automorphism
  GRP_ListMat_Vdiag_EXT_Isomorphism
  GRP_ListMat_Vdiag_EXT_Invariant
  GRP_LinPolytope_Automorphism
  GRP_LinPolytope_Isomorphism
  GRP_LinPolytopeIntegral_Isomorphism
  GRP_LinPolytopeIntegral_Automorphism
  GRP_LinPolytopeIntegral_Automorphism_RightCoset
  GRP_LinPolytope_Canonic
  GRP_LinPolytope_Invariant
  GRP_RuntimeOrbitSplitting
  GRP_OrbitSplittingPerfect
  GRP_MatrixGroupPermSimplification
  GRP_LinPolytope_Automorphism_GramMat
  GRP_LinPolytope_Isomorphism_GramMat
  GRP_IsomorphismReduction
  GRP_GroupAverage_Matrix
  GRP_GroupAverage_Vector
  GRP_LinearSpace_Stabilizer_DoubleCoset
  GRP_LinPolytopeIntegral_Automorphism_DoubleCoset
  GRP_LinPolytopeIntegral_Automorphism_DoubleCosetStabilizer
  GRP_MatrixGroupSimplification
  GRP_MatrixGroupSimplificationOnline
  GRP_MatrixGroupSimplificationOnlineOpt
)
foreach(target IN LISTS GROUP_PROGRAMS)
  add_polyhedral_binary(
    ${target}
    "${CMAKE_SOURCE_DIR}/src_group/${target}.cpp"
    polyhedral_group_support
  )
endforeach()

set(
  INDEFINITE_PROGRAMS
  INDEF_FORM_GetOrbit_IsotropicKplane
  INDEF_ApproximateOrbitRepresentative
  INDEF_FORM_GetOrbitRepresentative
  INDEF_FORM_TestEquivalence
  INDEF_FORM_AutomorphismGroup
  INDEF_FORM_TestEquivalenceVector
  INDEF_FORM_ApproxCanonicalForm
  INDEF_FORM_InvariantIsotropicPlane
  INDEF_FORM_StabilizerIsotropicPlane
)
foreach(target IN LISTS INDEFINITE_PROGRAMS)
  add_polyhedral_binary(
    ${target}
    "${CMAKE_SOURCE_DIR}/src_indefinite/${target}.cpp"
    polyhedral_indefinite_support
  )
endforeach()

set(
  ISOTROPY_PROGRAMS
  LATT_IndefiniteLLL
  LATT_IndefiniteReduction
  LATT_LLLreduceBasis
  LATT_DetMinimization
  LATT_FindIsotropic
  LATT_TestIsotropic
  LATT_FindPositiveVector
  LATT_lll
  CheckPositiveSemiDefinite
  LATT_DiophantApprox
  VectFamily_Reduction
)
foreach(target IN LISTS ISOTROPY_PROGRAMS)
  add_polyhedral_binary(
    ${target}
    "${CMAKE_SOURCE_DIR}/src_isotropy/${target}.cpp"
    polyhedral_isotropy_support
  )
endforeach()

set(
  PERFECT_PROGRAMS
  PERF_MPI_EnumeratePerfectCones
  IsEutactic
  PERF_SerialPerfectComputation
  PERF_IsBoundedFace
)
foreach(target IN LISTS PERFECT_PROGRAMS)
  add_polyhedral_binary(
    ${target}
    "${CMAKE_SOURCE_DIR}/src_perfect/${target}.cpp"
    polyhedral_perfect_support
    "${POLY_C_CDDLIB_MPQ_SOURCE}"
  )
endforeach()

set(
  SHORT_PROGRAMS
  SHORT_EnumerateCyclicCases
  SHORT_GetShortestVector
  SHORT_CheckPrimeRealizability
  SHORT_TestRealizability
  SHORT_ReduceVectorFamilyGAP
  SHORT_AutomorphismGroup
  SHORT_SplitVectorFamily
)
foreach(target IN LISTS SHORT_PROGRAMS)
  add_polyhedral_binary(
    ${target}
    "${CMAKE_SOURCE_DIR}/src_short/${target}.cpp"
    polyhedral_short_support
  )
endforeach()

set(
  LATT_PROGRAMS
  LATT_Canonicalize
  LATT_CanonicalizeMultiple
  LATT_CanonicalizeSymplectic
  LATT_Automorphism
  LATT_Isomorphism
  LATT_near
  LATT_GenerateCharacteristicVectorSet
  TSPACE_FileFormatConversion
  TSPACE_Stabilizer
  TSPACE_Equivalence
  TSPACE_IntegralSaturation
  TEST_EquiStabFamily
  LATT_ComputeShortestOrbits
  GRP_TestFiniteness
  LATT_ResolveModAction
  LATT_ComputeGroupModAction
)
foreach(target IN LISTS LATT_PROGRAMS)
  add_polyhedral_binary(
    ${target}
    "${CMAKE_SOURCE_DIR}/src_latt/${target}.cpp"
    polyhedral_latt_support
  )
endforeach()

set(
  DELAUNAY_PROGRAMS
  LATT_FindOneVertex
  LATT_MPI_Lattice_IsoDelaunayDomain
  LATT_SerialLattice_IsoDelaunayDomain
  LATT_MPI_ComputeDelaunay
  LATT_SerialComputeDelaunay
)
foreach(target IN LISTS DELAUNAY_PROGRAMS)
  add_polyhedral_binary(
    ${target}
    "${CMAKE_SOURCE_DIR}/src_delaunay/${target}.cpp"
    polyhedral_delaunay_support
  )
endforeach()

set(
  LORENTZIAN_PROGRAMS
  LORENTZ_FundDomain_AllcockEdgewalk
  LORENTZ_FundDomain_AllcockEdgewalk_Isomorphism
  LORENTZ_FundDomain_Vinberg
  LORENTZ_TwoDimAnisotropic_Allcock
  COXDYN_FindExtensions
  COXDYN_ComputeSymbol
  LORENTZ_ComputeStabilizer_Vertex
  COXDYN_GetFacetOneDomain
  COXDYN_FindExtensionsCoxMat
  LATT_GetIntegralMatricesPossibleOrders
  LORENTZ_PERF_Automorphism
  LORENTZ_PERF_Isomorphism
  LORENTZ_MPI_PerfectLorentzian
  LORENTZ_ReflectiveEdgewalk
)
foreach(target IN LISTS LORENTZIAN_PROGRAMS)
  add_polyhedral_binary(
    ${target}
    "${CMAKE_SOURCE_DIR}/src_lorentzian/${target}.cpp"
    polyhedral_lorentzian_support
    "${POLY_C_CDDLIB_MPQ_SOURCE}"
  )
endforeach()

set(
  POLY_PROGRAMS
  POLY_CreateAffineBasis
  POLY_cdd_skeletons
  POLY_SolutionMatNonnegativeComplete
  POLY_redundancy
  POLY_redundancyGroup
  POLY_IntegralPoints
  POLY_SkelettonClarkson
  POLY_FaceLatticeGen
  VectFamily_ColumnReduction
  POLY_lrs
  POLY_SmallPolytope
  POLY_LinearDetermineByInequalities
  POLY_IsPointedCone
  POLY_DirectFaceLattice
  POLY_SolutionMatNonnegative
  POLY_cdd_LinearProgramming
  POLY_dual_description
  POLY_dual_description_group
  POLY_GetFullRankFacetSet
  POLY_GeometricallyUniqueInteriorPoint
  POLY_lrs_volume
  POLY_lrs_triangulation
  POLY_TwoLaminations
  POLY_lrs_triang_facets
)
foreach(target IN LISTS POLY_PROGRAMS)
  add_polyhedral_binary(
    ${target}
    "${CMAKE_SOURCE_DIR}/src_poly/${target}.cpp"
    polyhedral_poly_support
    "${POLY_C_CDDLIB_MPQ_SOURCE}"
  )
endforeach()

set(
  CTYPE_PROGRAMS
  CTYP_PrepareInitialFile
  CTYP_MPI_Enumeration_c
  CTYP_ComputeInvariant
  CTYP_PrepareAdjacencyFile
  CTYP_MPI_EnumerationAdjacencies
  CTYP_CheckAdjacency
  CTYP_SearchMatrix
  CTYP_LookForNoFreeVector
  CTYP_MPI_AdjScheme
  CTYP_ComputeHashStat
  CTYP_ComputeInvariant_B
  NC_ComputeAverage
)
foreach(target IN LISTS CTYPE_PROGRAMS)
  add_polyhedral_binary(
    ${target}
    "${CMAKE_SOURCE_DIR}/src_ctype/${target}.cpp"
    polyhedral_ctype_support
    "${POLY_C_CDDLIB_MPQ_SOURCE}"
  )
endforeach()

set(
  SPARSE_SOLVER_PROGRAMS
  StandaloneSparseSolver
  StandaloneSparseSolver_NNZ
)
foreach(target IN LISTS SPARSE_SOLVER_PROGRAMS)
  add_polyhedral_binary(
    ${target}
    "${CMAKE_SOURCE_DIR}/src_sparse_solver/${target}.cpp"
    polyhedral_sparse_solver_support
  )
endforeach()

set(
  POLYDECOMP_PROGRAMS
  DEC_ComputeDecomposition
  DEC_TestUnionCones
  DEC_TestIntersectionMethods
)
foreach(target IN LISTS POLYDECOMP_PROGRAMS)
  add_polyhedral_binary(
    ${target}
    "${CMAKE_SOURCE_DIR}/src_polydecomp/${target}.cpp"
    polyhedral_polydecomp_support
    "${POLY_C_CDDLIB_MPQ_SOURCE}"
  )
endforeach()

add_polyhedral_binary(
  PolyGen_Difference
  "${CMAKE_SOURCE_DIR}/src_polygen/PolyGen_Difference.cpp"
  polyhedral_polygen_support
)

add_polyhedral_binary(
  RANKIN_Compute_k_min
  "${CMAKE_SOURCE_DIR}/src_rankin/RANKIN_Compute_k_min.cpp"
  polyhedral_rankin_support
)

add_polyhedral_binary(
  POINCARE_Initial_Computation
  "${CMAKE_SOURCE_DIR}/src_poincare_polyhedron/POINCARE_Initial_Computation.cpp"
  polyhedral_poincare_support
  "${POLY_C_CDDLIB_MPQ_SOURCE}"
)

set(
  ROBUST_PROGRAMS
  Robust_RandomEstimation
  Robust_RandomVertexEstimation
  Robust_GetPVoronoi
  Robust_PVoronoi_Stabilizer
  Robust_ExactRobustCoveringDensity
)
foreach(target IN LISTS ROBUST_PROGRAMS)
  add_polyhedral_binary(
    ${target}
    "${CMAKE_SOURCE_DIR}/src_robust_covering/${target}.cpp"
    polyhedral_robust_support
  )
endforeach()
