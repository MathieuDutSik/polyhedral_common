#PROGRAM		= POLY_ThreadedADM ComputeDelaunay
PROGRAM		= TEST_LinPolytope_Isomorphism TEST_LinPolytope_Automorphism \
                  TEST_PermMatrTest_PolyEquivalence TEST_PermMatrTest_PolyStabilizer \
                  TEST_LinPolytope_Canonic TEST_WeightMatrixLimited \
		  TEST_LinPolytope_WMatCanonic \
		  TEST_PolytopeIntegralPoints TEST_LinPolytopeIntegral_Automorphism \
		  POLY_CreateAffineBasis \
		  POLY_cdd_skeletons POLY_SolutionMatNonnegativeComplete \
		  POLY_redundancy POLY_redundancyClarkson \
		  POLY_redundancyClarksonGroup POLY_redundancy_Equivariant \
		  POLY_IntegralPoints POLY_SkelettonClarkson \
		  POLY_lcdd POLY_FaceLatticeGen \
		  POLY_lrs POLY_lrs_set_reduction \
		  POLY_LinearDetermineByInequalities \
		  POLY_LinearProgram \
                  GRP_TestEquivalenceSets \
                  GRP_OrbitSplitting GRP_VectorSplitting \
                  GRP_LinearSpace_Equivalence GRP_LinearSpace_Stabilizer \
                  GRP_ListMat_Diag_EXT_Automorphism GRP_ListMat_Diag_EXT_Isomorphism GRP_ListMat_Diag_EXT_Invariant \
                  GRP_LinPolytope_Automorphism GRP_LinPolytope_Isomorphism GRP_LinPolytopeIntegral_Isomorphism \
                  GRP_LinPolytope_Canonic GRP_LinPolytope_Invariant \
                  GRP_RuntimeOrbitSplitting GRP_OrbitSplittingPerfect \
		  POLY_SolutionMatNonnegative POLY_cdd_LinearProgramming \
		  POLY_dual_description POLY_dual_description_group \
		  POLY_ComputeEngelSymbol \
		  POLY_GetFullRankFacetSet POLY_sampling_facets \
		  POLY_IsomorphismReduction \
		  GRP_LinPolytope_Automorphism_GramMat GRP_LinPolytope_Isomorphism_GramMat \
		  POLY_FaceLatticeDirect

SRCFILES	= POLY_c_cddlib_mpq.cpp

OBJECTS		= $(SRCFILES:.cpp=.o) $(PROGRAM:%=%.o)

LINKOBJECTS	= $(SRCFILES:.cpp=.o)

#CC		= icpc -std=c++17 -Wall -g -O2
#CC		= clang++ -Weverything -std=c++17
#CC		= clang++ -g -Weverything -fsanitize=thread -fsanitize=undefined -fsanitize=bool -fsanitize=bounds -fsanitize=bounds -fsanitize=return -fsanitize=null -std=c++17
#CC		= clang++ -O2 -Weverything -std=c++17 -DDEBUG
#CC              = /usr/share/clang/scan-build-3.8/libexec/c++-analyzer -std=c++17 -O2
#CC              = /home/mathieu/opt/gcc-6.3.0/bin/g++ -std=c++17 -O3
#CC		= g++ -std=c++17 -Wall -Wextra -O2 -g -ftemplate-backtrace-limit=0 -fsanitize=float-divide-by-zero -fsignaling-nans
CC		= g++ -std=c++17 -Wall -Wextra -O3
#CC		= g++ -std=c++17 -Wall -Wextra -O3 -Wconversion
#CC		= g++ -std=c++17 -Wall -Wextra -O2 -g -static-libgcc -static-libstdc++

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_comb -I../src_graph -I../src_latt -I../src_lorentzian -I../permutalib/src -I.
CFLAGS 		= -I$(GMP_INCDIR) -I$(LIBBLISS_INCDIR) -I$(LIBSHORT_INCDIR) -I$(BOOST_INCDIR) -I$(MPFR_INCDIR) -I$(EIGEN_PATH) $(FLINT_INCLUDE) -I$(TBB_INCDIR) $(GSL_INCLUDE) $(MPREAL_INCLUDE) $(GLPK_INCLUDE) $(NAUTY_INCLUDE) $(CDDLIB_INCLUDE)
CHASHMAP        = -I../basic_common_cpp/sparse-map/include/tsl/ -I../basic_common_cpp/robin-map/include/tsl/ -I../basic_common_cpp/hopscotch-map/include/tsl/
CDEFS           = -DTIMINGS -DDEBUG_SPECIFIC

CTOTAL = $(CLOCAL) $(CFLAGS) $(CHASHMAP) $(CDEFS)


LIBS		= -lm $(GMP_CXX_LINK) -pthread $(LIBSHORT_LINK) $(LIBBLISS_LINK) $(MPFR_LINK) $(FLINT_LINK) $(BOOST_LINK) $(TBB_LINK) $(GSL_C_LINK) $(GLPK_LINK) $(NAUTY_LINK) $(CDDLIB_GMP_LINK)

%.o: %.cpp
	@echo
	@echo compile $<
	$(CC) $(OUTPUTFLAG) $@ $(CTOTAL) $(NOLINKFLAG) $<

all:	$(PROGRAM)

$(PROGRAM): 	$(OBJECTS)
		@echo
		@echo link $@
		$(CC) $(LINKOBJECTS) $(DEBUGFLAG) $@.o $(LIBS) -o $@

clean:;
	rm -f *.o
	rm -f $(PROGRAM)
