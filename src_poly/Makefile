#PROGRAM		= POLY_ThreadedADM ComputeDelaunay
PROGRAM		= Temp_StandaloneOrbitFacePolytopeLP \
                  TEST_LinPolytope_Isomorphism \
                  TEST_LinPolytope_Automorphism \
                  TEST_LinPolytope_Canonic \
		  POLY_CreateAffineBasis \
		  POLY_cdd_skeletons \
		  POLY_redundancy \
		  POLY_redundancyClarkson \
		  POLY_redundancyClarksonBlocks \
		  POLY_IntegralPoints \
		  TEST_PolytopeIntegralPoints \
		  POLY_lcdd \
		  POLY_lrs \
		  POLY_lrs_set \
		  POLY_lrs_set_reduction \
		  POLY_lp_cdd \
		  POLY_LinearDetermineByInequalities \
		  POLY_LinearProgram \
                  GRP_TestEquivalenceSets \
                  GRP_OrbitSplitting \
                  GRP_VectorSplitting \
                  GRP_ListMat_Subset_EXT_Automorphism \
                  GRP_ListMat_Subset_EXT_Isomorphism \
                  GRP_ListMat_Subset_EXT_Invariant \
                  GRP_LinPolytope_Automorphism \
                  GRP_LinPolytope_Isomorphism \
                  GRP_LinPolytope_Canonic \
                  GRP_OrbitSplittingPerfect \
		  POLY_ComputeEngelSymbol \
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
#CC		= g++ -std=c++17 -Wall -g -O2 -ftemplate-backtrace-limit=0 -fsanitize=float-divide-by-zero -fsignaling-nans
CC		= g++ -std=c++17 -Wall -O3
#CC		= g++ -std=c++17 -Wall -Wextra -O3
#CC		= g++ -std=c++17 -Wall -O3 -Wconversion
#CC		= g++ -std=c++17 -Wall -Wextra -O3
#CC		= g++ -std=c++17 -Wall -g -O2 -static-libgcc -static-libstdc++

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_comb -I../src_graph -I../src_latt -I.
CFLAGS 		= -I$(GMP_INCDIR)  -I$(PERMLIB_INCDIR) $(PERMUTALIB_INCLUDE) -I$(LIBBLISS_INCDIR) -I$(LIBSHORT_INCDIR) -I$(BOOST_INCDIR) -I$(MPFR_INCDIR) -I$(EIGEN_PATH) -I$(FLINT_INCDIR) -I$(TBB_INCDIR) $(GSL_INCLUDE) $(MPREAL_INCLUDE) $(GLPK_INCLUDE) $(NAUTY_INCLUDE) $(CDDLIB_INCLUDE)
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
