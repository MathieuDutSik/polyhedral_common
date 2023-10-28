PROGRAM		= TEST_PolytopeIntegralPoints \
		  POLY_CreateAffineBasis \
		  POLY_cdd_skeletons POLY_SolutionMatNonnegativeComplete \
		  POLY_redundancy POLY_redundancyGroup \
		  POLY_IntegralPoints POLY_SkelettonClarkson \
		  POLY_FaceLatticeGen VectFamily_ColumnReduction \
		  POLY_lrs POLY_SmallPolytope \
		  POLY_LinearDetermineByInequalities \
		  POLY_IsPointedCone \
		  POLY_SolutionMatNonnegative POLY_cdd_LinearProgramming \
		  POLY_dual_description POLY_dual_description_group \
		  POLY_ComputeEngelSymbol \
		  POLY_GetFullRankFacetSet POLY_sampling_facets \
		  VectFamily_Reduction

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

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_comb -I../src_graph -I../src_group -I../src_latt -I../src_lorentzian -I../permutalib/src -I.
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
