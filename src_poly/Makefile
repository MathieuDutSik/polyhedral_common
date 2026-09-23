PROGRAM		= POLY_CreateAffineBasis \
		  POLY_cdd_skeletons POLY_SolutionMatNonnegativeComplete \
		  POLY_redundancy POLY_redundancyGroup \
		  POLY_IntegralPoints POLY_SkelettonClarkson \
		  VectFamily_ColumnReduction \
		  POLY_lrs POLY_SmallPolytope \
		  POLY_LinearDetermineByInequalities \
		  POLY_IsPointedCone \
		  POLY_SolutionMatNonnegative POLY_cdd_LinearProgramming \
		  POLY_dual_description \
		  POLY_GetFullRankFacetSet  \
		  POLY_GeometricallyUniqueInteriorPoint \
		  POLY_lrs_volume POLY_lrs_triangulation POLY_TwoLaminations POLY_lrs_triang_facets

SRCFILES	=

OBJECTS		= $(SRCFILES:.cpp=.o) $(PROGRAM:%=%.o)

LINKOBJECTS	= $(SRCFILES:.cpp=.o)

#CC		= g++ -std=c++20 -Wall -Wextra -O2 -g -ftemplate-backtrace-limit=0 -fsanitize=float-divide-by-zero -fsignaling-nans
CC		= g++ -std=c++20 -Wall -Wextra -O3
#CC		= g++ -std=c++20 -Wall -Wextra -O3 -Wconversion
#CC		= g++ -std=c++20 -Wall -Wextra -O2 -g -static-libgcc -static-libstdc++

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_comb -I../src_graph -I../src_group -I../src_isotropy -I../permutalib/src
CFLAGS 		= -I$(GMP_INCDIR) -I$(LIBBLISS_INCDIR) -I$(BOOST_INCDIR) -isystem $(EIGEN_PATH) $(NAUTY_INCLUDE)
CHASHMAP        = -I../basic_common_cpp/sparse-map/include/tsl/ -I../basic_common_cpp/robin-map/include/tsl/ -I../basic_common_cpp/hopscotch-map/include/tsl/

#CDEFS           = -DTIMINGS -DDEBUG_SPECIFIC
CDEFS           = -DGMP_POOL $(CHOICE_COMPILATION)

CTOTAL = $(CLOCAL) $(CFLAGS) $(CHASHMAP) $(CDEFS)


LIBS		= -lm $(GMP_CXX_LINK) -pthread $(LIBBLISS_LINK) $(BOOST_LINK) $(NAUTY_LINK)

# Optional flint arithmetic: build with `make ENABLE_FLINT_SUPPORT=1`.
# FLINT_INCDIR / FLINT_LINK can be set in the environment; the defaults
# assume flint is installed next to gmp. This replaces the older
# unconditional FLINT_INCLUDE / FLINT_LINK hook, which was empty unless
# the environment happened to define it.
ifeq ($(ENABLE_FLINT_SUPPORT),1)
FLINT_INCDIR	?= $(GMP_INCDIR)
FLINT_LINK	?= -lflint
CTOTAL		+= -DENABLE_FLINT_SUPPORT -I$(FLINT_INCDIR)
LIBS		+= $(FLINT_LINK)
endif

# Optional boost multiprecision arithmetics (the gmp_boost / multi_boost
# entries, and the mpq_rational / cpp_rational / mpz_int / cpp_int ones where
# the arithmetics are named by type). They duplicate what gmp already does and
# exist as a consistency check, so they are not compiled by default: build with
# `make ENABLE_BOOST_TYPES=1` to get them. Boost.Multiprecision is header only,
# so this only adds a define.
ifeq ($(ENABLE_BOOST_TYPES),1)
CTOTAL		+= -DENABLE_BOOST_TYPES
endif

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
