PROGRAM		= LATT_lll VectFamily_Reduction TEST_MinkowskiReduction TEST_ReductionBenchmark TEST_SeysenReduction TEST_DeepLLL TEST_BKZ TEST_SlideReduction LATT_Canonicalize LATT_CanonicalizeMultiple LATT_CanonicalizeSymplectic LATT_Automorphism LATT_Isomorphism LATT_near LATT_GenerateCharacteristicVectorSet TSPACE_FileFormatConversion TSPACE_Stabilizer TSPACE_Equivalence TEST_EquiStabFamily LATT_ComputeShortestOrbits GRP_TestFiniteness LATT_ResolveModAction LATT_ComputeGroupModAction LATT_RankinMinimum LATT_RankinSublattices

OBJECTS		= $(PROGRAM:%=%.o)

#CC		= g++ -std=c++20 -Wall -Wextra -g -O2 -ftemplate-backtrace-limit=0 -fsanitize=float-divide-by-zero -fsignaling-nans -fsanitize=thread -fsanitize=undefined -fsanitize=bool -fsanitize=bounds
#CC		= g++ -std=c++20 -Wall -Wextra -g -O2 -ftemplate-backtrace-limit=0
CC		= g++ -std=c++20 -Wall -Wextra -g -O3
#   APPLE CLANG
#CC              = g++ -std=c++20 -Wall -Wextra -g -O3 -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=bounds -fsanitize=alignment -Wstrict-aliasing
#   GCC
#CC             = g++ -std=c++20 -Wall -Wextra -g -O3 -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=bounds -fsanitize=bounds-strict -fsanitize=alignment -fsanitize=leak -Wstrict-aliasing

#CC		= g++ -std=c++20 -Wall -Wextra -g -O2 -static-libgcc -static-libstdc++

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_graph -I../src_comb -I../src_group -I../src_poly -I../permutalib/src -I../src_lorentzian -I../src_dualdesc -I../src_isotropy -I../src_short -I../src_enum_schemes -I.
CFLAGS 		= -I$(GMP_INCDIR) -I$(LIBBLISS_INCDIR) -I$(BOOST_INCDIR) -isystem $(EIGEN_PATH) $(NAUTY_INCLUDE) $(shell pkg-config --cflags ompi)
CHASHMAP        = -I../basic_common_cpp/sparse-map/include/tsl/ -I../basic_common_cpp/robin-map/include/tsl/ -I../basic_common_cpp/hopscotch-map/include/tsl/

#CDEFS           = -DDEBUG -DSANITY_CHECK -DTIMINGS
CDEFS           = -DGMP_POOL $(CHOICE_COMPILATION)

CTOTAL          = $(CLOCAL) $(CFLAGS) $(CHASHMAP) $(CDEFS)

LIBS		= -lm $(GMP_CXX_LINK) -pthread $(LIBBLISS_LINK) $(BOOST_LINK) $(NAUTY_LINK) $(shell pkg-config --libs ompi) $(MPI_LINK_CPP)

# Optional flint arithmetic: build with `make ENABLE_FLINT_SUPPORT=1`.
# FLINT_INCDIR / FLINT_LINK can be set in the environment; the defaults
# assume flint is installed next to gmp.
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
		$(CC) $(DEBUGFLAG) $@.o $(LIBS) -o $@

clean:;
	rm -f *.o
	rm -f $(PROGRAM)
