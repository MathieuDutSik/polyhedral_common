PROGRAM		= INDEF_FORM_GetOrbit_IsotropicKplane INDEF_ApproximateOrbitRepresentative INDEF_FORM_GetOrbitRepresentative INDEF_FORM_TestEquivalence INDEF_FORM_AutomorphismGroup INDEF_FORM_TestEquivalenceVector INDEF_FORM_ApproxCanonicalForm INDEF_FORM_InvariantIsotropicPlane INDEF_FORM_StabilizerIsotropicPlane

OBJECTS		= $(PROGRAM:%=%.o)

#CC		= g++ -std=c++20 -Wall -Wextra -O3 -g -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=bounds -fsanitize=alignment
CC		= g++ -std=c++20 -Wall -Wextra -O3 -g

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_comb -I../src_graph -I../src_group -I../src_latt -I../src_poly -I../src_dualdesc -I../src_isotropy -I../src_lorentzian -I../permutalib/src -I../src_enum_schemes -I.
CHASHMAP        = -I../basic_common_cpp/sparse-map/include/tsl/ -I../basic_common_cpp/robin-map/include/tsl/ -I../basic_common_cpp/hopscotch-map/include/tsl/
CFLAGS 		= -I$(GMP_INCDIR) -I$(BOOST_INCDIR) -isystem $(EIGEN_PATH) $(NAUTY_INCLUDE)

#CDEFS           = -DDEBUG -DTIMINGS -DSANITY_CHECK
CDEFS           = -DGMP_POOL $(CHOICE_COMPILATION)

CTOTAL          = $(CLOCAL) $(CFLAGS) $(CHASHMAP) $(CDEFS)

LIBS		= -lm $(GMP_CXX_LINK) -pthread $(BOOST_LINK) $(NAUTY_LINK)

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
