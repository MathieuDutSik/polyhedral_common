PROGRAM		= LATT_IndefiniteLLL LATT_IndefiniteReduction LATT_LLLreduceBasis LATT_DetMinimization LATT_FindIsotropic LATT_TestIsotropic LATT_FindPositiveVector LATT_lll CheckPositiveSemiDefinite LATT_DiophantApprox
#PROGRAM		= TestPermutationSignCanonic
#PROGRAM		= IndefiniteLLL
#PROGRAM		= IndefiniteReduction

OBJECTS		= $(PROGRAM:%=%.o)

#CC		= g++ -std=c++20 -Wall -Wextra -g -O2 -ftemplate-backtrace-limit=0 -fsanitize=float-divide-by-zero -fsignaling-nans
#CC		= g++ -std=c++20 -Wall -Wextra -O3 -g -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=bounds -fsanitize=bounds-strict -fsanitize=alignment
#CC		= g++ -std=c++20 -Wall -Wextra -g -O2 -static-libgcc -static-libstdc++
CC		= g++ -std=c++20 -Wall -Wextra -O3 -g

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_comb -I../src_graph -I.
CFLAGS 		= -I$(GMP_INCDIR) -I$(BOOST_INCDIR) -isystem $(EIGEN_PATH) $(NAUTY_INCLUDE)

#CDEFS           = -DDEBUG -DTIMINGS -DSANITY_CHECK
CDEFS           = -DGMP_POOL $(CHOICE_COMPILATION)

CTOTAL          = $(CLOCAL) $(CFLAGS) $(CDEFS)

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
