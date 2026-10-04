PROGRAM		= IGUSA_EnumerateVertices IGUSA_FacetIncidence IGUSA_FacetEnumeration

OBJECTS		= $(PROGRAM:%=%.o)

CC		= g++ -std=c++20 -Wall -Wextra -O3

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_comb -I../src_graph -I../src_group -I../src_latt -I../src_poly -I../permutalib/src -I../src_lorentzian -I../src_short -I../src_isotropy -I../src_dualdesc -I../src_enum_schemes -I../src_polydecomp -I../src_sparse_solver -I../src_perfect -I../src_milp -I.
CFLAGS 		= -I$(GMP_INCDIR) -I$(LIBBLISS_INCDIR) -I$(BOOST_INCDIR) -isystem $(EIGEN_PATH) $(NAUTY_INCLUDE)
CHASHMAP        = -I../basic_common_cpp/sparse-map/include/tsl/ -I../basic_common_cpp/robin-map/include/tsl/ -I../basic_common_cpp/hopscotch-map/include/tsl/

#CDEFS           = -DDEBUG -DSANITY_CHECK -DTIMINGS
CDEFS           = -DGMP_POOL $(CHOICE_COMPILATION)

CTOTAL          = $(CLOCAL) $(CFLAGS) $(CHASHMAP) $(CDEFS)

LIBS		= -lm $(GMP_CXX_LINK) -pthread $(LIBBLISS_LINK) $(BOOST_LINK) $(NAUTY_LINK)

# Optional flint arithmetic: build with `make ENABLE_FLINT_SUPPORT=1`.
# FLINT_INCDIR / FLINT_LINK can be set in the environment; the defaults
# assume flint is installed next to gmp.
ifeq ($(ENABLE_FLINT_SUPPORT),1)
FLINT_INCDIR	?= $(GMP_INCDIR)
FLINT_LINK	?= -lflint
CTOTAL		+= -DENABLE_FLINT_SUPPORT -I$(FLINT_INCDIR)
LIBS		+= $(FLINT_LINK)
endif

# Optional SCIP integer programming: build with `make ENABLE_SCIP_SUPPORT=1`.
# SCIP_INCDIR / SCIP_LINK can be set in the environment; the defaults
# assume SCIP is installed next to gmp. Without it the integer programs are
# solved by the exact branch and bound of src_milp/integer_linear_programming.h
ifeq ($(ENABLE_SCIP_SUPPORT),1)
SCIP_INCDIR	?= $(GMP_INCDIR)
SCIP_LINK	?= -lscip
CTOTAL		+= -DENABLE_SCIP_SUPPORT -I$(SCIP_INCDIR)
LIBS		+= $(SCIP_LINK)
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
