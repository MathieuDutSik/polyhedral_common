PROGRAM		= Robust_RandomEstimation Robust_InitialVoronoiData Robust_ExactRobustCoveringDensity

OBJECTS		= $(PROGRAM:%=%.o)

CC		= g++ -std=c++17 -Wall -Wextra -g -O3
#   APPLE CLANG
#CC              = g++ -std=c++17 -Wall -Wextra -g -O3 -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=bounds -fsanitize=alignment -Wstrict-aliasing
#   GCC
#CC             = g++ -std=c++17 -Wall -Wextra -g -O3 -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=bounds -fsanitize=bounds-strict -fsanitize=alignment -fsanitize=leak -Wstrict-aliasing

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_graph -I../src_comb -I../src_group -I../src_poly -I../src_isotropy -I../permutalib/src -I../src_latt -I../src_dualdesc -I../src_short -I../src_lorentzian -I../src_delaunay -I../src_enum_schemes -I.
CFLAGS 		= -I$(GMP_INCDIR) -I$(LIBBLISS_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH) $(GLPK_INCLUDE) $(NAUTY_INCLUDE) $(shell pkg-config --cflags ompi)
CHASHMAP        = -I../basic_common_cpp/sparse-map/include/tsl/ -I../basic_common_cpp/robin-map/include/tsl/ -I../basic_common_cpp/hopscotch-map/include/tsl/

#CDEFS           = -DDEBUG -DTIMINGS -DSANITY_CHECK
CDEFS           = $(CHOICE_COMPILATION)

CTOTAL          = $(CLOCAL) $(CFLAGS) $(CHASHMAP) $(CDEFS)

LIBS		= -lm $(GMP_CXX_LINK) -pthread $(LIBBLISS_LINK) $(BOOST_LINK) $(GLPK_LINK) $(NAUTY_LINK) $(shell pkg-config --libs ompi) $(MPI_LINK_CPP)

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
