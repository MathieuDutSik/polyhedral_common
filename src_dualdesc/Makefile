PROGRAM		= POLY_CheckCanonicalGAP POLY_SerialDualDesc POLY_MPI_DualDesc POLY_RunTheBank POLY_DatabaseRestructuration POLY_FindNeededBalinski POLY_EvaluateDualDesc POLY_ReadPartialEnum Read_MPI_sizes POLY_UpgradeDatabase POLY_EvaluateBalinski POLY_sampling_facets POLY_DirectSerialDualDesc

SRCFILES	= POLY_c_cddlib_mpq.cpp
OBJECTS		= $(SRCFILES:.cpp=.o) $(PROGRAM:%=%.o)
LINKOBJECTS	= $(SRCFILES:.cpp=.o)

CC		= g++ -std=c++17 -Wall -Wextra -O3 -g
#CC		= g++ -std=c++17 -Wall -Wextra -O1 -g
#  GCC
#CC		= g++ -std=c++17 -Wall -Wextra -O3 -g -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=bounds -fsanitize=bounds-strict -fsanitize=alignment
#  CLANG
#CC		= g++ -std=c++17 -Wall -Wextra -O3 -g -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=bounds -fsanitize=alignment
#CC		= g++ -std=c++17 -Wall -Wextra -O3 -g -ftemplate-backtrace-limit=0
#CC		= g++ -std=c++17 -Wall -Wextra -O3 -g -pg
#CC		= clang++ -Wall -Wextra -std=c++17
#CC		= g++ -std=c++17 -Wall -Wextra -g -O2 -static-libgcc -static-libstdc++

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

#COPTIONS        = -DPERMUTALIB_BLOCKING_SANITY_CHECK
#COPTIONS        = -DTIMINGS

ENV_OPTIONS     = $(DUALDESC_COMPILATION_OPTIONS)

#COPTIONS        = -DDEBUG -DTIMINGS -DTRACK_DATABASE -DLRS_PRINT_ANALYSIS
COPTIONS        =

CHASHMAP        = -I../basic_common_cpp/sparse-map/include/tsl/ -I../basic_common_cpp/robin-map/include/tsl/ -I../basic_common_cpp/hopscotch-map/include/tsl/
CLOCAL          = -I../src_basic -I../basic_common_cpp   -I../src_number -I../src_matrix -I../src_comb -I../src_graph -I../src_group -I../src_poly -I../src_isotropy -I../permutalib/src -I.
CFLAGS 		= -I$(GMP_INCDIR) -I$(LIBBLISS_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH) $(shell pkg-config --cflags ompi) $(CDDLIB_INCLUDE) $(GLPK_INCLUDE) $(NAUTY_INCLUDE) -DEIGEN_NO_DEBUG

CTOTAL		= $(CHASHMAP) $(CLOCAL) $(CFLAGS) $(COPTIONS) $(ENV_OPTIONS)

LIBS		= -lm $(GMP_CXX_LINK) -pthread $(LIBBLISS_LINK) $(BOOST_LINK) $(shell pkg-config --libs ompi) $(MPI_LINK_CPP) $(CDDLIB_GMP_LINK) $(GLPK_LINK) $(NAUTY_LINK)


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
