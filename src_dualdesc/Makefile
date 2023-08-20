#PROGRAM		= POLY_CheckCanonicalGAP POLY_SerialDualDesc POLY_MPI_DualDesc POLY_RunTheBank Test_face Test_bool POLY_DatabaseRestructuration POLY_FindNeededBalinski POLY_EvaluateDualDesc
#PROGRAM		= POLY_CheckCanonicalGAP
PROGRAM		= POLY_SerialDualDesc
#PROGRAM		= POLY_SerialDualDesc POLY_MPI_DualDesc
#PROGRAM		= POLY_MPI_DualDesc
#PROGRAM		= POLY_gather_vectface
#PROGRAM		= POLY_FindNeededBalinski POLY_SerialDualDesc
#PROGRAM		= POLY_FindNeededBalinski
#PROGRAM		= POLY_MPI_DualDesc POLY_SerialDualDesc POLY_RunTheBank
#PROGRAM		= POLY_RunTheBank
#PROGRAM		= Test_face Test_bool
#PROGRAM		= POLY_DatabaseRestructuration
#PROGRAM		= POLY_EvaluateDualDesc

SRCFILES	= POLY_c_cddlib_mpq.cpp
OBJECTS		= $(SRCFILES:.cpp=.o) $(PROGRAM:%=%.o)
LINKOBJECTS	= $(SRCFILES:.cpp=.o)

#CC		= icpc -std=c++17 -Wall -g -O2
#CC		= clang++ -g -Wall -Wextra -std=c++17
#CC		= clang++ -g -Wall -Wextra -std=c++17 -fsanitize=thread -fsanitize=undefined -fsanitize=bool -fsanitize=bounds -fsanitize=return -fsanitize=null
#CC		= clang++ -g -Wall -Wextra -std=c++17
#CC		= clang++ -g -Wall -Wextra -std=c++17 -fsanitize=bounds -fsanitize-memory-use-after-dtor -fsanitize-memory-track-origins
#CC		= clang++ -O2 -Wall -Wextra -std=c++17 -DDEBUG
#CC              = /usr/share/clang/scan-build-3.8/libexec/c++-analyzer -std=c++17 -O2
#CC              = /home/mathieu/opt/gcc-6.3.0/bin/g++ -std=c++17 -O3
#CC		= g++ -std=c++17 -Wall -Wextra -g -O2 -ftemplate-backtrace-limit=0 -fsanitize=float-divide-by-zero -fsignaling-nans
#CC		= g++ -std=c++17 -Wall -Wextra -O3 -g
#CC		= g++ -std=c++17 -Wall -Wextra -O1 -g
#CC		= g++ -std=c++17 -Wall -Wextra -Wconversion -O3 -g
#CC		= g++-11 -std=c++17 -Wall -Wextra -O3 -g
#CC		= g++ -std=c++17 -Wall -Wextra -Wconversion -O3 -g
#CC		= g++ -std=c++17 -Wall -Wextra -O3 -ggdb
#  GCC
#CC		= g++ -std=c++17 -Wall -Wextra -O3 -g -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=bounds -fsanitize=bounds-strict -fsanitize=alignment
#  CLANG
CC		= g++ -std=c++17 -Wall -Wextra -O3 -g -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=bounds -fsanitize=alignment
#CC		= g++ -std=c++17 -Wall -Wextra -O3 -g -ftemplate-backtrace-limit=0
#CC		= g++ -std=c++17 -Wall -Wextra -O3 -g -pg
#CC		= clang++ -Wall -Wextra -std=c++17
#CC		= g++ -std=c++17 -Wall -Wextra -g -O2 -static-libgcc -static-libstdc++

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

#COPTIONS        = -DPERMUTALIB_BLOCKING_SANITY_CHECK
#COPTIONS        = -DTIMINGS

ENV_OPTIONS     = $(DUALDESC_COMPILATION_OPTIONS)
#COPTIONS        = -DTIMINGS -DPRINT_DOUBLE_COSETS_TEST_PROBLEM -DDEBUG
#COPTIONS        = -DDEBUG
COPTIONS        = -DDEBUG -DTIMINGS -DTRACK_DATABASE -DLRS_PRINT_ANALYSIS
CHASHMAP        = -I../basic_common_cpp/sparse-map/include/tsl/ -I../basic_common_cpp/robin-map/include/tsl/ -I../basic_common_cpp/hopscotch-map/include/tsl/
CLOCAL          = -I../src_basic -I../basic_common_cpp   -I../src_number -I../src_matrix -I../src_comb -I../src_graph -I../src_latt -I../src_poly -I../src_perfect -I../permutalib/src -I../src_lorentzian -I.
CFLAGS 		= -I$(GMP_INCDIR) -I$(PERMLIB_INCDIR) -I$(LIBBLISS_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH) $(MPI_INCLUDE_CPP) $(CDDLIB_INCLUDE) $(GLPK_INCLUDE) $(NAUTY_INCLUDE) -DEIGEN_NO_DEBUG

CTOTAL		= $(CHASHMAP) $(CLOCAL) $(CFLAGS) $(COPTIONS) $(ENV_OPTIONS)

LIBS		= -lm $(GMP_CXX_LINK) -pthread $(LIBBLISS_LINK) $(BOOST_LINK) $(MPI_LINK_CPP) $(CDDLIB_GMP_LINK) $(GLPK_LINK) $(NAUTY_LINK)


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
