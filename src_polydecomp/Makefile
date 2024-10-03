PROGRAM		= DEC_ComputeDecomposition DEC_TestUnionCones DEC_TestIntersectionMethods

SRCFILES        = POLY_c_cddlib_mpq.cpp
OBJECTS		= $(SRCFILES:.cpp=.o) $(PROGRAM:%=%.o)
LINKOBJECTS     = $(SRCFILES:.cpp=.o)

#CC		= clang++ -Wall -std=c++17
#CC		= clang++ -g -Weverything -fsanitize=thread -fsanitize=undefined -fsanitize=bool -fsanitize=bounds -fsanitize=bounds -fsanitize=return -fsanitize=null -std=c++17
#CC		= clang++ -O2 -Weverything -std=c++17 -DDEBUG
#CC              = /usr/share/clang/scan-build-3.8/libexec/c++-analyzer -std=c++17 -O2
#CC              = /home/mathieu/opt/gcc-6.3.0/bin/g++ -std=c++17 -O3
#CC		= g++ -std=c++17 -Wall -Wextra -O2 -g -ftemplate-backtrace-limit=0 -fsanitize=float-divide-by-zero -fsignaling-nans
#CC		= g++ -std=c++17 -Wall -Wextra -O3 -g -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=bounds -fsanitize=bounds-strict -fsanitize=alignment
CC		= g++ -std=c++17 -Wall -Wextra -O3 -g
#CC		= g++ -std=c++17 -Wall -Wextra -O2 -g -static-libgcc -static-libstdc++

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_comb -I../src_graph -I../src_group -I../src_latt -I../src_poly -I../src_dualdesc -I../permutalib/src -I../src_lorentzian -I../src_isotropy -I.
CDEFS           = -DPERMUTALIB_BLOCKING_SANITY_CHECK -DDEBUG_POLYEDRAL_DECOMPOSITION
#CDEFS           = -DPERMUTALIB_BLOCKING_SANITY_CHECK -DDEBUG_MATRIX_GROUP
CFLAGS 		= -I$(GMP_INCDIR) -I$(LIBBLISS_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH) $(MPI_INCLUDE_CPP) $(NAUTY_INCLUDE) $(GLPK_INCLUDE) $(CDDLIB_INCLUDE) $(MPI_INCLUDE_CPP)
CHASHMAP        = -I../basic_common_cpp/sparse-map/include/tsl/ -I../basic_common_cpp/robin-map/include/tsl/ -I../basic_common_cpp/hopscotch-map/include/tsl/
CTOTAL          = $(CLOCAL) $(CFLAGS) $(CHASHMAP) $(CDEFS)

LIBS		= -lm $(GMP_CXX_LINK) -pthread $(LIBBLISS_LINK) $(MPFR_LINK) $(NAUTY_LINK) $(GLPK_LINK) $(CDDLIB_GMP_LINK) $(MPI_LINK_CPP)


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
