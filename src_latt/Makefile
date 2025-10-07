PROGRAM		= LATT_Canonicalize LATT_CanonicalizeMultiple LATT_CanonicalizeSymplectic LATT_Automorphism LATT_Isomorphism LATT_near LATT_MPI_ComputeDelaunay LATT_MPI_Lattice_IsoDelaunayDomain LATT_GenerateCharacteristicVectorSet TSPACE_FileFormatConversion TSPACE_Stabilizer TSPACE_Equivalence TEST_EquiStabFamily LATT_ComputeShortestOrbits LATT_FindOneVertex LATT_SerialLattice_IsoDelaunayDomain LATT_SerialComputeDelaunay

OBJECTS		= $(PROGRAM:%=%.o)

#CC		= icpc -std=c++17 -Wall -g -O2
#CC		= clang++ -Weverything -std=c++17
#CC		= clang++ -g -Weverything -fsanitize=thread -fsanitize=undefined -fsanitize=bool -fsanitize=bounds -fsanitize=bounds -fsanitize=return -fsanitize=null -std=c++17
#CC		= clang++ -O2 -Weverything -std=c++17 -DDEBUG
#CC              = /usr/share/clang/scan-build-3.8/libexec/c++-analyzer -std=c++17 -O2
#CC              = /home/mathieu/opt/gcc-6.3.0/bin/g++ -std=c++17 -O3
#CC		= g++ -std=c++17 -Wall -Wextra -g -O2 -ftemplate-backtrace-limit=0 -fsanitize=float-divide-by-zero -fsignaling-nans -fsanitize=thread -fsanitize=undefined -fsanitize=bool -fsanitize=bounds
#CC		= g++ -std=c++17 -Wall -Wextra -g -O2 -ftemplate-backtrace-limit=0
CC		= g++ -std=c++17 -Wall -Wextra -g -O3
#   APPLE CLANG
#CC              = g++ -std=c++17 -Wall -Wextra -g -O3 -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=bounds -fsanitize=alignment -Wstrict-aliasing
#   GCC
#CC             = g++ -std=c++17 -Wall -Wextra -g -O3 -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=bounds -fsanitize=bounds-strict -fsanitize=alignment -fsanitize=leak -Wstrict-aliasing

#CC		= g++ -std=c++17 -Wall -Wextra -g -O2 -static-libgcc -static-libstdc++

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_graph -I../src_comb -I../src_group -I../src_poly -I../permutalib/src -I../src_lorentzian -I../src_dualdesc -I../src_isotropy -I../src_short -I../src_enum_schemes -I.
CFLAGS 		= -I$(GMP_INCDIR) -I$(LIBBLISS_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH) $(GLPK_INCLUDE) $(NAUTY_INCLUDE) $(shell pkg-config --cflags ompi)
CHASHMAP        = -I../basic_common_cpp/sparse-map/include/tsl/ -I../basic_common_cpp/robin-map/include/tsl/ -I../basic_common_cpp/hopscotch-map/include/tsl/
#CDEFS           = -DDEBUG_CANONIC -DTIMINGS
#CDEFS           = -DCHECK_SHVEC
CDEFS           =
#CDEFS           = -DDEBUG -DSANITY_CHECK -DTIMINGS
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
