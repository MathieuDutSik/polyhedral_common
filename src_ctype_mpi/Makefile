#PROGRAM		= CTYP_MPI_Enumeration_c CTYP_MergeLogs CTYP_MakeInitialFile CTYP_CheckCanonical
#PROGRAM		= CTYP_MPI_Enumeration_c CTYP_MergeLogs  CTYP_MakeInitialFile
#PROGRAM		= CTYP_MPI_Enumeration_c CTYP_MPI_Enumeration
#PROGRAM		= CTYP_MPI_Enumeration_c
#PROGRAM		= CTYP_MPI_Enumeration_debug
#PROGRAM		= CTYP_MPI_Enumeration_c CTYP_PrepareInitialFile CTYP_ComputeMaxCoefficients CTYP_ComputeInvariant
#PROGRAM		= CTYP_ComputeInvariant
#PROGRAM		= CTYP_PrepareAdjacencyFile
PROGRAM		= CTYP_PrepareInitialFile CTYP_MPI_Enumeration_c CTYP_ComputeInvariant CTYP_PrepareAdjacencyFile CTYP_MPI_EnumerationAdjacencies CTYP_CheckAdjacency CTYP_SearchMatrix
#PROGRAM		= CTYP_PrepareInitialFile
#PROGRAM		= CTYP_MPI_EnumerationAdjacencies
#PROGRAM		= CTYP_CheckAdjacency
#PROGRAM		= CTYP_TimingRandomAccess
#PROGRAM		= CTYP_PrepareInitialFile
#PROGRAM		= CTYP_MakeInitialFile
#PROGRAM		= CTYP_CheckCanonical

OBJECTS		= $(PROGRAM:%=%.o)

#CC		= icpc -std=c++17 -Wall -g -O2
#CC		= clang++ -Weverything -std=c++17
#CC		= clang++ -g -Weverything -std=c++17 -fsanitize=thread -fsanitize=undefined -fsanitize=bool -fsanitize=bounds -fsanitize=bounds -fsanitize=return -fsanitize=null
#CC		= clang++ -g -Weverything -std=c++17
#CC		= clang++ -g -Weverything -std=c++17 -fsanitize=bounds -fsanitize-memory-use-after-dtor -fsanitize-memory-track-origins
#CC		= clang++ -O2 -Weverything -std=c++17 -DDEBUG
#CC              = /usr/share/clang/scan-build-3.8/libexec/c++-analyzer -std=c++17 -O2
#CC              = /home/mathieu/opt/gcc-6.3.0/bin/g++ -std=c++17 -O3
#CC		= g++ -std=c++17 -Wall -g -O2 -ftemplate-backtrace-limit=0 -fsanitize=float-divide-by-zero -fsignaling-nans
#CC		= g++ -std=c++17 -Wall -O3 -g
CC		= g++ -std=c++17 -Wall -O3 -ggdb
#CC		= g++ -std=c++17 -Wall -O3 -g -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=bounds -fsanitize=bounds-strict -fsanitize=alignment
#CC		= g++ -std=c++17 -Wall -g -O2 -static-libgcc -static-libstdc++

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../basic_common_cpp   -I../src_number -I../src_matrix -I../src_comb -I../src_graph -I../src_latt -I../src_poly -I../src_perfect -I.
CFLAGS 		= -I$(GMP_INCDIR) -I$(PERMLIB_INCDIR) -I$(LIBBLISS_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH) $(MPI_INCLUDE_CPP) $(CDDLIB_INCLUDE) $(NAUTY_INCLUDE) $(NETCDF_CXX_ALLINC) $(MPI_CXX_INCLUDE)
CHASHMAP        = -I../basic_common_cpp/sparse-map/include/tsl/ -I../basic_common_cpp/robin-map/include/tsl/ -I../basic_common_cpp/hopscotch-map/include/tsl/
CTOTAL          = $(CLOCAL) $(CFLAGS) $(CHASHMAP)

LIBS		= -lm $(GMP_CXX_LINK) -pthread $(LIBBLISS_LINK) $(BOOST_LINK) $(MPI_LINK_CPP) $(CDDLIB_DOUBLE_LINK) $(NAUTY_LINK) $(NETCDF_CXX_LINK) $(MPI_CXX_LINK)


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
