#PROGRAM		= CTYP_MPI_Enumeration CTYP_MergeLogs CTYP_MakeInitialFile CTYP_CheckCanonical
#PROGRAM		= CTYP_MPI_Enumeration CTYP_MergeLogs CTYP_MakeInitialFile
#PROGRAM		= CTYP_MergeLogs
#PROGRAM		= CTYP_MPI_Enumeration
PROGRAM		= CTYP_MPI_Enumeration_c
#PROGRAM		= CTYP_debug
#PROGRAM		= CTYP_MPI_Enumeration_debug
#PROGRAM		= CTYP_MakeInitialFile
#PROGRAM		= CTYP_CheckCanonical
#PROGRAM		= CTYP_CheckCanonical CTYP_MakeInitialFile
#PROGRAM		= CTYP_CheckCanonical CTYP_MakeInitialFile CTYP_MergeLogs

OBJECTS		= $(PROGRAM:%=%.o)

#CC		= icpc -std=c++11 -Wall -g -O2
#CC		= clang++ -Weverything -std=c++11
#CC		= clang++ -g -Weverything -std=c++11 -fsanitize=thread -fsanitize=undefined -fsanitize=bool -fsanitize=bounds -fsanitize=bounds -fsanitize=return -fsanitize=null
#CC		= clang++ -g -Weverything -std=c++11
#CC		= clang++ -g -Weverything -std=c++11 -fsanitize=bounds -fsanitize-memory-use-after-dtor -fsanitize-memory-track-origins
#CC		= clang++ -O2 -Weverything -std=c++11 -DDEBUG
#CC              = /usr/share/clang/scan-build-3.8/libexec/c++-analyzer -std=c++11 -O2
#CC              = /home/mathieu/opt/gcc-6.3.0/bin/g++ -std=c++11 -O3
#CC		= g++ -std=c++11 -Wall -g -O2 -ftemplate-backtrace-limit=0 -fsanitize=float-divide-by-zero -fsignaling-nans
#CC		= g++ -std=c++11 -Wall -O3 -g
CC		= g++ -std=c++11 -Wall -O3 -ggdb
#CC		= g++ -std=c++11 -Wall -O3 -g -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=bounds -fsanitize=bounds-strict -fsanitize=alignment
#CC		= g++ -std=c++11 -Wall -g -O2 -static-libgcc -static-libstdc++

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_comb -I../src_graph -I../src_latt -I../src_poly -I../src_perfect -I.


CFLAGS 		= -I$(GMP_INCDIR) -I$(PERMLIB_INCDIR) -I$(LIBBLISS_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH) $(MPI_INCLUDE_CPP) $(CDDLIB_INCLUDE) $(NAUTY_INCLUDE)


LIBS		= -lm $(GMP_CXX_LINK) -pthread $(LIBBLISS_LINK) $(BOOST_LINK) $(MPI_LINK_CPP) $(CDDLIB_LINK) $(NAUTY_LINK)


%.o: %.cpp
	@echo
	@echo compile $<
	$(CC) $(OUTPUTFLAG) $@ $(CLOCAL) $(CFLAGS) $(NOLINKFLAG) $<

all:	$(PROGRAM)

$(PROGRAM): 	$(OBJECTS)
		@echo
		@echo link $@
		$(CC) $(DEBUGFLAG) $@.o $(LIBS) -o $@

clean:;
	rm -f *.o
	rm -f $(PROGRAM)
