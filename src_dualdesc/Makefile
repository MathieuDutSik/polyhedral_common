#PROGRAM		= DD_ConvertGAP_to_NC POLY_CheckCanonicalGAP POLY_NetcdfToGAP POLY_SerialDualDesc ConvertNC_to_NC_FB_FF Test_face Test_bool
#PROGRAM		= DD_ConvertGAP_to_NC
#PROGRAM		= POLY_CheckCanonicalGAP
#PROGRAM		= POLY_NetcdfToGAP
PROGRAM		= POLY_SerialDualDesc
#PROGRAM		= POLY_RunTheBank
#PROGRAM		= ConvertNC_to_NC_FB_FF
#PROGRAM		= Test_face Test_bool

SRCFILES	= POLY_c_cddlib_mpq.cpp
OBJECTS		= $(SRCFILES:.cpp=.o) $(PROGRAM:%=%.o)
LINKOBJECTS	= $(SRCFILES:.cpp=.o)

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
#CC		= g++ -std=c++17 -Wall -Wconversion -O3 -g
#CC		= g++-11 -std=c++17 -Wall -O3 -g
#CC		= g++ -std=c++17 -Wall -Wconversion -O3 -g
#CC		= g++ -std=c++17 -Wall -O3 -ggdb
#CC		= g++ -std=c++17 -Wall -O3 -g -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=bounds -fsanitize=bounds-strict -fsanitize=alignment
CC		= g++ -std=c++17 -Wall -O3 -g
#CC		= g++ -pg -std=c++17 -Wall -O3 -g
#CC		= clang++ -Weverything -std=c++17
#CC		= g++ -std=c++17 -Wall -g -O2 -static-libgcc -static-libstdc++

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

#COPTIONS        = -DPERMUTALIB_BLOCKING_SANITY_CHECK
COPTIONS        = -DTIMINGS
CHASHMAP        = -I../basic_common_cpp/sparse-map/include/tsl/ -I../basic_common_cpp/robin-map/include/tsl/ -I../basic_common_cpp/hopscotch-map/include/tsl/
CLOCAL          = -I../src_basic -I../basic_common_cpp   -I../src_number -I../src_matrix -I../src_comb -I../src_graph -I../src_latt -I../src_poly -I../src_perfect -I.
CFLAGS 		= -I$(GMP_INCDIR) -I$(PERMLIB_INCDIR) $(PERMUTALIB_INCLUDE) -I$(LIBBLISS_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH) $(MPI_INCLUDE_CPP) $(CDDLIB_INCLUDE) $(GLPK_INCLUDE) $(NAUTY_INCLUDE) $(NETCDF_CXX_ALLINC) $(MPI_CXX_INCLUDE) -DEIGEN_NO_DEBUG

CTOTAL		= $(CHASHMAP) $(CLOCAL) $(CFLAGS) $(COPTIONS)

LIBS		= -lm $(GMP_CXX_LINK) -pthread $(LIBBLISS_LINK) $(BOOST_LINK) $(MPI_LINK_CPP) $(CDDLIB_GMP_LINK) $(GLPK_LINK) $(NAUTY_LINK) $(NETCDF_CXX_LINK) $(MPI_CXX_LINK)


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
