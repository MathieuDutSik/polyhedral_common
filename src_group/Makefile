PROGRAM		= TEST_LinPolytope_Isomorphism TEST_LinPolytope_Automorphism \
                  TEST_LinPolytope_Canonic TEST_WeightMatrixLimited \
		  TEST_LinPolytope_WMatCanonic \
		  TEST_LinPolytopeIntegral_Automorphism \
                  GRP_OrbitSplitting GRP_VectorSplitting \
                  GRP_LinearSpace_Equivalence GRP_LinearSpace_Stabilizer \
                  GRP_LinearSpace_Stabilizer_RightCoset \
                  GRP_ListMat_Vdiag_EXT_Automorphism GRP_ListMat_Vdiag_EXT_Isomorphism GRP_ListMat_Vdiag_EXT_Invariant \
                  GRP_LinPolytope_Automorphism GRP_LinPolytope_Isomorphism \
                  GRP_LinPolytopeIntegral_Isomorphism GRP_LinPolytopeIntegral_Automorphism \
                  GRP_LinPolytopeIntegral_Automorphism_RightCoset \
                  GRP_LinPolytope_Canonic GRP_LinPolytope_Invariant \
                  GRP_RuntimeOrbitSplitting GRP_OrbitSplittingPerfect GRP_MatrixGroupSimplification \
		  GRP_LinPolytope_Automorphism_GramMat GRP_LinPolytope_Isomorphism_GramMat \
		  GRP_IsomorphismReduction GRP_GroupAverage_Matrix GRP_GroupAverage_Vector \
		  GRP_LinearSpace_Stabilizer_DoubleCoset GRP_LinPolytopeIntegral_Automorphism_DoubleCoset \
		  GRP_LinPolytopeIntegral_Automorphism_DoubleCosetStabilizer \
		  GRP_MatrixGroupPermSimplification GRP_MatrixGroupPermSimplificationOnline \
		  GRP_MatrixGroupPermSimplificationOnlineOpt

OBJECTS		= $(PROGRAM:%=%.o)

#CC		= icpc -std=c++17 -Wall -g -O2
#CC		= clang++ -Weverything -std=c++17
#CC		= clang++ -g -Weverything -fsanitize=thread -fsanitize=undefined -fsanitize=bool -fsanitize=bounds -fsanitize=bounds -fsanitize=return -fsanitize=null -std=c++17
#CC		= clang++ -O2 -Weverything -std=c++17 -DDEBUG
#CC              = /usr/share/clang/scan-build-3.8/libexec/c++-analyzer -std=c++17 -O2
#CC              = /home/mathieu/opt/gcc-6.3.0/bin/g++ -std=c++17 -O3
#CC		= g++ -std=c++17 -Wall -Wextra -O2 -g -ftemplate-backtrace-limit=0 -fsanitize=float-divide-by-zero -fsignaling-nans
CC		= g++ -std=c++17 -Wall -Wextra -O3
#CC		= g++ -std=c++17 -Wall -Wextra -O3 -Wconversion
#CC		= g++ -std=c++17 -Wall -Wextra -O2 -g -static-libgcc -static-libstdc++

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_comb -I../src_graph -I../src_latt -I../src_lorentzian -I../src_isotropy -I../permutalib/src -I../src_poly -I.
CFLAGS 		= -I$(GMP_INCDIR) -I$(LIBBLISS_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH) $(NAUTY_INCLUDE)
CHASHMAP        = -I../basic_common_cpp/sparse-map/include/tsl/ -I../basic_common_cpp/robin-map/include/tsl/ -I../basic_common_cpp/hopscotch-map/include/tsl/
#CDEFS           = -DTIMINGS
#CDEFS           = -DTIMINGS -DDEBUG
#CDEFS           = -DTIMINGS -DUSE_BLISS
CDEFS           = -DUSE_BLISS -DDEBUG_DOUBLE_COSET_ENUM

CTOTAL = $(CLOCAL) $(CFLAGS) $(CHASHMAP) $(CDEFS)


LIBS		= -lm $(GMP_CXX_LINK) -pthread $(LIBBLISS_LINK) $(BOOST_LINK) $(NAUTY_LINK)

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
