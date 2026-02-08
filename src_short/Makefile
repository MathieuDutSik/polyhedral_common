PROGRAM		= SHORT_EnumerateCyclicCases \
                  SHORT_GetShortestVector \
                  SHORT_CheckPrimeRealizability \
                  SHORT_TestRealizability \
                  SHORT_ReduceVectorFamilyGAP SHORT_AutomorphismGroup \
                  SHORT_SplitVectorFamily

OBJECTS		= $(PROGRAM:%=%.o)

#CC		= g++ -std=c++17 -Wall -Wextra -O2 -g -ftemplate-backtrace-limit=0 -fsanitize=float-divide-by-zero -fsignaling-nans
CC		= g++ -std=c++17 -Wall -Wextra -O3
#CC		= g++ -std=c++17 -Wall -Wextra -O2 -g -static-libgcc -static-libstdc++

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_comb -I../src_graph -I../src_group -I../src_latt -I../src_poly -I../permutalib/src -I../src_lorentzian -I../src_isotropy -I.
CFLAGS 		= -I$(GMP_INCDIR) -I$(LIBBLISS_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH) $(FLINT_INCLUDE) $(GLPK_INCLUDE) $(NAUTY_INCLUDE)
CHASHMAP        = -I../basic_common_cpp/sparse-map/include/tsl/ -I../basic_common_cpp/robin-map/include/tsl/ -I../basic_common_cpp/hopscotch-map/include/tsl/

#CDEFS 		= -DDEBUG -DTIMINGS -DSANITY_CHECK
CDEFS 		=

CTOTAL          = $(CLOCAL) $(CFLAGS) $(CHASHMAP) $(CDEFS)

LIBS		= -lm $(GMP_CXX_LINK) -pthread $(LIBBLISS_LINK) $(FLINT_LINK) $(BOOST_LINK) $(GLPK_LINK) $(NAUTY_LINK)

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
