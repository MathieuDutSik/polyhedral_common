PROGRAM		= CODE_Analysis

SRCFILES	= POLY_c_cddlib_mpq.cpp
OBJECTS		= $(SRCFILES:.cpp=.o) $(PROGRAM:%=%.o)
LINKOBJECTS	= $(SRCFILES:.cpp=.o)

CC		= g++ -std=c++17 -Wall -Wextra -g -O3

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_comb -I../src_graph -I../src_group -I../src_poly -I../src_dualdesc -I../permutalib/src -I.
#CFLAGS 		= -I$(GMP_INCDIR) -I$(LIBBLISS_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH) $(NAUTY_INCLUDE) $(GLPK_INCLUDE) $(CDDLIB_INCLUDE) $(MPREAL_INCLUDE)
CFLAGS 		= -I$(GMP_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH) $(NAUTY_INCLUDE)
CHASHMAP        = -I../basic_common_cpp/sparse-map/include/tsl/ -I../basic_common_cpp/robin-map/include/tsl/ -I../basic_common_cpp/hopscotch-map/include/tsl/
#CDEFS           = -DTIMINGS -DDEBUG_SPECIFIC -DDEBUG -DDEBUG_MATRIX_GROUP -DPERMUTALIB_BLOCKING_SANITY_CHECK
#CDEFS           = -DDEBUG_SPECIFIC -DDEBUG -DDEBUG_MATRIX_GROUP -DPERMUTALIB_BLOCKING_SANITY_CHECK
#CDEFS           = -DTIMINGS
#CDEFS           = -DTIMINGS -DDEBUG
CDEFS           =
CTOTAL          = $(CLOCAL) $(CFLAGS) $(CHASHMAP) $(CDEFS)




#LIBS		= -lm $(GMP_CXX_LINK) -pthread $(LIBBLISS_LINK) $(MPFR_LINK) $(BOOST_LINK) $(NAUTY_LINK) $(GLPK_LINK) $(CDDLIB_GMP_LINK)
LIBS		= -lm $(GMP_CXX_LINK) -pthread $(BOOST_LINK) $(NAUTY_LINK)


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
