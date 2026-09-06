PROGRAM		= LATT_ZeroOneSolutions MILP_SystemSymmetry MILP_AffineSymmetry MILP_ZeroOneLattice

OBJECTS		= $(PROGRAM:%=%.o)

CC		= g++ -std=c++20 -Wall -Wextra -O3

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CHASHMAP        = -I../basic_common_cpp/sparse-map/include/tsl/ -I../basic_common_cpp/robin-map/include/tsl/ -I../basic_common_cpp/hopscotch-map/include/tsl/
CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_comb -I../src_graph -I../src_group -I../permutalib/src -I../src_poly -I../src_dualdesc -I../src_latt -I../src_isotropy -I.
COPTIONS        = -DGMP_POOL $(CHOICE_COMPILATION)

CFLAGS 		= -I$(GMP_INCDIR) -I$(LIBBLISS_INCDIR) -I$(BOOST_INCDIR) -isystem $(EIGEN_PATH) $(NAUTY_INCLUDE)

CTOTAL          = $(CLOCAL) $(CFLAGS) $(CHASHMAP) $(COPTIONS)
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
