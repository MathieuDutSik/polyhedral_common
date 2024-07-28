#PROGRAM		= INDEF_ApproximateOrbitRepresentative
PROGRAM		= TEST_EquiStabFamily


OBJECTS		= $(PROGRAM:%=%.o)

CC		= g++ -std=c++17 -Wall -Wextra -O3 -g

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_comb -I../src_graph -I../src_group -I../src_latt -I../src_poly -I../src_dualdesc -I../src_isotropy -I../src_lorentzian -I../permutalib/src -I.
CHASHMAP        = -I../basic_common_cpp/sparse-map/include/tsl/ -I../basic_common_cpp/robin-map/include/tsl/ -I../basic_common_cpp/hopscotch-map/include/tsl/

CFLAGS 		= -I$(GMP_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH) $(NAUTY_INCLUDE)
CTOTAL          = $(CLOCAL) $(CFLAGS) $(CHASHMAP)

LIBS		= -lm $(GMP_CXX_LINK) -pthread $(MPFR_LINK) $(BOOST_LINK) $(NAUTY_LINK)


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
