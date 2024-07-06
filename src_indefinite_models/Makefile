PROGRAM		= INDEF_ApproximateOrbitRepresentative


OBJECTS		= $(PROGRAM:%=%.o)

CC		= g++ -std=c++17 -Wall -Wextra -O3 -g

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_comb -I../src_graph -I../src_group -I../src_latt -I../src_poly -I../src_dualdesc -I.
CFLAGS 		= -I$(GMP_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH) $(NAUTY_INCLUDE)
CTOTAL          = $(CLOCAL) $(CFLAGS)

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
