PROGRAM		= POLYNORM_Packing POLYNORM_Covering POLYNORM_TestCovering

OBJECTS		= $(PROGRAM:%=%.o)

CC		= g++ -std=c++20 -Wall -Wextra -g -O3

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_graph -I../src_comb -I../src_group -I../src_poly -I../permutalib/src -I../src_lorentzian -I../src_dualdesc -I../src_isotropy -I../src_short -I../src_enum_schemes -I../src_latt -I.
CFLAGS 		= -I$(GMP_INCDIR) -I$(LIBBLISS_INCDIR) -I$(BOOST_INCDIR) -isystem $(EIGEN_PATH) $(NAUTY_INCLUDE)
CHASHMAP        = -I../basic_common_cpp/sparse-map/include/tsl/ -I../basic_common_cpp/robin-map/include/tsl/ -I../basic_common_cpp/hopscotch-map/include/tsl/

#CDEFS           = -DDEBUG_POLYNORM_COVERING -DSANITY_CHECK_POLYNORM_COVERING -DTIMINGS_POLYNORM_COVERING
CDEFS           = -DGMP_POOL $(CHOICE_COMPILATION)

CTOTAL          = $(CLOCAL) $(CFLAGS) $(CHASHMAP) $(CDEFS)

LIBS		= -lm $(GMP_CXX_LINK) -pthread $(LIBBLISS_LINK) $(BOOST_LINK) $(NAUTY_LINK)

%.o: %.cpp
	@echo
	@echo compile $<
	$(CC) $(OUTPUTFLAG) $@ $(CTOTAL) $(NOLINKFLAG) $<

all:	$(PROGRAM)

$(PROGRAM): 	$(OBJECTS)
		@echo
		@echo link $@
		$(CC) $@.o $(LIBS) -o $@

clean:;
	rm -f *.o
	rm -f $(PROGRAM)

# The programs are thin drivers over the headers, so a header edit must
# force a rebuild.
$(OBJECTS): PolyNorm_Basic.h PolyNorm_Packing.h PolyNorm_Covering.h
