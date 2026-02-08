PROGRAM		= RANKIN_Compute_k_min
#PROGRAM		= LORENTZ_ComputeDomain
#PROGRAM		= Test

OBJECTS		= $(PROGRAM:%=%.o)

#CC		= g++ -std=c++17 -Wall -Wextra -g -fsanitize=address
CC		= g++ -std=c++17 -Wall -Wextra -g -O2
#CC		= g++ -std=c++17 -Wall -Wextra -g -O3
#CC		= g++ -std=c++17 -Wall -Wextra -g -O3 -fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow -fno-sanitize=null -fno-sanitize=alignment

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_isotropy -I../src_comb -I../src_graph -I../src_group -I../src_latt -I../src_poly -I../src_dualdesc -I../src_short -I../permutalib/src -I.
#CFLAGS 		= -I$(GMP_INCDIR) -I$(LIBBLISS_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH) $(NAUTY_INCLUDE) $(GLPK_INCLUDE) $(CDDLIB_INCLUDE)
CFLAGS 		= -I$(GMP_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH) $(NAUTY_INCLUDE)
CHASHMAP        = -I../basic_common_cpp/sparse-map/include/tsl/ -I../basic_common_cpp/robin-map/include/tsl/ -I../basic_common_cpp/hopscotch-map/include/tsl/

#CDEFS           = -DSANITY_CHECK -DDEBUG -DTIMINGS
CDEFS           =

CTOTAL          = $(CLOCAL) $(CFLAGS) $(CHASHMAP) $(CDEFS)

#LIBS		= -lm $(GMP_CXX_LINK) -pthread $(LIBBLISS_LINK) $(BOOST_LINK) $(NAUTY_LINK) $(GLPK_LINK) $(CDDLIB_GMP_LINK)
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
