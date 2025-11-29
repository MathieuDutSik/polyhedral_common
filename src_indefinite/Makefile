#PROGRAM		= INDEF_ApproximateOrbitRepresentative
#PROGRAM		= TEST_EquiStabFamily
PROGRAM		= INDEF_FORM_GetOrbit_IsotropicKplane INDEF_ApproximateOrbitRepresentative INDEF_FORM_GetOrbitRepresentative INDEF_FORM_TestEquivalence INDEF_FORM_AutomorphismGroup INDEF_FORM_TestEquivalenceVector INDEF_FORM_ApproxCanonicalForm INDEF_FORM_InvariantIsotropicPlane INDEF_FORM_StabilizerIsotropicPlane

OBJECTS		= $(PROGRAM:%=%.o)

CC		= g++ -std=c++17 -Wall -Wextra -O3 -g
#CC		= g++ -std=c++17 -Wall -Wextra -O3 -g -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=bounds -fsanitize=alignment

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_comb -I../src_graph -I../src_group -I../src_latt -I../src_poly -I../src_dualdesc -I../src_isotropy -I../src_lorentzian -I../permutalib/src -I../src_enum_schemes -I.
CHASHMAP        = -I../basic_common_cpp/sparse-map/include/tsl/ -I../basic_common_cpp/robin-map/include/tsl/ -I../basic_common_cpp/hopscotch-map/include/tsl/
CFLAGS 		= -I$(GMP_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH) $(NAUTY_INCLUDE)

#CDEFS           = -DDEBUG_LATTICE_STAB_EQUI_CAN -DDEBUG_EQUI_STAB_MEMOIZATION
#CDEFS           = -DDEBUG -DTIMINGS -DSANITY_CHECK -DPERMUTALIB_BLOCKING_SANITY_CHECK
#CDEFS           = -DSANITY_CHECK
#CDEFS           = -DSANITY_CHECK -DPERMUTALIB_BLOCKING_SANITY_CHECK
#CDEFS           = -DDEBUG_INDEFINITE_COMBINED_ALGORITHMS -DDEBUG_LORENTZIAN_PERFECT -DDEBUG_APPROXIMATE_MODELS -DTIMINGS
CDEFS           =

CTOTAL          = $(CLOCAL) $(CFLAGS) $(CHASHMAP) $(CDEFS)

LIBS		= -lm $(GMP_CXX_LINK) -pthread $(BOOST_LINK) $(NAUTY_LINK)


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
