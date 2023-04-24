#PROGRAM		= test
PROGRAM		= GRP_ListMat_Subset_EXT_Automorphism GRP_ListMat_Subset_EXT_Isomorphism sv_exact GRP_ListMat_Subset_EXT_Invariant IndefiniteReduction POLY_dual_description_group POLY_sampling_facets POLY_cdd_lp2
#PROGRAM		= GRP_ListMat_Subset_EXT_Automorphism
#PROGRAM		= sv_exact
#PROGRAM		= POLY_cdd_lp2

OBJECTS		= $(PROGRAM:%=%.o)

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_comb -I../src_graph -I../src_latt -I../src_poly -I../src_indefinite -I../src_lorentzian -I../permutalib/src_gap -I../src_short -I../src_dualdesc -I.
CFLAGS 		= -I$(GMP_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH) $(NAUTY_INCLUDE) -DOSCAR_USE_BOOST_GMP_BINDINGS -DDISABLE_MPQ_CLASS
CHASHMAP	= -I../basic_common_cpp/sparse-map/include/tsl/ -I../basic_common_cpp/robin-map/include/tsl/ -I../basic_common_cpp/hopscotch-map/include/tsl/
CTOTAL          = -std=c++17 $(CLOCAL) $(CFLAGS) $(CHASHMAP)

LIBS		= -lm $(GMP_C_LINK) -pthread $(BOOST_LINK) $(NAUTY_LINK)

%.o: %.cpp
	@echo
	@echo compile $<
	$(CXX) $(OUTPUTFLAG) $@ $(CTOTAL) $(NOLINKFLAG) $<

all:	$(PROGRAM)

$(PROGRAM): 	$(OBJECTS)
		@echo
		@echo link $@
		$(CXX) $(DEBUGFLAG) $@.o $(LIBS) -o $@

clean:;
	rm -f *.o
	rm -f $(PROGRAM)
