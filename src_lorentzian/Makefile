PROGRAM		= LORENTZ_FundDomain_AllcockEdgewalk LORENTZ_FundDomain_AllcockEdgewalk_Isomorphism LORENTZ_FundDomain_Vinberg LORENTZ_TwoDimAnisotropic_Allcock COXDYN_FindExtensions COXDYN_ComputeSymbol LORENTZ_ComputeStabilizer_Vertex COXDYN_GetFacetOneDomain COXDYN_FindExtensionsCoxMat LATT_GetIntegralMatricesPossibleOrders LORENTZ_PERF_Automorphism LORENTZ_PERF_Isomorphism LORENTZ_MPI_PerfectLorentzian LORENTZ_ReflectiveEdgewalk
#PROGRAM		= COXDYN_FindExtensions
#PROGRAM		= COXDYN_ComputeSymbol
#PROGRAM		= LORENTZ_TwoDimAnisotropic_Allcock
#PROGRAM		= LORENTZ_FundDomain_Vinberg
#PROGRAM		= LORENTZ_FundDomain_AllcockEdgewalk LORENTZ_FundDomain_AllcockEdgewalk_Isomorphism
#PROGRAM		= LORENTZ_FundDomain_AllcockEdgewalk COXDYN_FindExtensions
#PROGRAM		= LORENTZ_FundDomain_AllcockEdgewalk LORENTZ_FundDomain_Vinberg
#PROGRAM		= LATT_GetIntegralMatricesPossibleOrders
#PROGRAM		= LORENTZ_FundDomain_AllcockEdgewalk LORENTZ_ComputeStabilizer_Vertex
#PROGRAM		= LORENTZ_ComputeStabilizer_Vertex
#PROGRAM		= LORENTZ_FundDomain_AllcockEdgewalk
#PROGRAM		= LORENTZ_ComputeRoots_Vertex
#PROGRAM		= COXDYN_FindExtensionsCoxMat
#PROGRAM		= LORENTZ_SingleEdgewalk
#PROGRAM		= COXDYN_GetFacetOneDomain

#PROGRAM		= LORENTZ_ComputeDomain
#PROGRAM		= Test

SRCFILES	= POLY_c_cddlib_mpq.cpp
OBJECTS		= $(SRCFILES:.cpp=.o) $(PROGRAM:%=%.o)
LINKOBJECTS	= $(SRCFILES:.cpp=.o)

#CC		= icpc -std=c++17 -Wall -g -O2
#CC		= clang++ -Weverything -std=c++17
#CC		= clang++ -Wall -std=c++17 -g -fsanitize=memory
#CC		= clang++ -g -Wall -Wextra -fsanitize=thread -fsanitize=undefined -fsanitize=bool -fsanitize=bounds -fsanitize=return -fsanitize=null -std=c++17
#CC		= clang++ -g -Wall -Wextra -std=c++17
#CC		= clang++ -O2 -Weverything -std=c++17 -DDEBUG
#CC              = /usr/share/clang/scan-build-3.8/libexec/c++-analyzer -std=c++17 -O2
#CC              = /home/mathieu/opt/gcc-6.3.0/bin/g++ -std=c++17 -O3
#CC		= g++ -std=c++17 -Wall -Wextra -g -O2 -ftemplate-backtrace-limit=0 -fsanitize=float-divide-by-zero -fsignaling-nans
#CC		= g++ -std=c++17 -Wall -Wextra -g -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=bounds -fsanitize=bounds-strict -fsanitize=alignment -fsanitize-undefined-trap-on-error -fstack-protector-all  -fstack-protector-strong -fstack-protector-explicit  -fstack-check
#CC		= g++ -std=c++17 -Wall -Wextra -g -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=bounds -fsanitize=bounds-strict -fsanitize=alignment -fsanitize-undefined-trap-on-error -fstack-protector-all  -fstack-protector-strong -fstack-protector-explicit
#CC		= g++ -std=c++17 -Wall -Wextra -g -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=alignment -fsanitize-undefined-trap-on-error -fstack-protector-all  -fsanitize=pointer-compare -fsanitize-address-use-after-scope -fcf-protection=full
#CC		= g++ -std=c++17 -Wall -Wextra -g -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=alignment -fsanitize-undefined-trap-on-error -fstack-protector-all  -fsanitize=pointer-compare -fsanitize-address-use-after-scope -fcf-protection=full -fstack-check -fno-stack-clash-protection -fsanitize=leak -fsanitize-recover -fbounds-check
#CC		= g++ -std=c++17 -Wall -Wextra -O3 -g -fsanitize=address -fsanitize=null
#CC		= g++ -std=c++17 -Wall -Wextra -O3 -g
CC		= g++ -std=c++17 -Wall -Wextra -g -O3
#CC		= g++ -std=c++17 -Wall -Wextra -g -O3 -fsanitize=address -fsanitize=undefined -fno-sanitize-recover=all -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow -fno-sanitize=null -fno-sanitize=alignment
#CC		= g++ -std=c++17 -Wall -Wextra -g -fsanitize=address -fno-omit-frame-pointer
#CC		= g++ -std=c++17 -Wall -Wextra -g -fsanitize=memory
#CC		= g++ -std=c++17 -Wall -Wextra -g -fstack-check -fno-stack-clash-protection
#CC		= g++ -std=c++17 -Wall -Wextra -g -O2 -static-libgcc -static-libstdc++

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_isotropy -I../src_comb -I../src_graph -I../src_group -I../src_latt -I../src_poly -I../src_dualdesc -I../src_short -I../permutalib/src -I.
#CFLAGS 		= -I$(GMP_INCDIR) -I$(LIBBLISS_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH) $(NAUTY_INCLUDE) $(GLPK_INCLUDE) $(CDDLIB_INCLUDE)
CFLAGS 		= -I$(GMP_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH) $(NAUTY_INCLUDE) $(MPI_INCLUDE_CPP)
CHASHMAP        = -I../basic_common_cpp/sparse-map/include/tsl/ -I../basic_common_cpp/robin-map/include/tsl/ -I../basic_common_cpp/hopscotch-map/include/tsl/
#CDEFS           = -DTIMINGS -DDEBUG_SPECIFIC -DDEBUG -DDEBUG_MATRIX_GROUP -DPERMUTALIB_BLOCKING_SANITY_CHECK
#CDEFS           = -DDEBUG_SPECIFIC -DDEBUG -DDEBUG_MATRIX_GROUP -DPERMUTALIB_BLOCKING_SANITY_CHECK
#CDEFS           = -DTIMINGS
#CDEFS           = -DSANITY_CHECK -DPRINT_SYMBOL_INFORMATION -DDEBUG_MATRIX_GROUP -DDEBUG_LORENTZIAN_LINALG
#CDEFS           = -DSANITY_CHECK -DPRINT_SYMBOL_INFORMATION
#CDEFS           = -DSANITY_CHECK -DPRINT_SYMBOL_INFORMATION -DMATRIX_GROUP_DIAGNOSTICS
#CDEFS           = -DSANITY_CHECK -DPRINT_SYMBOL_INFORMATION -DMATRIX_GROUP_DIAGNOSTICS -DDEBUG_MATRIX_GROUP -DDEBUG_LORENTZIAN_LINALG -DTIMINGS -DDEBUG_POSITIVITY -DTRACK_INFOS_LOG
#CDEFS           =
CDEFS           = -DSANITY_CHECK -DDEBUG -DTIMINGS
#CDEFS           = -DSANITY_CHECK -DPRINT_SYMBOL_INFORMATION -DDEBUG_LORENTZIAN_LINALG -DDEBUG_POSITIVITY -DTRACK_INFOS_LOG -DDEBUG_MATRIX_GROUP -DCHECK_LORENTZIAN_STAB_EQUIV -DDEBUG_LORENTZIAN_STAB_EQUIV -DTIMINGS -DDEBUG_EDGEWALK_GENERIC -DPERMUTALIB_BLOCKING_SANITY_CHECK -DDEBUG_GET_CANONICALIZED_RECORD
CTOTAL          = $(CLOCAL) $(CFLAGS) $(CHASHMAP) $(CDEFS)




#LIBS		= -lm $(GMP_CXX_LINK) -pthread $(LIBBLISS_LINK) $(BOOST_LINK) $(NAUTY_LINK) $(GLPK_LINK) $(CDDLIB_GMP_LINK)
LIBS		= -lm $(GMP_CXX_LINK) -pthread $(BOOST_LINK) $(NAUTY_LINK) $(MPI_LINK_CPP)


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
