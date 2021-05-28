PROGRAM		= sv_exact LATT_lll CheckPositiveSemiDefinite LATT_canonicalize LATT_canonicalizeMultiple LATT_canonicalizeSymplectic LATT_DiophantApprox
#PROGRAM		= CheckPositiveSemiDefinite
#PROGRAM		= LATT_canonicalizeMultiple
#PROGRAM		= LATT_canonicalizeSymplectic


OBJECTS		= $(PROGRAM:%=%.o)

#CC		= icpc -std=c++17 -Wall -g -O2
#CC		= clang++ -Weverything -std=c++17
#CC		= clang++ -g -Weverything -fsanitize=thread -fsanitize=undefined -fsanitize=bool -fsanitize=bounds -fsanitize=bounds -fsanitize=return -fsanitize=null -std=c++17
#CC		= clang++ -O2 -Weverything -std=c++17 -DDEBUG
#CC              = /usr/share/clang/scan-build-3.8/libexec/c++-analyzer -std=c++17 -O2
#CC              = /home/mathieu/opt/gcc-6.3.0/bin/g++ -std=c++17 -O3
#CC		= g++ -std=c++17 -Wall -g -O2 -ftemplate-backtrace-limit=0 -fsanitize=float-divide-by-zero -fsignaling-nans -fsanitize=thread -fsanitize=undefined -fsanitize=bool -fsanitize=bounds
#CC		= g++ -std=c++17 -Wall -g -O2 -ftemplate-backtrace-limit=0
CC		= g++ -std=c++17 -Wall -O3
#CC		= g++ -std=c++17 -Wall -g -O2 -static-libgcc -static-libstdc++

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_graph -I../src_comb -I../src_poly -I.


CFLAGS 		= -I$(GMP_INCDIR)  -I$(PERMLIB_INCDIR) -I$(LIBBLISS_INCDIR) -I$(LIBSHORT_INCDIR) -I$(BOOST_INCDIR) -I$(MPFR_INCDIR) -I$(EIGEN_PATH) -I$(GSL_INCDIR) $(MPREAL_INCLUDE) $(GLPK_INCLUDE) $(NAUTY_INCLUDE)

LIBS		= -lm $(GMP_CXX_LINK) -pthread $(LIBSHORT_LINK) $(LIBBLISS_LINK) $(MPFR_LINK) $(BOOST_LINK) $(GSL_C_LINK) $(GLPK_LINK) $(NAUTY_LINK)

%.o: %.cpp
	@echo
	@echo compile $<
	$(CC) $(OUTPUTFLAG) $@ $(CLOCAL) $(CFLAGS) $(NOLINKFLAG) $<

all:	$(PROGRAM)

$(PROGRAM): 	$(OBJECTS)
		@echo
		@echo link $@
		$(CC) $(DEBUGFLAG) $@.o $(LIBS) -o $@

clean:;
	rm -f *.o
	rm -f $(PROGRAM)
