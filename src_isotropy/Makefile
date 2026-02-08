PROGRAM		= LATT_IndefiniteLLL LATT_IndefiniteReduction LATT_LLLreduceBasis LATT_DetMinimization LATT_FindIsotropic LATT_TestIsotropic LATT_FindPositiveVector LATT_lll CheckPositiveSemiDefinite LATT_DiophantApprox VectFamily_Reduction
#PROGRAM		= TestPermutationSignCanonic
#PROGRAM		= IndefiniteLLL
#PROGRAM		= IndefiniteReduction

OBJECTS		= $(PROGRAM:%=%.o)

#CC		= g++ -std=c++17 -Wall -Wextra -g -O2 -ftemplate-backtrace-limit=0 -fsanitize=float-divide-by-zero -fsignaling-nans
#CC		= g++ -std=c++17 -Wall -Wextra -O3 -g -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=bounds -fsanitize=bounds-strict -fsanitize=alignment
CC		= g++ -std=c++17 -Wall -Wextra -O3 -g
#CC		= g++ -std=c++17 -Wall -Wextra -g -O2 -static-libgcc -static-libstdc++

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_comb -I../src_graph -I.
CFLAGS 		= -I$(GMP_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH) $(NAUTY_INCLUDE)

#CDEFS           = -DDEBUG -DTIMINGS -DSANITY_CHECK
CDEFS           = $(CHOICE_COMPILATION)

CTOTAL          = $(CLOCAL) $(CFLAGS) $(CDEFS)

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
