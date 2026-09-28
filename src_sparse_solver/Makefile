PROGRAM		= StandaloneSparseSolver StandaloneSparseSolver_NNZ

OBJECTS		= $(PROGRAM:%=%.o)

#CC		= icpc -std=c++20 -Wall -g -O2
#CC		= clang++ -Weverything -std=c++20
#CC		= clang++ -Weverything -fsanitize=thread -fsanitize=undefined -fsanitize=bool -fsanitize=bounds -fsanitize=bounds -fsanitize=return -fsanitize=null -std=c++20
#CC		= g++ -std=c++20 -Wall -Wextra -O2 -g -ftemplate-backtrace-limit=0 -fsanitize=float-divide-by-zero -fsignaling-nans
CC		= g++ -std=c++20 -Wall -Wextra -O3
#CC		= g++ -std=c++20 -Wall -Wextra -O2 -g -static-libgcc -static-libstdc++

NOLINKFLAG 	= -c
CDEPS           = -MMD -MP
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix
CFLAGS 		= -I$(BOOST_INCDIR) -isystem $(EIGEN_PATH) -I$(GMP_INCDIR)

#CDEFS           = -DDEBUG -DTIMINGS -DSANITY_CHECK
CDEFS           = -DGMP_POOL $(CHOICE_COMPILATION)

CTOTAL          = $(CLOCAL) $(CFLAGS) $(CDEFS)

LIBS		= -lm -pthread $(GMP_CXX_LINK)

%.o: %.cpp
	@echo
	@echo compile $<
	$(CC) $(OUTPUTFLAG) $@ $(CTOTAL) $(CDEPS) $(NOLINKFLAG) $<

all:	$(PROGRAM)

$(PROGRAM): 	$(OBJECTS)
		@echo
		@echo link $@
		$(CC) $(DEBUGFLAG) $@.o $(LIBS) -o $@

clean:;
	rm -f *.o
	rm -f *.d
	rm -f $(PROGRAM)

# The header dependencies, written by the compiler as it compiles (-MMD),
# so that a change to a header rebuilds the objects that include it.
-include $(OBJECTS:.o=.d)
