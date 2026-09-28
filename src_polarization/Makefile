#PROGRAM	= hot polarization_gaussian test_lambda
#PROGRAM	= polarization_gaussian
#PROGRAM	= hot_param
#PROGRAM	= test_hessian
PROGRAM	= test_hessian_quadratic test_hot test_hot_param polarization_gaussian test_lambda

OBJECTS		= $(PROGRAM:%=%.o)

#CC		= icpc -std=c++11 -Wall -g -O2
#CC		= clang++ -Weverything -std=c++11
CC		= g++ -std=c++20 -Wall -Wextra -g -O2

NOLINKFLAG 	= -c
CDEPS           = -MMD -MP
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_number -I../src -I../src_matrix -I../src_basic -I../src_comb

CFLAGS 		= $(AUTODIFF_INCLUDE) -isystem $(EIGEN_PATH) -DGMP_POOL

LIBS		= -lm $(GMP_CXX_LINK) -pthread $(LIBBLISS_LINK)

%.o: %.cpp
	@echo
	@echo compile $<
	$(CC) $(OUTPUTFLAG) $@ $(CLOCAL) $(CFLAGS) $(CDEPS) $(NOLINKFLAG) $<

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
