#PROGRAM	= hot polarization_gaussian test_lambda
#PROGRAM	= polarization_gaussian
#PROGRAM	= hot_param
#PROGRAM	= test_hessian
PROGRAM	= test_hessian_quadratic test_hot test_hot_param polarization_gaussian test_lambda

OBJECTS		= $(PROGRAM:%=%.o)

#CC		= icpc -std=c++11 -Wall -g -O2
#CC		= clang++ -Weverything -std=c++11
CC		= g++ -std=c++17 -Wall -Wextra -g -O2

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_number -I../src -I../src_matrix -I../src_basic -I../src_comb

CFLAGS 		= $(AUTODIFF_INCLUDE) -I$(EIGEN_PATH)

LIBS		= -lm $(GMP_CXX_LINK) -pthread $(LIBBLISS_LINK)

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
