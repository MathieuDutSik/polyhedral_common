PROGRAM		= PERIODIC_NonLatticePacking

OBJECTS		= $(PROGRAM:%=%.o)

CC		= g++ -std=c++20 -Wall -Wextra -g -O3

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

QHULL_INCDIR ?= /opt/homebrew/include
QHULL_LINK   ?= -L/opt/homebrew/lib -lqhullstatic_r

CFLAGS 		= -I$(QHULL_INCDIR) -isystem $(EIGEN_PATH)

CTOTAL          = $(CFLAGS)

LIBS		= -lm $(QHULL_LINK)

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

# header dependency: the program is a thin driver over these headers, so a
# header edit must force a rebuild (make does not scan includes on its own).
$(OBJECTS): NonLatticePacking.h ../src_delaunay/JointCoveringDouble.h
