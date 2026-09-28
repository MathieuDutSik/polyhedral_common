PROGRAM		= FindNegativeVect

OBJECTS		= $(PROGRAM:%=%.o)

CC		= gcc -Wall -O3 -g

NOLINKFLAG 	= -c
CDEPS           = -MMD -MP
OUTPUTFLAG 	= -o


CFLAGS 		= $(GSL_INCLUDE)

CTOTAL		= $(CFLAGS)

LIBS		= -lm $(GSL_C_LINK)


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
