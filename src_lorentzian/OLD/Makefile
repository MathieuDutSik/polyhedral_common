PROGRAM		= FindNegativeVect

OBJECTS		= $(PROGRAM:%=%.o)

CC		= gcc -Wall -O3 -g

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o


CFLAGS 		= $(GSL_INCLUDE)

CTOTAL		= $(CFLAGS)

LIBS		= -lm $(GSL_C_LINK)


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
