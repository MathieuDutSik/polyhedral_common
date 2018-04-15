/* panic.c  definition of PANIC                               */
/* Version July 11, 2005                                      */
/* Copyright: Frank Vallentin 2005, frank.vallentin@gmail.com */


#include "panic.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

void panic (char *filename,
	    int line)
{
  fprintf(stderr, "\nPanic in line %d of file %s\n", line, filename);
  perror("Unexpected library error");
  abort();
}
