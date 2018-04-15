#ifndef PANIC_H
#define PANIC_H

#ifdef __cplusplus
extern "C" {
#endif

/* panic.h  definition of PANIC                               */
/* Version July 11, 2005                                      */
/* Copyright: Frank Vallentin 2005, frank.vallentin@gmail.com */

#define PANIC (panic(__FILE__,__LINE__))

extern void panic (char *filename,
		   int line);

#ifdef __cplusplus
}
#endif

#endif
