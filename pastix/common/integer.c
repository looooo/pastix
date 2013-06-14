/* Copyright 2004,2007-2009 ENSEIRB, INRIA & CNRS
**
** This file is part of the Scotch software package for static mapping,
** graph partitioning and sparse matrix ordering.
**
** This software is governed by the CeCILL-C license under French law
** and abiding by the rules of distribution of free software. You can
** use, modify and/or redistribute the software under the terms of the
** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
** URL: "http://www.cecill.info".
**
** As a counterpart to the access to the source code and rights to copy,
** modify and redistribute granted by the license, users are provided
** only with a limited warranty and the software's author, the holder of
** the economic rights, and the successive licensors have only limited
** liability.
**
** In this respect, the user's attention is drawn to the risks associated
** with loading, using, modifying and/or developing or reproducing the
** software by the user in light of its specific status of free software,
** that may mean that it is complicated to manipulate, and that also
** therefore means that it is reserved for developers and experienced
** professionals having in-depth computer knowledge. Users are therefore
** encouraged to load and test the software's suitability as regards
** their requirements in conditions enabling the security of their
** systems and/or data to be ensured and, more generally, to use and
** operate it in the same conditions as regards security.
**
** The fact that you are presently reading this means that you have had
** knowledge of the CeCILL-C license and that you accept its terms.
*/
/************************************************************/
/**                                                        **/
/**   NAME       : integer.c                               **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module handles the generic integer **/
/**                type.                                   **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 07 sep 1998     **/
/**                                 to     22 sep 1998     **/
/**                # Version 0.1  : from : 07 jan 2002     **/
/**                                 to     17 jan 2003     **/
/**                # Version 1.0  : from : 23 aug 2005     **/
/**                                 to   : 19 dec 2006     **/
/**                # Version 2.0  : from : 26 feb 2008     **/
/**                                 to   : 26 feb 2008     **/
/**                # Version 5.1  : from : 09 nov 2008     **/
/**                                 to   : 21 jan 2009     **/
/**                                                        **/
/************************************************************/
#include "common.h"
#include <ctype.h>
#include <limits.h>
#include <time.h>

/********************************/
/*                              */
/* Basic routines for fast I/O. */
/*                              */
/********************************/

/* Fast read for pastix_int_t values.
** It returns:
** - 1  : on success.
** - 0  : on error.
*/

int
intLoad (FILE * const         stream,               /*+ Stream to read from     +*/
	 pastix_int_t * const valptr)               /*+ Area where to put value +*/
{
    int                 sign;                       /* Sign flag      */
    int                 car;                        /* Character read */
    pastix_int_t        val;                        /* Value          */

    sign = 0;                                       /* Assume positive constant     */
    for ( ; ; ) {                                   /* Consume whitespaces and sign */
	car = getc (stream);
	if (isspace (car))
	    continue;
	if ((car >= '0') && (car <= '9'))
	    break;
	if (car == '-') {
	    sign = 1;
	    car  = getc (stream);
	    break;
	}
	if (car == '+') {
	    car = getc (stream);
	    break;
	}
	return (0);
    }
    if ((car < '0') || (car > '9'))                 /* If first char is non numeric */
	return (0);                                 /* Then it is an error          */
    val = car - '0';                                /* Get first digit              */
    for ( ; ; ) {
	car = getc (stream);
	if ((car < '0') || (car > '9')) {
	    ungetc (car, stream);
	    break;
	}
	val = val * 10 + (car - '0');               /* Accumulate digits */
    }
    *valptr = (sign != 0) ? (- val) : val;          /* Set result */

    return (1);
}

/* Write routine for pastix_int_t values.
** It returns:
** - 1  : on success.
** - 0  : on error.
*/

int
intSave (FILE * const       stream,               /*+ Stream to write to +*/
	 const pastix_int_t val)                  /*+ Value to write     +*/
{
    return ((fprintf (stream, "%ld", (long) val) == EOF) ? 0 : 1);
}

/**********************************/
/*                                */
/* Permutation building routines. */
/*                                */
/**********************************/

/* This routine fills an array with
** consecutive pastix_int_t values, in
** ascending order.
** It returns:
** - VOID  : in all cases.
*/

void
intAscn (pastix_int_t * const permtab,              /*+ Permutation array to build +*/
	 const pastix_int_t   permnbr,              /*+ Number of entries in array +*/
	 const pastix_int_t   baseval)              /*+ Base value                 +*/
{
    pastix_int_t *permtax;
    pastix_int_t  permnum;
    pastix_int_t  permnnd;

    for (permnum = baseval, permnnd = baseval + permnbr, permtax = permtab - baseval;
	 permnum < permnnd; permnum ++)
	permtax[permnum] = permnum;
}

/* This routine computes a random permutation
** of an array of pastix_int_t values.
** It returns:
** - VOID  : in all cases.
*/

void
intPerm (pastix_int_t * const permtab,              /*+ Permutation array to build +*/
	 const pastix_int_t   permnbr)              /*+ Number of entries in array +*/
{
    pastix_int_t *permptr;
    pastix_int_t  permrmn;

    for (permptr = permtab, permrmn = permnbr;      /* Perform random permutation */
	 permrmn > 0; permptr ++, permrmn --) {
	pastix_int_t permnum;
	pastix_int_t permtmp;

	permnum          = intRandVal (permrmn);    /* Select index to swap       */
	permtmp          = permptr[0];              /* Swap it with current index */
	permptr[0]       = permptr[permnum];
	permptr[permnum] = permtmp;
    }
}

/********************/
/*                  */
/* Random routines. */
/*                  */
/********************/

static volatile int intrandflag = 0;              /*+ Flag set if generator already initialized +*/
static unsigned int intrandseed = 1;              /*+ Random seed                               +*/

/* This routine initializes the pseudo-random
** generator if necessary. In order for multi-sequential
** programs to have exactly the same behavior on any
** process, the random seed does not depend on process
** rank. This routine is not really thread-safe, so it
** should not be called concurrently when it has never
** been initialized before.
** It returns:
** - VOID  : in all cases.
*/

void
intRandInit (void)
{
    if (intrandflag == 0) {                         /* If generator not yet initialized */
#if ! ((defined COMMON_DEBUG) || (defined COMMON_RANDOM_FIXED_SEED) || (defined SCOTCH_DETERMINISTIC))
	intrandseed = time (NULL);                    /* Set random seed if needed */
#endif /* ((defined COMMON_DEBUG) || (defined COMMON_RANDOM_FIXED_SEED)) */
#ifdef COMMON_RANDOM_RAND
	srand (intrandseed);
#else /* COMMON_RANDOM_RAND */
	srandom (intrandseed);
#endif /* COMMON_RANDOM_RAND */
	intrandflag = 1;                              /* Generator has been initialized */
    }
}

/* This routine reinitializes the pseudo-random
** generator to its initial value. This routine
** is not thread-safe.
** It returns:
** - VOID  : in all cases.
*/

void
intRandReset (void)
{
    if (intrandflag != 0) {                         /* Keep seed computed during first initialization */
#ifdef COMMON_RANDOM_RAND
	srand (intrandseed);
#else /* COMMON_RANDOM_RAND */
	srandom (intrandseed);
#endif /* COMMON_RANDOM_RAND */
    }
    else
	intRandInit ();
}

/*********************/
/*                   */
/* Sorting routines. */
/*                   */
/*********************/

/* This routine sorts an array of
** pastix_int_t values by ascending order
** on their first value, used as key.
** It returns:
** - VOID  : in all cases.
*/

#define INTSORTNAME                 intSort1asc1
#define INTSORTSIZE                 (sizeof (pastix_int_t))
#define INTSORTSWAP(p,q)            do {			\
	pastix_int_t t;						\
	t = *((pastix_int_t *) (p));				\
	*((pastix_int_t *) (p)) = *((pastix_int_t *) (q));	\
	*((pastix_int_t *) (q)) = t;				\
    } while (0)
#define INTSORTCMP(p,q)             (*((pastix_int_t *) (p)) < *((pastix_int_t *) (q)))
#include "common_sort.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP

/* This routine sorts an array of pairs of
** pastix_int_t values by ascending order on their
** first value, used as key.
** It returns:
** - VOID  : in all cases.
*/

#define INTSORTNAME                 intSort2asc1
#define INTSORTSIZE                 (2 * sizeof (pastix_int_t))
#define INTSORTSWAP(p,q)            do {				\
	pastix_int_t t, u;						\
	t = *((pastix_int_t *) (p));					\
	u = *((pastix_int_t *) (p) + 1);				\
	*((pastix_int_t *) (p)) = *((pastix_int_t *) (q));		\
	*((pastix_int_t *) (p) + 1) = *((pastix_int_t *) (q) + 1);	\
	*((pastix_int_t *) (q)) = t;					\
	*((pastix_int_t *) (q) + 1) = u;				\
    } while (0)
#define INTSORTCMP(p,q)             (*((pastix_int_t *) (p)) < *((pastix_int_t *) (q)))
#include "common_sort.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP

/* This routine sorts an array of 3-tuples of
** pastix_int_t values by ascending order on their
** first value, used as key.
** It returns:
** - VOID  : in all cases.
*/

#define INTSORTNAME                 intSort3asc1
#define INTSORTSIZE                 (3 * sizeof (pastix_int_t))
#define INTSORTSWAP(p,q)            do {				\
	pastix_int_t t, u, v;						\
	t = *((pastix_int_t *) (p));					\
	u = *((pastix_int_t *) (p) + 1);				\
	v = *((pastix_int_t *) (p) + 2);				\
	*((pastix_int_t *) (p)) = *((pastix_int_t *) (q));		\
	*((pastix_int_t *) (p) + 1) = *((pastix_int_t *) (q) + 1);	\
	*((pastix_int_t *) (p) + 2) = *((pastix_int_t *) (q) + 2);	\
	*((pastix_int_t *) (q)) = t;					\
	*((pastix_int_t *) (q) + 1) = u;				\
	*((pastix_int_t *) (q) + 2) = v;				\
    } while (0)
#define INTSORTCMP(p,q)             (*((pastix_int_t *) (p)) < *((pastix_int_t *) (q)))
#include "common_sort.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP

/* This routine sorts an array of pairs of
** pastix_int_t values by ascending order on both
** of their values, used as primary and
** secondary keys.
** It returns:
** - VOID  : in all cases.
*/

#define INTSORTNAME                 intSort2asc2
#define INTSORTSIZE                 (2 * sizeof (pastix_int_t))
#define INTSORTSWAP(p,q)            do {				\
	pastix_int_t t, u;						\
	t = *((pastix_int_t *) (p));					\
	u = *((pastix_int_t *) (p) + 1);				\
	*((pastix_int_t *) (p)) = *((pastix_int_t *) (q));		\
	*((pastix_int_t *) (p) + 1) = *((pastix_int_t *) (q) + 1);	\
	*((pastix_int_t *) (q)) = t;					\
	*((pastix_int_t *) (q) + 1) = u;				\
    } while (0)
#define INTSORTCMP(p,q)             ((*((pastix_int_t *) (p)) < *((pastix_int_t *) (q))) || ((*((pastix_int_t *) (p)) == *((pastix_int_t *) (q))) && (*((pastix_int_t *) (p) + 1) < *((pastix_int_t *) (q) + 1))))
#include "common_sort.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP
