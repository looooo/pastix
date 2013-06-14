/*
 *  File: integer.h
 *
 *  Authors:
 *    Mathieu  Faverge    - faverge@labri.fr
 *    Xavier   LACOSTE    - lacoste@labri.fr
 *    Pierre   RAMET      - ramet@labri.fr
 *
 *  Dates:
 *    Version 0.0 - from 08 may 1998
 *                  to   08 jan 2001
 *    Version 1.0 - from 06 jun 2002
 *                  to   06 jun 2002
 */
#ifndef _INTEGER_H_
#define _INTEGER_H_

#ifndef MIN
#  define MIN(x,y) (((x)<(y))?(x):(y))
#endif

#ifndef MAX
#  define MAX(x,y) (((x)<(y))?(y):(x))
#endif

int          intLoad     (FILE * const, pastix_int_t * const);
int          intSave     (FILE * const, const pastix_int_t);
void         intAscn     (pastix_int_t * restrict const, const pastix_int_t, const pastix_int_t);
void         intPerm     (pastix_int_t * restrict const, const pastix_int_t);
void         intRandInit (void);
pastix_int_t intRandVal  (pastix_int_t);
void         intSort1asc1(void * const, const pastix_int_t);
void         intSort2asc1(void * const, const pastix_int_t);
void         intSort2asc2(void * const, const pastix_int_t);

#endif /* _INTEGER_H_ */
