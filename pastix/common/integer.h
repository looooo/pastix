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
void         intSort1asc1(void * const, const pastix_int_t);
void         intSort2asc1(void * const, const pastix_int_t);
void         intSort2asc2(void * const, const pastix_int_t);

/*
  Function: qsortIntFloatAsc

  Sort 2 arrays simultaneously, the first array is an
  array of pastix_int_t and used as key for sorting.
  The second array is an array of PASTIX_FLOAT.

  Parameters:
  pbase       - Array of pointers to the first element of each array to sort.
  total_elems - Number of element in each array.

  Returns:
  Nothing

*/
void qsortIntFloatAsc(void ** const pbase,
                      const pastix_int_t     total_elems);

/*
  Function: qsort2IntFloatAsc

  Sort 3 arrays simultaneously, the first array is an
  array of pastix_int_t and used as primary key for sorting.
  The second array is an other array of pastix_int_t used
  as secondary key.
  The third array is an array of PASTIX_FLOAT.

  Parameters:
  pbase       - Array of pointers to the first element of each array to sort.
  total_elems - Number of element in each array.

  Returns:
  Nothing

*/
void qsort2IntFloatAsc(void ** const pbase,
                       const pastix_int_t     total_elems);


/*
  Function: qsort2IntAsc

  Sort 2 arrays simultaneously, the first array is an
  array of pastix_int_t and used as primary key for sorting.
  The second array is an other array of pastix_int_t used
  as secondary key.

  Parameters:
  pbase       - Array of pointers to the first element of each array to sort.
  total_elems - Number of element in each array.

  Returns:
  Nothing

*/
void qsort2IntAsc(void ** const pbase,
                  const pastix_int_t     total_elems);

/*
  Function: qsort2SmallIntAsc

  Sort 2 arrays simultaneously, the first array is an
  array of integers (int) and used as primary key for sorting.
  The second array is an other array of int used
  as secondary key.

  Parameters:
  pbase       - Array of pointers to the first element of each array to sort.
  total_elems - Number of element in each array.

  Returns:
  Nothing

*/
void qsort2SmallIntAsc(void ** const pbase,
                       const pastix_int_t     total_elems);


pastix_int_t
pastix_intset_union(       pastix_int_t  n1,
                     const pastix_int_t *set1,
                           pastix_int_t  n2,
                     const pastix_int_t *set2,
                           pastix_int_t *set );

static inline pastix_int_t pastix_imin( pastix_int_t a, pastix_int_t b) {
    return ( a < b ) ? a : b;
}

static inline pastix_int_t pastix_imax( pastix_int_t a, pastix_int_t b) {
    return ( a > b ) ? a : b;
}

#endif /* _INTEGER_H_ */
