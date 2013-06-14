/************************************************************/
/**                                                        **/
/**   NAME       : queue.h                                 **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : queue of pastix_int_t that sorts elements        **/
/**                in ascending way according to a         **/
/**                pastix_float_t key                               **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     08 sep 1998     **/
/**                                                        **/
/************************************************************/

#ifndef QUEUE_H
#define QUEUE_H

/*
**  The type and structure definitions.
*/

typedef struct Queue_ {
  pastix_int_t        size;                  /*+ Allocated memory size             +*/ 
  pastix_int_t        used;                  /*+ Number of element in the queue    +*/
  pastix_int_t    *   elttab;                /*+ Array of the element              +*/
  double *   keytab;                /*+ Array of keys                     +*/
  pastix_int_t    *   keytab2;               /*+ Another array of keys             +*/
} Queue;


#define static

int     queueInit       (Queue *, pastix_int_t size);
void    queueExit       (Queue *);
Queue * queueCopy       (Queue *dst, Queue *src);
void    queueAdd        (Queue *, pastix_int_t, double);
void    queueAdd2       (Queue *, pastix_int_t, double, pastix_int_t);
pastix_int_t     queueGet        (Queue *);
pastix_int_t     queueSize       (Queue *);
void    queueClear      (Queue *);
pastix_int_t     queueRead       (Queue *);
pastix_int_t     queueGet2       (Queue *, double *, pastix_int_t *);
int     queuePossess    (Queue *, pastix_int_t);
void    queuePrint      (Queue *);

static pastix_int_t compWith2keys(Queue *, pastix_int_t, pastix_int_t);
#undef static
#endif
