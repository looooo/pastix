/************************************************************/
/**                                                        **/
/**   NAME       : queue.h                                 **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : queue of PASTIX_INT that sorts elements        **/
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
  PASTIX_INT        size;                  /*+ Allocated memory size             +*/ 
  PASTIX_INT        used;                  /*+ Number of element in the queue    +*/
  PASTIX_INT    *   elttab;                /*+ Array of the element              +*/
  double *   keytab;                /*+ Array of keys                     +*/
  PASTIX_INT    *   keytab2;               /*+ Another array of keys             +*/
} Queue;


#define static

int     queueInit       (Queue *, PASTIX_INT size);
void    queueExit       (Queue *);
Queue * queueCopy       (Queue *dst, Queue *src);
void    queueAdd        (Queue *, PASTIX_INT, double);
void    queueAdd2       (Queue *, PASTIX_INT, double, PASTIX_INT);
PASTIX_INT     queueGet        (Queue *);
PASTIX_INT     queueSize       (Queue *);
void    queueClear      (Queue *);
PASTIX_INT     queueRead       (Queue *);
PASTIX_INT     queueGet2       (Queue *, double *, PASTIX_INT *);
int     queuePossess    (Queue *, PASTIX_INT);
void    queuePrint      (Queue *);

static PASTIX_INT compWith2keys(Queue *, PASTIX_INT, PASTIX_INT);
#undef static
#endif
