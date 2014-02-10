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

typedef struct pastix_queue_item_s {
    double          key1;               /*+ Key 1 of the element             +*/
    double          key2;               /*+ Key 2 of the element             +*/
    /* void           *eltptr;             /\*+ Pointer to the element           +*\/ */
    pastix_int_t    eltptr;                 /*+ Pointer to the element           +*/
} pastix_queue_item_t;

typedef struct pastix_queue_s {
    pastix_int_t         size;          /*+ Allocated memory size             +*/
    pastix_int_t         used;          /*+ Number of element in the queue    +*/
    pastix_queue_item_t *elttab;        /*+ Array of the element              +*/
} pastix_queue_t;

int          pqueueInit( pastix_queue_t *, pastix_int_t);
void         pqueueExit( pastix_queue_t *);
pastix_int_t pqueueSize( pastix_queue_t *);
void         pqueueClear(pastix_queue_t *);
/*void         pqueueCopy(pastix_queue_t *, pastix_queue_t *);*/
void         pqueuePush2(pastix_queue_t *, pastix_int_t, double, double);
pastix_int_t pqueueRead (pastix_queue_t *);
pastix_int_t pqueuePop2 (pastix_queue_t *, double *, double *);

static inline void
pqueuePush1(pastix_queue_t *q, pastix_int_t elt, double key1) {
    return pqueuePush2( q, elt, key1, 0 );
}

static inline pastix_int_t
pqueuePop(pastix_queue_t *q){
    return pqueuePop2(q, NULL, NULL);
}

static inline pastix_int_t
pqueuePop1(pastix_queue_t *q, double *key1){
    return pqueuePop2(q, key1, NULL);
}

/* int     queuePossess    (Queue *, pastix_int_t); */
/* void    queuePrint      (Queue *); */

/* static pastix_int_t compWith2keys(Queue *, pastix_int_t, pastix_int_t); */

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
