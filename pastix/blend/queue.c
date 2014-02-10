/*
  File: Queue.c

  Operations on the queue structure.

*/
#include <stdio.h>

#include "common.h"
#include "queue.h"

static inline int
pqueueItemComparison(pastix_queue_t *q,
                     pastix_int_t    elt1,
                     pastix_int_t    elt2)
{
    /* if elt1 < elt2 return 1  */
    /* if elt1 = elt2 return 0  */
    /* if elt1 > elt2 return 0 */

    pastix_queue_item_t *item1 = q->elttab + elt1;
    pastix_queue_item_t *item2 = q->elttab + elt2;

    if ( item1->key1 == item2->key1)
        return item1->key2 < item2->key2;
    else
        return item1->key1 < item2->key1;
}

int
pqueueInit(pastix_queue_t *q,
           pastix_int_t    size)
{
    q->size = size;
    q->used = 0;
    if (q->size != 0)
    {
        MALLOC_INTERN(q->elttab, size, pastix_queue_item_t);
    }
    else
    {
        q->elttab  = NULL;
    }
    return PASTIX_SUCCESS;
}

void
pqueueExit(pastix_queue_t *q)
{
    if(q->size != 0)
    {
        memFree_null(q->elttab);
    }
    q->size = 0;
}

pastix_int_t
pqueueSize(pastix_queue_t *q)
{
  return q->used;
}

void
pqueueClear(pastix_queue_t *q)
{
  q->used = 0;
}

void
pqueuePush2(pastix_queue_t *q,
            void *elt,
            double key1,
            double key2)
{
    pastix_int_t i, hi;

    /* Allocate more space if necessary */
    if(q->size == q->used)
    {
        pastix_queue_item_t *tmp;
        tmp = q->elttab;

        assert( (q->size == 0) || (tmp != NULL) );

        /* OIMBE Realloc ?? */
        MALLOC_INTERN(q->elttab, q->size*2+1, pastix_queue_item_t);
        memcpy(q->elttab, tmp, q->size * sizeof(pastix_queue_item_t));

        q->size = q->size*2 +1;
        if (tmp != NULL)
            memFree_null(tmp);
    }

    q->elttab[q->used].key1   = key1;
    q->elttab[q->used].key2   = key2;
    q->elttab[q->used].eltptr = elt;
    q->used++;

    i = q->used - 1;
    hi= (i+1)/2-1;

    while( (i > 0) &&
           pqueueItemComparison(q, i, hi) )
    {
        pastix_queue_item_t swap = q->elttab[i];

        q->elttab[i ] = q->elttab[hi];
        q->elttab[hi] = swap;

        i = hi+1; hi = (i+1)/2-1;
    }
}

void *
pqueueRead(pastix_queue_t *q)
{
    return q->elttab[0].eltptr;
}

void *
pqueuePop2(pastix_queue_t *q, double *key1, double*key2)
{
    pastix_int_t i, j;
    void *return_elt;

    if (q->used == 0)
        return (void*)-1;

    return_elt = q->elttab[0].eltptr;
    if (key1 != NULL) *key1 = q->elttab[0].key1;
    if (key2 != NULL) *key2 = q->elttab[0].key2;

    q->elttab[0] = q->elttab[q->used-1];
    q->used--;

    i = 1;
    while(i <= (q->used/2))
    {
        if( (2*i == q->used)
            || pqueueItemComparison(q, 2*i-1, 2*i))     /*(q->keytab[2*i-1] < q->keytab[2*i]))*/
        {
            j = 2*i;
        }
        else
        {
            j = 2*i+1;
        }
        if (!pqueueItemComparison(q, i-1, j-1))         /*(q->keytab[i-1] >= q->keytab[j-1])*/
        {
            pastix_queue_item_t swap;

            swap           = q->elttab[i-1];
            q->elttab[i-1] = q->elttab[j-1];
            q->elttab[j-1] = swap;

            i = j;
        }
        else
            break;
    }
    return return_elt;
}

/*
  Function: QueueInit

  Allocate the queue array

  Parameters:
    q    - The queue to initialize
    size - The initial size of the queue.

  Return:
    PASTIX_SUCCESS - If all goes well.
*/
int queueInit(Queue *q,
              pastix_int_t    size)
{
    q->size = size;
    q->used = 0;
    if (q->size != 0)
      {
        MALLOC_INTERN(q->elttab,  size, pastix_int_t);
        MALLOC_INTERN(q->keytab,  size, double);
        MALLOC_INTERN(q->keytab2, size, pastix_int_t);
      }
    else
      {
        q->elttab  = NULL;
        q->keytab  = NULL;
        q->keytab2 = NULL;
      }
    return PASTIX_SUCCESS;
}
/*
  Function: queueCopy

  Perform a copy of a queue.

  Parameters:
    dst - Destination of the copy.
    src - Source of the copy.

  Returns:
    The copy address or NULL if the destination
    or the source is NULL.
*/
Queue * queueCopy(Queue *dst,
                  Queue *src)
{
  if(src == NULL || dst == NULL)
    return NULL;
  memcpy(dst, src, sizeof(Queue));
  MALLOC_INTERN(dst->elttab,  src->size, pastix_int_t);
  memcpy(dst->elttab, src->elttab, src->size * sizeof(pastix_int_t));

  MALLOC_INTERN(dst->keytab,  src->size, double);
  memcpy(dst->keytab, src->keytab, src->size * sizeof(double));

  MALLOC_INTERN(dst->keytab2, src->size, pastix_int_t);
  memcpy(dst->keytab2, src->keytab2, src->size * sizeof(pastix_int_t));

  return dst;
}

/*
  Function: queueExit

  Free a queue structure.

  Parameters:
    q - The queue to free.
*/
void queueExit(Queue *q)
{
  if(q->size != 0)
    {
      memFree_null(q->elttab);
      memFree_null(q->keytab);
      memFree_null(q->keytab2);
    }
  q->size = 0;
}

/*
  Function: queueAdd

  Add an element in a queue, following one
  double value as key.

  The element has to be a positive integer.

  Parameters:
    q   - the queue to fill.
    elt - the element to add.
    key - The key associated to the element.
*/
void queueAdd(Queue *q,
              pastix_int_t    elt,
              double key)
{
  queueAdd2(q, elt, key, 0);
}

/*
  Function: queueAdd2

  Add an element in a queue, following one
  double value and an integer value as keys.

  The element has to be a positive integer.

  Parameters:
    q    - the queue to fill.
    elt  - the element to add.
    key  - The double key associated to the element.
    key2 - The integer key associated to the element.
*/
void queueAdd2(Queue *q,
               pastix_int_t    elt,
               double key,
               pastix_int_t    key2)
{
    pastix_int_t i;
    pastix_int_t   swap_elt;
    double swap_key;
    pastix_int_t   swap_key2;
    pastix_int_t *tmp;
    double * tmp2;
    pastix_int_t   * tmp3;


    /** Allocate more space if necessary **/
    if(q->size == q->used)
        {
            tmp  = q->elttab;
            tmp2 = q->keytab;
            tmp3 = q->keytab2;
            /* OIMBE Realloc ?? */
            MALLOC_INTERN(q->elttab, q->size*2+1, pastix_int_t);
            memcpy(q->elttab, tmp, q->size * sizeof(pastix_int_t));

            MALLOC_INTERN(q->keytab, q->size*2+1, double);
            memcpy(q->keytab, tmp2, q->size * sizeof(double));

            MALLOC_INTERN(q->keytab2, q->size*2+1, pastix_int_t);
            memcpy(q->keytab2, tmp3, q->size * sizeof(pastix_int_t));

            q->size = q->size*2 +1;
            if (tmp != NULL)
              memFree_null(tmp);
            if (tmp2 != NULL)
              memFree_null(tmp2);
            if (tmp3 != NULL)
              memFree_null(tmp3);
        }

    q->elttab[q->used] = elt;
    q->keytab[q->used] = key;
    q->keytab2[q->used] = key2;
    q->used++;
    i = q->used;

    while( (i>1)
           &&  compWith2keys(q, i-1, i/2-1))
        {
            swap_elt = q->elttab[i-1];
            swap_key = q->keytab[i-1];
            swap_key2 = q->keytab2[i-1];
            q->elttab[i-1] = q->elttab[i/2-1];
            q->keytab[i-1] = q->keytab[i/2-1];
            q->keytab2[i-1] = q->keytab2[i/2-1];

            q->elttab[i/2-1] = swap_elt;
            q->keytab[i/2-1] = swap_key;
            q->keytab2[i/2-1] = swap_key2;
            i=i/2;
        }
}

/*
  Function: queueGet

  Get next element of the queue and
  remove it from the queue.

  Parameters:
    q - The queue from which user wants an element.

  Returns:
    The element if it was found, or -1 if it wasn't.
*/
pastix_int_t queueGet(Queue *q)
{
    pastix_int_t i, j;
    pastix_int_t return_elt;
    pastix_int_t swap_elt;
    double swap_key;
    pastix_int_t   swap_key2;

    if (q->used == 0)
      return -1;

    return_elt = q->elttab[0];

    q->elttab[0]  = q->elttab[q->used-1];
    q->keytab[0]  = q->keytab[q->used-1];
    q->keytab2[0] = q->keytab2[q->used-1];
    q->used--;

    i=1;

    while(i <= (q->used/2))
        {
            if( (2*i == q->used)
                || compWith2keys(q, 2*i-1, 2*i))     /*(q->keytab[2*i-1] < q->keytab[2*i]))*/
                {
                    j = 2*i;
                }
            else
                {
                    j = 2*i+1;
                }
            if (!compWith2keys(q, i-1, j-1))         /*(q->keytab[i-1] >= q->keytab[j-1])*/
                {
                    swap_elt = q->elttab[i-1];
                    swap_key = q->keytab[i-1];
                    swap_key2 = q->keytab2[i-1];

                    q->elttab[i-1]  = q->elttab[j-1];
                    q->keytab[i-1]  = q->keytab[j-1];
                    q->keytab2[i-1] = q->keytab2[j-1];

                    q->elttab[j-1] = swap_elt;
                    q->keytab[j-1] = swap_key;
                    q->keytab2[j-1] = swap_key2;

                    i=j;
                }
            else
                break;
        }
    return return_elt;
}


/*
  Function: queueSize

  Compute the size of a queue.

  Parameters:
    q - the queue.

  Returns:
    The size of the queue.
*/
pastix_int_t queueSize(Queue *q)
{
  return q->used;
}


void queueClear(Queue *q)
{
  q->used = 0;
}


/*
  Function: queueRead

  Read the next element that 'll be given by queueGet
  but not suppress it from the queue

  Parameters:
    q - The queue.

  Returns:
    The next element.
*/
pastix_int_t queueRead(Queue *q)
{
  return q->elttab[0];
}

/*
  Function: compwith2keys

  Compare 2 elements following their two keys.

  Parameters:
    q    - compare two elements
    elt1 - index of the first element in the queue.
    elt2 - index of the second element in the queue.
*/
pastix_int_t compWith2keys(Queue *q,
                  pastix_int_t    elt1,
                  pastix_int_t    elt2)
{
  /* if elt1 < elt2 return 1  */
  /* if elt1 = elt2 return 0  */
  /* if elt1 > elt2 return 0 */

  if(q->keytab[elt1] < q->keytab[elt2])
    return 1;
  if(q->keytab[elt1] > q->keytab[elt2])
    return 0;
  if(q->keytab2[elt1] < q->keytab2[elt2])
    return 1;
  return 0;
}
/*
  Function: queueGet2

  Get next element of the queue and remove it from the queue.

  Parameters:
    q    - The queue.
    key  - The first key (double) of the element.
    key2 - The second key (integer) of the element.

  Returns:
    The element, or -1 if not found.
*/
pastix_int_t queueGet2(Queue  *q,
              double *key,
              pastix_int_t    *key2)
{
    pastix_int_t i, j;
    pastix_int_t return_elt;
    pastix_int_t swap_elt;
    double swap_key;
    pastix_int_t   swap_key2;

    if (q->used == 0)
      return -1;

    return_elt = q->elttab[0];
    if (key  != NULL) (*key)  = q->keytab[0];
    if (key2 != NULL) (*key2) = q->keytab2[0];

    q->elttab[0]  = q->elttab[q->used-1];
    q->keytab[0]  = q->keytab[q->used-1];
    q->keytab2[0] = q->keytab2[q->used-1];
    q->used--;

    i=1;

    while(i <= (q->used/2))
        {
            if( (2*i == q->used)
                || compWith2keys(q, 2*i-1, 2*i))     /*(q->keytab[2*i-1] < q->keytab[2*i]))*/
                {
                    j = 2*i;
                }
            else
                {
                    j = 2*i+1;
                }
            if (!compWith2keys(q, i-1, j-1))         /*(q->keytab[i-1] >= q->keytab[j-1])*/
                {
                    swap_elt = q->elttab[i-1];
                    swap_key = q->keytab[i-1];
                    swap_key2 = q->keytab2[i-1];

                    q->elttab[i-1]  = q->elttab[j-1];
                    q->keytab[i-1]  = q->keytab[j-1];
                    q->keytab2[i-1] = q->keytab2[j-1];

                    q->elttab[j-1] = swap_elt;
                    q->keytab[j-1] = swap_key;
                    q->keytab2[j-1] = swap_key2;

                    i=j;
                }
            else
                break;
        }

    return return_elt;
}
/*
  Function: queuePosses

  Check if an element belongs to a queue.

  Parameters:
    q   - The queue.
    elt - The searched element.

  Returns:
    API_YES - if the element was found.
    API_NO  - if the element was not found.
*/
int queuePossess(Queue * q,
                 pastix_int_t     elt ){
  pastix_int_t i;
  for (i = 0; i < q->used; i++)
    if (q->elttab[i] == elt)
      return API_YES;

  return API_NO;
}

/*
  Function: queueInit

  Print a queue entries to standar error output.

  Parameters:
    q - The queue.
*/
void queuePrint(Queue *q)
{
  pastix_int_t i;
  fprintf(stderr, "Queue :");
  for (i = 0; i < q->used; i++)
    fprintf(stderr, "(%d %f %d) ",
            (int)q->elttab[i],
            (double)q->keytab[i],
            (int)q->keytab2[i] );
  fprintf(stderr, "\n");
}
