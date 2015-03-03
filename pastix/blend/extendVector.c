#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "common.h"
#include "extendVector.h"

pastix_int_t * extendint_Init(ExtendVectorINT *vec, pastix_int_t size)
{
    vec->vecsize = size;
    vec->eltnbr  = 0;
    vec->inttab  = NULL;
    MALLOC_INTERN(vec->inttab, size, pastix_int_t);
    return vec->inttab;
}


void extendint_Exit(ExtendVectorINT *vec)
{
    if(vec->inttab != NULL)
        memFree_null(vec->inttab);
    /*memFree_null(vec);*/
}

void extendint_Add(ExtendVectorINT *vec, pastix_int_t elt)
{
    vec->inttab[vec->eltnbr] = elt;
    extendint_incr(vec);
}

pastix_int_t extendint_Size(ExtendVectorINT *vec)
{
  return vec->eltnbr;
}

pastix_int_t extendint_Read(ExtendVectorINT *vec, pastix_int_t eltnum)
{
  ASSERT(eltnum <= vec->eltnbr,MOD_BLEND);
  return vec->inttab[eltnum];
}



void extendint_ToSize(pastix_int_t size, ExtendVectorINT *vec)
{
    extendint_Clear(vec);

    if(size <= vec->vecsize)  /* there 's enough space */
        return;


    if(vec->inttab != NULL)
        memFree_null(vec->inttab);

    MALLOC_INTERN(vec->inttab, size, pastix_int_t);
    vec->vecsize = size;
}

void extendint_incr(ExtendVectorINT *vec)
{
    vec->eltnbr++;
    /** if the vector is not big enough, make it bigger !! **/
    if(!(vec->eltnbr < vec->vecsize))
        {
            pastix_int_t *tmp;
            tmp = vec->inttab;
            /* add memory space */
            MALLOC_INTERN(vec->inttab, vec->vecsize + vec->vecsize/2 +1, pastix_int_t);
            memcpy(vec->inttab, tmp, sizeof(pastix_int_t)*vec->eltnbr);
            vec->vecsize = vec->vecsize + vec->vecsize/2 +1;
            memFree_null(tmp);
        }
}

void extendint_Clear(ExtendVectorINT *vec)
{
    vec->eltnbr = 0;
}
