#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "common_pastix.h"
#include "extendVector.h"

PASTIX_INT * extendint_Init(ExtendVectorINT *vec, PASTIX_INT size)
{
    vec->vecsize = size;
    vec->eltnbr  = 0;
    vec->inttab  = NULL;
    MALLOC_INTERN(vec->inttab, size, PASTIX_INT);
    return vec->inttab;
}

    
void extendint_Exit(ExtendVectorINT *vec)
{
    if(vec->inttab != NULL)
	memFree_null(vec->inttab);
    /*memFree_null(vec);*/
}

void extendint_Add(ExtendVectorINT *vec, PASTIX_INT elt)
{
    vec->inttab[vec->eltnbr] = elt;
    extendint_incr(vec);
}

PASTIX_INT extendint_Size(ExtendVectorINT *vec)
{
  return vec->eltnbr;
}

PASTIX_INT extendint_Read(ExtendVectorINT *vec, PASTIX_INT eltnum)
{
  ASSERT(eltnum <= vec->eltnbr,MOD_BLEND);
  return vec->inttab[eltnum];
}



void extendint_ToSize(PASTIX_INT size, ExtendVectorINT *vec)
{
    extendint_Clear(vec);

    if(size <= vec->vecsize)  /* there 's enough space */
	return;
    

    if(vec->inttab != NULL)   
	memFree_null(vec->inttab);

    MALLOC_INTERN(vec->inttab, size, PASTIX_INT);
    vec->vecsize = size;
}
    
void extendint_incr(ExtendVectorINT *vec)
{
    vec->eltnbr++;
    /** if the vector is not big enough, make it bigger !! **/
    if(!(vec->eltnbr < vec->vecsize))
	{
	    PASTIX_INT *tmp;
	    tmp = vec->inttab;
	    /* add memory space */
	    MALLOC_INTERN(vec->inttab, vec->vecsize + vec->vecsize/2 +1, PASTIX_INT);
	    memcpy(vec->inttab, tmp, sizeof(PASTIX_INT)*vec->eltnbr);
	    vec->vecsize = vec->vecsize + vec->vecsize/2 +1;
	    memFree_null(tmp);
	}
}

void extendint_Clear(ExtendVectorINT *vec)
{
    vec->eltnbr = 0;
}
