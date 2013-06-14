#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "common_pastix.h"
/* #include "symbol.h" */

#include "sparRow.h"
#include "sort_row.h"

void sort_row(csptr P)
{
  /***********************************************************************************************/
  /* This function sort all the row of a symbolic matrix (just the pattern) by ascending indices */
  /* This is done in place                                                                       */
  /***********************************************************************************************/
  PASTIX_INT i;
  
  for(i=0;i<P->n;i++)
    if(P->nnzrow[i] > 1)
      intSort1asc1(P->ja[i], P->nnzrow[i]);
}
