#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <stdint.h>

#ifdef FORCE_NOMPI
#include "pastix_nompi.h"
#else
#include <mpi.h>
#endif
#include <stdlib.h>

#ifdef X_ARCHsun
#include <inttypes.h>
#endif

#include "pastix.h"
#include "common_drivers.h"

#ifdef __INTEL_COMPILER
/* Ignore icc remark : "external declaration in primary source file" */
#pragma warning(disable:1419)
#endif

/* Trouver une solution plus propre pour cette 
   déclarer fonction interne de la libpastix */

#include<ctype.h>

/*
  Function: myupcase

  Rewrites *s* to upper case.
  
  Parameters: 
    s - string to rexwrite in upcase.
*/
void myupcase(char *S)
{
  pastix_int_t iter=0;

  while (S[iter] != '\0')
    {
      S[iter] = (char)toupper(S[iter]);
      iter++;
    }
}


/*
  Function:  mysubstr

  Copy len element, from *S[pos]* in *s*.
  
  Parameters:
    s   - destination
    S   - Source
    pos - sarting position
    len - Size of the string to copy.
*/
void mysubstr(char *s, const char *S, const pastix_int_t pos, const pastix_int_t len)
{
  pastix_int_t iter;
  for (iter=0; iter<len;iter++)
    {
      s[iter] = S[pos+iter];
    }
  s[len] = '\0';
}

/*
  Function:  mysubstr2

  Copy the number placed between a and b in fmt.
  
  Parameters:
    fmt - String in which there is a and b
    a   - first element
    b   - last element
    val - the integer between a and b

*/
void mysubstr2(const char *fmt, const char a, const char b, pastix_int_t *val)
{
  char * posb = strchr(fmt,b);
  char * posa = strchr(fmt,a);
  int len = posb - posa - 1;
  char *tmp = (char *) malloc(len+1);
  if (tmp == NULL)
    {
      fprintf(stderr, "mysubstr2 : Not enough memory for tmp\n");
      EXIT(MOD_SI,OUTOFMEMORY_ERR);
    }
  mysubstr(tmp, fmt, strchr(fmt, a) - fmt +1, len);
  *val = atoi(tmp);
  memFree_null(tmp);
}
