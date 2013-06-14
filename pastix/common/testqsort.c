#include "common_pastix.h"
#include "errors.h"

#define NB_ELEMENTS 10000
int main(int argc, char ** argv)
{
  int i,j;
  PASTIX_INT   * index;
  PASTIX_FLOAT * vals;
  PASTIX_INT   * index2;
  PASTIX_FLOAT * vals2;
  Clock  clk;
  void  * ptr[2];

  if (NULL == (index  = memAlloc(sizeof(PASTIX_INT)*NB_ELEMENTS)))
    MALLOC_ERROR("index");
  if (NULL == (index2  = memAlloc(sizeof(PASTIX_INT)*NB_ELEMENTS)))
    MALLOC_ERROR("index2");
  if (NULL == (vals  = memAlloc(sizeof(PASTIX_FLOAT)*NB_ELEMENTS)))
    MALLOC_ERROR("vals");
  if (NULL == (vals2  = memAlloc(sizeof(PASTIX_FLOAT)*NB_ELEMENTS)))
    MALLOC_ERROR("vals2");

  memcpy(index2, index, sizeof(PASTIX_INT)*NB_ELEMENTS);
  memcpy(vals2, vals, sizeof(PASTIX_FLOAT)*NB_ELEMENTS);

  clockInit(&clk);
  clockStart(&clk);
  ptr[0] = index;
  ptr[1] = vals; 
  qsortIntFloatAsc(ptr,NB_ELEMENTS);

  clockStop (&clk);
  for (i = 0; i < NB_ELEMENTS -1; i++)
    if (index[i] > index[i+1])
      {
	errorPrint("Mauvais tri 1 : mal ordonné");
	return EXIT_FAILURE;
      }
  
  for (i = 0; i < NB_ELEMENTS; i++)
    {
      j = 0;
      while (index[j] != index2[i] && j < NB_ELEMENTS)
	j++;
      while (index[j] == index2[i] && vals[j] != vals2[i] && j < NB_ELEMENTS)
	j++;

      if (j == NB_ELEMENTS)
	{
	  errorPrint("Mauvais tri 1: élément disparu");
	  return EXIT_FAILURE;
	}
      if (index[j] != index2[i] || vals[j] != vals2[i])
	{
	  errorPrint("Mauvais tri 1 : mauvaise recopie de valeur associee a l'index");
	  return EXIT_FAILURE;
	}
    }
  fprintf(stdout,"tri de %ld elements aleatoires correct en %.3g s\n", (long)NB_ELEMENTS, (float)clockVal(&clk));


  /* Sort sorted elements */
  for (i = 0; i < NB_ELEMENTS; i++)
    {
      index[i] = i;
      vals[i] = (PASTIX_FLOAT)random();
    }

  memcpy(index2, index, sizeof(PASTIX_INT)*NB_ELEMENTS);
  memcpy(vals2, vals, sizeof(PASTIX_FLOAT)*NB_ELEMENTS);
  clockInit(&clk);
  clockStart(&clk);
  ptr[0] = index;
  ptr[1] = vals; 
  qsortIntFloatAsc(ptr,NB_ELEMENTS);
  clockStop (&clk);
  for (i = 0; i < NB_ELEMENTS -1; i++)
    if (index[i] > index[i+1])
      {
	errorPrint("Mauvais tri 2 : mal ordonné");
	return EXIT_FAILURE;
      }
  
  for (i = 0; i < NB_ELEMENTS; i++)
    {
      j = 0;
      while (index[j] != index2[i] && j < NB_ELEMENTS)
	j++;
      while (index[j] == index2[i] && vals[j] != vals2[i] && j < NB_ELEMENTS)
	j++;

      if (j == NB_ELEMENTS)
	{
	  errorPrint("Mauvais tri 2 : élément disparu");
	  return EXIT_FAILURE;
	}
      if (index[j] != index2[i] || vals[j] != vals2[i])
	{
	  errorPrint("Mauvais tri 2 : mauvaise recopie de valeur associee a l'index");
	  return EXIT_FAILURE;
	}
    }
  fprintf(stdout,"tri de %ld elements triés correct en %.3g s \n", (long)NB_ELEMENTS, (float)clockVal(&clk));

  memFree_null(index);
  memFree_null(index2);
  memFree_null(vals);
  memFree_null(vals2);

  return EXIT_SUCCESS;
}

