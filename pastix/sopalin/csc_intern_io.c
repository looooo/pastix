/*
  File: csc_intern_io.c

  Functions to save or load internal CSC in binary or ascii mode.
  
*/
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#include "common.h"
#include "csc.h"
#include "csc_intern_io.h"


#include "ftgt.h"
#include "updown.h"
#include "queue.h"
#include "bulles.h"
#include "solver.h"
/*
  Function: CscSave

  Writes on disk an internal CSCd in text format.

  Format is :
  
  > CSC_FNBR(cscptr)
  > CSC_COLNBR(cscptr,iter)    ! iter = 0 to CSC_FNBR(cscptr) - 1
  > CSC_COL(cscptr,iter,iter2) ! iter2 = 0 to CSC_COLNBR(cscptr,iter)
  > ...
  > CSC_ROW(cscptr,iter) ! For all rows and values (iter)
  > CSC_VAL(cscptr,iter)

  Parameters :
    cscprt - the internal CSCd structure to save.
    stream - the FILE to write into, open in write mode.
*/
PASTIX_INT CscSave(const CscMatrix * const cscptr, 
	    FILE            * const stream)
{
  PASTIX_INT iter=0;
  PASTIX_INT iter2=0;
  PASTIX_INT valnbr;
  PASTIX_INT o=0;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscSave \n");
#endif

  fprintf(stream, "%ld\n", (long)CSC_FNBR(cscptr));

  for (iter=0; iter < CSC_FNBR(cscptr); iter++)
    {
      fprintf(stream, "%ld\n", (long)CSC_COLNBR(cscptr,iter));

      for (iter2=0; iter2<CSC_COLNBR(cscptr,iter)+1; iter2++)
	{
	  fprintf(stream, "%ld\n", (long)CSC_COL(cscptr,iter,iter2));
	}
    }

  if (CSC_FNBR(cscptr) > 0)
    {
      valnbr = CSC_VALNBR(cscptr);
      
      for (iter=0; iter < valnbr ;iter++)
	{
	  fprintf(stream, "%ld\n", (long)CSC_ROW(cscptr,iter));
#ifdef CPLX
	  fprintf(stream, "%e %e\n", creal(CSC_VAL(cscptr,iter)), cimag(CSC_VAL(cscptr,iter)));
#else
	  fprintf(stream, "%e\n", CSC_VAL(cscptr,iter));
#endif
	}
    }

#ifdef CSC_LOG
  fprintf(stdout, "<- CscSave \n");
#endif
  
  return o;
}

#define CscSaveIJV PASTIX_EXTERN_F(CscSaveIJV)
PASTIX_INT CscSaveIJV(const CscMatrix * const cscptr,
	       const SolverMatrix     *solvmtx,
	       PASTIX_INT                    *l2g,
	       PASTIX_INT                    *peritab,
	       PASTIX_INT                     dof,
	       FILE            * const stream)
{
  PASTIX_INT itercblk;
  PASTIX_INT itercol;
  PASTIX_INT iterval, indcol, colstart, colend;
  PASTIX_INT indcblk;
  PASTIX_INT o=0;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscSave \n");
#endif

  for (itercblk=0; itercblk<CSC_FNBR(cscptr); itercblk++)
    {
      indcblk = solvmtx->cblktab[itercblk].fcolnum;
      
      for (itercol=0; itercol<CSC_COLNBR(cscptr,itercblk); itercol++)
	{
	  colstart = CSC_COL(cscptr,itercblk,itercol);
	  colend   = CSC_COL(cscptr,itercblk,itercol+1);
	  indcol   = indcblk+itercol;
	      
	  for (iterval=colstart; iterval<colend; iterval++)
	    {
	      fprintf(stream, "%ld ", (long)( peritab[(CSC_ROW(cscptr,iterval) - 
						       CSC_ROW(cscptr,iterval)%dof)/dof]*dof + 1
					      + CSC_ROW(cscptr,iterval)%dof));
	      fprintf(stream, "%ld ", (long)((l2g[peritab[(indcol- indcol%dof)/dof]]-1)*dof+1+indcol%dof));
#ifdef CPLX
	      fprintf(stream, "%e %e\n", creal(CSC_VAL(cscptr,iterval)), cimag(CSC_VAL(cscptr,iterval)));
#else
	      fprintf(stream, "%e\n", CSC_VAL(cscptr,iterval));
#endif
	    }
	}
    }

#ifdef CSC_LOG
  fprintf(stdout, "<- CscSaveIJV \n");
#endif
  
  return o;
}
/*
  Function: CscBSave

  Writes on disk an internal CSCd in binary format.

  Parameters :
    cscprt - the internal CSCd structure to save.
    stream - the FILE to write into, open in write mode.
*/
PASTIX_INT CscBSave(const CscMatrix * const cscptr, 
	     FILE            * const stream)
{
  PASTIX_INT iter=0;
  PASTIX_INT valnbr;
  PASTIX_INT o=0;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscBSave \n");
#endif
  
  fwrite(&(CSC_FNBR(cscptr)),sizeof(PASTIX_INT), 1, stream);

  for (iter=0; iter < CSC_FNBR(cscptr); iter++)
    {
      fwrite(&(CSC_COLNBR(cscptr,iter)), sizeof(PASTIX_INT), 1, stream); 

      fwrite(CSC_COLTAB(cscptr,iter), sizeof(PASTIX_INT),
	     (CSC_COLNBR(cscptr,iter)+1), stream);
    }

  valnbr = CSC_VALNBR(cscptr);

  for (iter=0; iter < valnbr ;iter++)
    {
      fwrite(&(CSC_ROW(cscptr,iter)), sizeof(PASTIX_INT), 1, stream);
      fwrite(&(CSC_VAL(cscptr,iter)), sizeof(PASTIX_FLOAT), 1,stream);
    }

#ifdef CSC_LOG
  fprintf(stdout, "<- CscBSave \n");
#endif
  
  return o;
}

/* 
   Function: CscLoad

   Reads an internal CSCd from disk.

   Format is :
   
   > CSC_FNBR(cscptr)
   > CSC_COLNBR(cscptr,iter)    ! iter = 0 to CSC_FNBR(cscptr) - 1
   > CSC_COL(cscptr,iter,iter2) ! iter2 = 0 to CSC_COLNBR(cscptr,iter)
   > ...
   > CSC_ROW(cscptr,iter) ! For all rows and values (iter)
   > CSC_VAL(cscptr,iter)

   Parameters :
     cscprt - the internal CSCd structure to load.
     stream - the FILE to write into, open in read mode. 
*/
PASTIX_INT CscLoad(CscMatrix * cscptr, 
	    FILE      * stream)
{
  PASTIX_INT iter=0;
  PASTIX_INT iter2=0;
  PASTIX_INT valnbr;
  long temp;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscLoad \n");
#endif
  
  if (1 != fscanf(stream, "%ld\n", &temp)){
    errorPrint("CSC badly formated");
    return EXIT_FAILURE;
  }

  CSC_FNBR(cscptr) = temp;
  
  MALLOC_INTERN(CSC_FTAB(cscptr), CSC_FNBR(cscptr), CscFormat);
  
  for (iter=0; iter < CSC_FNBR(cscptr); iter++)
    {
      if (1 != fscanf(stream, "%ld\n", &temp)){
	errorPrint("CSC badly formated");
	return EXIT_FAILURE;
      }
  
      CSC_COLNBR(cscptr,iter)=temp;
      
      MALLOC_INTERN(CSC_COLTAB(cscptr,iter), 
		    CSC_COLNBR(cscptr,iter)+1, 
		    PASTIX_INT);
      
      for (iter2=0; iter2<CSC_COLNBR(cscptr,iter)+1; iter2++)
	{
	  if (1 != fscanf(stream, "%ld\n", &temp)){
	    errorPrint("CSC badly formated");
	    return EXIT_FAILURE;
	  }
	  
	  CSC_COL(cscptr,iter,iter2)=temp;
	}
      
    }
  
  valnbr = CSC_VALNBR(cscptr);

  MALLOC_INTERN(CSC_ROWTAB(cscptr), valnbr, PASTIX_INT);
  MALLOC_INTERN(CSC_VALTAB(cscptr), valnbr, PASTIX_FLOAT);

  for (iter=0; iter < valnbr; iter++)
    {
      if (1 != fscanf(stream, "%ld\n", &temp)){
	errorPrint("CSC badly formated");
	return EXIT_FAILURE;
      }

      CSC_ROW(cscptr,iter)=temp;
#ifdef CPLX
      {
	double tempreal, tempimag;
	if (2 != fscanf(stream, "%lf %lf\n", &tempreal, &tempimag)){
	  errorPrint("CSC badly formated");
	  return EXIT_FAILURE;
	}

#if (defined X_ARCHalpha_compaq_osf1)
	CSC_VAL(cscptr,iter) = PASTIX_FLOAT (tempreal, tempimag);
#else
	CSC_VAL(cscptr,iter) = (PASTIX_FLOAT) tempreal+( (PASTIX_FLOAT) tempimag)*I;
#endif
      }
#else
      if (1 != fscanf(stream, "%lf\n", (double *)&(CSC_VAL(cscptr,iter)))){
	errorPrint("CSC badly formated");
	return EXIT_FAILURE;
      }
      
#endif
    }

#ifdef CSC_LOG
  fprintf(stdout, "<- CscLoad \n");
#endif
  
  return 0;
}

/*
  Function: CscBLoad
  
  Loads an internal CSCd from a file saved in binary mode.

  Parameters :
    cscprt - the internal CSCd structure to load.
    stream - the FILE to write into, open in read mode. 
*/
PASTIX_INT CscBLoad(CscMatrix * cscptr, 
	     FILE      * stream)
{
  PASTIX_INT iter=0;
  PASTIX_INT valnbr;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscBLoad \n");
#endif

  PASTIX_FREAD(&(CSC_FNBR(cscptr)), sizeof(PASTIX_INT), 1, stream);

  MALLOC_INTERN(CSC_FTAB(cscptr), CSC_FNBR(cscptr), CscFormat);

  for (iter=0; iter < CSC_FNBR(cscptr); iter++)
    {
      PASTIX_FREAD(&(CSC_COLNBR(cscptr,iter)), sizeof(PASTIX_INT), 1, stream);

      MALLOC_INTERN(CSC_COLTAB(cscptr,iter),
		    CSC_COLNBR(cscptr,iter)+1, 
		    PASTIX_INT);

      PASTIX_FREAD(CSC_COLTAB(cscptr,iter), sizeof(PASTIX_INT),
	    (CSC_COLNBR(cscptr,iter)+1), stream);
    }

  valnbr = CSC_VALNBR(cscptr);

  MALLOC_INTERN(CSC_ROWTAB(cscptr), valnbr, PASTIX_INT);
  MALLOC_INTERN(CSC_VALTAB(cscptr), valnbr, PASTIX_FLOAT);

  for (iter=0; iter < valnbr; iter++)
    {
      PASTIX_FREAD(&(CSC_ROW(cscptr,iter)), sizeof(PASTIX_INT), 1, stream);
      PASTIX_FREAD(&(CSC_VAL(cscptr,iter)), sizeof(PASTIX_FLOAT), 1, stream);
    }

#ifdef CSC_LOG
  fprintf(stdout, "<- CscBLoad \n");
#endif
  return 0;
}
