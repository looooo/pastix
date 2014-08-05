/**
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/
/*
 * File: debug_dump.c
 *
 * Functions to dump informations on disk.
 */
#include "common.h"
#include "order.h"
#include "z_csc.h"
#include "symbol.h"
#include "z_ftgt.h"
#include "z_updown.h"
#include "queue.h"
#include "bulles.h"
#include "z_solver.h"
#include "sopalin_acces.h"

#ifdef DEBUG_RAFF
#define CSC_LOG
#endif

/*
 * Function: z_dump1
 *
 * Dumps ord->permtab on disk.
 *
 * Format:
 *   > i -> permtab[i]
 *
 * parameters:
 *   ord    - Order structure to print permtab from.
 *   stream - FILE *, opened in write mode, in which permtab will be writen.
 *   colnbr - Number of elements in permtab.
 */
void z_dump1(Order *ord,
           FILE  *stream,
           pastix_int_t    colnbr)
{
  pastix_int_t iter;

#ifdef CSC_LOG
  fprintf(stdout, "-> z_dump1 \n");
#endif

  for (iter=0; iter<colnbr; iter++)
  {
    fprintf(stream, "%ld -> %ld\n", (long)iter, (long)ord->permtab[iter]);
  }

#ifdef CSC_LOG
  fprintf(stdout, "<- z_dump1 \n");
#endif
}

/*
 Function: z_dump2

 Prints internal CSCd, in (i,j,v) format, in a file.

 Parameters:
 datacode - z_SolverMatrix.
 stream   - FILE * opened in write mode.

 */
void z_dump2(const z_SolverMatrix * datacode,
           const z_CscMatrix    * cscmtx,
           pastix_complex64_t              *trandcsc,
           FILE               *stream)
{
  /* Transforme une csc en ij valeur */
  pastix_int_t itercblk;
  pastix_int_t itercoltab;
  pastix_int_t iterval;
  pastix_int_t itercol;

#ifdef CSC_LOG
  fprintf(stdout, "-> z_dump2 \n");
#endif

  for (itercblk = 0; itercblk < CSC_FNBR(cscmtx); itercblk++)
  {
    itercol = SYMB_FCOLNUM(itercblk);
    for (itercoltab=0;
         itercoltab < cscmtx->cscftab[itercblk].colnbr;
         itercoltab++)
    {
      for (iterval = cscmtx->cscftab[itercblk].coltab[itercoltab];
           iterval < cscmtx->cscftab[itercblk].coltab[itercoltab+1];
           iterval++)
      {
        if (ABS_FLOAT(cscmtx->valtab[iterval]) != 0)
        {
          if (trandcsc == NULL)
          {
#ifdef TYPE_COMPLEX
            fprintf(stream, "%ld %ld (%10g,%10g)\n",
                    (long)itercol, (long)cscmtx->rowtab[iterval],
                    creal(cscmtx->valtab[iterval]), cimag(cscmtx->valtab[iterval]));
#else
            fprintf(stream, "%ld %ld %10g\n",
                    (long)itercol, (long)cscmtx->rowtab[iterval],
                    cscmtx->valtab[iterval]);
#endif
          }
          else
          {
#ifdef TYPE_COMPLEX
            fprintf(stream, "%ld %ld (%10g,%10g) (%10g,%10g)\n",
                    (long)itercol, (long)cscmtx->rowtab[iterval],
                    creal(cscmtx->valtab[iterval]), cimag(cscmtx->valtab[iterval]),
                    creal(trandcsc[iterval]), cimag(trandcsc[iterval]));
#else
            fprintf(stream, "%ld %ld %10g %10g\n",
                    (long)itercol, (long)cscmtx->rowtab[iterval],
                    cscmtx->valtab[iterval],
                    trandcsc[iterval]);
#endif
          }
        }
      }
      itercol++;
    }
  }

#ifdef CSC_LOG
  fprintf(stdout, "<- z_dump2 \n");
#endif
}

/*
 * Function: z_dump3
 *
 * Prints z_solver matrix informations, in (i,j,v) format, in a file,
 * for LLt or LDLt decomposition.
 *
 * Parameters:
 *   datacode - z_SolverMatrix.
 *   stream   - FILE * opened in write mode.
 */
void z_dump3(const z_SolverMatrix *datacode,
           FILE               *stream)
{
  pastix_int_t itercblk;
  pastix_int_t itercol;
  pastix_int_t iterbloc;
  pastix_int_t iterrow;
  pastix_int_t coefindx;
  /*   z_SolverMatrix * datacode = sopalin_data->datacode; */
#ifdef CSC_LOG
  fprintf(stdout, "-> z_dump3 \n");
#endif
  for (itercblk = 0; itercblk < SYMB_CBLKNBR; itercblk++)
  {
    for (itercol = SYMB_FCOLNUM(itercblk);
         itercol < SYMB_LCOLNUM(itercblk)+1;
         itercol++)
    {
      /* bloc diag */
      iterbloc = SYMB_BLOKNUM(itercblk);

      coefindx = SOLV_COEFIND(iterbloc);

      coefindx +=
        (itercol-SYMB_FCOLNUM(itercblk))*
        SOLV_STRIDE(itercblk);

         for (iterrow = SYMB_FROWNUM(iterbloc);
              iterrow < SYMB_LROWNUM(iterbloc)+1;
           iterrow++)
      {
        if ((ABS_FLOAT(SOLV_COEFTAB(itercblk)[coefindx]) != 0) &&
            (itercol<=iterrow))
        {
#ifdef TYPE_COMPLEX
          fprintf(stream, "%ld %ld (%13e,%13e)\n",
                  (long)itercol, (long)iterrow,
                  creal(SOLV_COEFTAB(itercblk)[coefindx]), cimag(SOLV_COEFTAB(itercblk)[coefindx]));
#else
          fprintf(stream, "%ld %ld %13e\n",
                  (long)itercol, (long)iterrow,
                  SOLV_COEFTAB(itercblk)[coefindx]);
#endif
        }
        coefindx++;
      }

      /* extra diag bloc */
         for (iterbloc = SYMB_BLOKNUM(itercblk)+1;
              iterbloc < SYMB_BLOKNUM(itercblk+1);
              iterbloc++)
      {
        coefindx = SOLV_COEFIND(iterbloc);

        coefindx +=
          (itercol-SYMB_FCOLNUM(itercblk))*
          datacode->cblktab[itercblk].stride;

        for (iterrow = SYMB_FROWNUM(iterbloc);
             iterrow < SYMB_LROWNUM(iterbloc)+1;
             iterrow++)
        {
          if (ABS_FLOAT(SOLV_COEFTAB(itercblk)[coefindx]) != 0)
          {
#ifdef TYPE_COMPLEX
            fprintf(stream, "%ld %ld (%13e,%13e)\n",
                    (long)itercol, (long)iterrow,
                    creal(SOLV_COEFTAB(itercblk)[coefindx]),cimag(SOLV_COEFTAB(itercblk)[coefindx]));
#else
            fprintf(stream, "%ld %ld %13e\n",
                    (long)itercol, (long)iterrow,
                    SOLV_COEFTAB(itercblk)[coefindx]);
#endif
          }
          /*
           if (SOLV_UCOEFTAB(itercblk)[coefindx] != 0)
           {
           #ifdef TYPE_COMPLEX
           fprintf(stream, "%ld %ld (%13e,%13e)\n",
           (long)iterrow, (long)itercol,
           creal(SOLV_UCOEFTAB(itercblk)[coefindx]),cimag(SOLV_UCOEFTAB(itercblk)[coefindx]));
           #else
           fprintf(stream, "%ld %ld %13e\n",
           (long)iterrow, (long)itercol,
           SOLV_UCOEFTAB(itercblk)[coefindx]);
           #endif
           }*/

          coefindx++;
        }
      }
    }
  }

#ifdef CSC_LOG
  fprintf(stdout, "<- z_dump3 \n");
#endif
}



/*
 * Function: z_dump3_LU
 *
 * Prints z_solver matrix informations, in (i,j,v) format, in a file,
 * for LU decomposition.
 *
 * Parameters:
 *   datacode - z_SolverMatrix.
 *   streamL  - FILE * opened in write mode.
 *   streamU  - FILE * opened in write mode.
 */
void z_dump3_LU(const z_SolverMatrix * datacode,
                FILE               * streamL,
                FILE               * streamU)
{
  pastix_int_t itercblk;
  pastix_int_t itercol;
  pastix_int_t iterbloc;
  pastix_int_t iterrow;
  pastix_int_t coefindx;

#ifdef CSC_LOG
  fprintf(stdout, "-> z_dump3 (LU)\n");
#endif
  for (itercblk = 0; itercblk < SYMB_CBLKNBR; itercblk++)
  {

    for (itercol = SYMB_FCOLNUM(itercblk);
         itercol < SYMB_LCOLNUM(itercblk)+1;
         itercol++)
    {
      /* bloc diag */
      iterbloc = SYMB_BLOKNUM(itercblk);

      coefindx = SOLV_COEFIND(iterbloc);

      coefindx +=
        (itercol-SYMB_FCOLNUM(itercblk))*
        datacode->cblktab[itercblk].stride;

      for (iterrow = SYMB_FROWNUM(iterbloc);
           iterrow < SYMB_LROWNUM(iterbloc)+1;
           iterrow++)
      {
        /* for L */
        if ((ABS_FLOAT(SOLV_COEFTAB(itercblk)[coefindx]) != 0) &&
            (itercol<=iterrow))
        {
#ifdef TYPE_COMPLEX
          fprintf(streamL, "%ld %ld (%13e,%13e)\n",
                  (long)itercol, (long)iterrow,
                  creal(SOLV_COEFTAB(itercblk)[coefindx]), cimag(SOLV_COEFTAB(itercblk)[coefindx]));
#else
          fprintf(streamL, "%ld %ld %13e\n",
                  (long)itercol, (long)iterrow,
                  SOLV_COEFTAB(itercblk)[coefindx]);
#endif
        }

        /* for U */
        if ((ABS_FLOAT(SOLV_UCOEFTAB(itercblk)[coefindx]) != 0) &&
            (itercol<iterrow)) /* attention ... ici strict */
        {
#ifdef TYPE_COMPLEX
          fprintf(streamU, "%ld %ld (%13e,%13e)\n",
                  (long)iterrow, (long)itercol,
                  creal(SOLV_UCOEFTAB(itercblk)[coefindx]), cimag(SOLV_UCOEFTAB(itercblk)[coefindx]));
#else
          fprintf(streamU, "%ld %ld %13e\n",
                  (long)iterrow, (long)itercol,
                  SOLV_UCOEFTAB(itercblk)[coefindx]);
#endif
        }
        coefindx++;
      }

      /* extra diag bloc */
      for (iterbloc = SYMB_BLOKNUM(itercblk)+1;
           iterbloc < SYMB_BLOKNUM(itercblk+1);
           iterbloc++)
      {
        coefindx = SOLV_COEFIND(iterbloc);

        coefindx +=
          (itercol-SYMB_FCOLNUM(itercblk))*
          datacode->cblktab[itercblk].stride;

        for (iterrow = SYMB_FROWNUM(iterbloc);
             iterrow < SYMB_LROWNUM(iterbloc)+1;
             iterrow++)
        {
          /* for L */
          if (ABS_FLOAT(SOLV_COEFTAB(itercblk)[coefindx]) != 0)
          {
#ifdef TYPE_COMPLEX
            fprintf(streamL, "%ld %ld (%13e,%13e)\n",
                    (long)itercol, (long)iterrow,
                    creal(SOLV_COEFTAB(itercblk)[coefindx]),cimag(SOLV_COEFTAB(itercblk)[coefindx]));
#else
            fprintf(streamL, "%ld %ld %13e\n",
                    (long)itercol, (long)iterrow,
                    SOLV_COEFTAB(itercblk)[coefindx]);
#endif
          }

          /* for U */
          if (ABS_FLOAT(SOLV_UCOEFTAB(itercblk)[coefindx]) != 0)
          {
#ifdef TYPE_COMPLEX
            fprintf(streamU, "%ld %ld (%13e,%13e)\n",
                    (long)iterrow, (long)itercol,
                    creal(SOLV_UCOEFTAB(itercblk)[coefindx]),cimag(SOLV_UCOEFTAB(itercblk)[coefindx]));

#else
            fprintf(streamU, "%ld %ld %13e\n",
                    (long)iterrow, (long)itercol,
                    SOLV_UCOEFTAB(itercblk)[coefindx]);
#endif
          }

          coefindx++;
        }
      }
    }
  }
#ifdef CSC_LOG
  fprintf(stdout, "<- z_dump3 (LU)\n");
#endif
}

/*
 * Function: z_dump4
 *
 * Writes column blocks and blocs dimension in a file.
 *
 * Parameters:
 *   datacode - z_SolverMatrix containing informations about blocs
 *   stream   - FILE * opened in write mode.
 */
void z_dump4(const z_SolverMatrix *datacode,
           FILE               *stream)
{
  pastix_int_t itercblk;
  pastix_int_t iterbloc;
  pastix_int_t itercolc;
  pastix_int_t itercola;

#ifdef CSC_LOG
  fprintf(stdout, "-> z_dump4 \n");
#endif
  for (itercblk=0; itercblk<SYMB_CBLKNBR; itercblk++)
  {
    fprintf(stream, "cblk %ld: %ld - %ld\n", (long)itercblk,
            (long)SYMB_FCOLNUM(itercblk),
            (long)SYMB_LCOLNUM(itercblk));

    for (iterbloc = SYMB_BLOKNUM(itercblk);
         iterbloc < SYMB_BLOKNUM(itercblk+1);
         iterbloc++)
    {
      fprintf(stream, "bloc %ld: %ld - %ld\n", (long)iterbloc,
              (long)SYMB_FROWNUM(iterbloc),
              (long)SYMB_LROWNUM(iterbloc));
    }
  }

  itercolc = 0;
  itercola = -1;

  /* cblk en continu */
  for (itercblk=0; itercblk<SYMB_CBLKNBR; itercblk++)
  {
    if (itercola == -1)
    {
      itercola = SYMB_FCOLNUM(itercblk);
      itercolc = SYMB_LCOLNUM(itercblk)+1;
    }
    else
    {
      if (itercolc == SYMB_FCOLNUM(itercblk))
        itercolc = SYMB_LCOLNUM(itercblk)+1;
      else
      {
        fprintf(stream, "col : %ld - %ld\n",
                (long)itercola, (long)(itercolc-1));
        itercola = SYMB_FCOLNUM(itercblk);
        itercolc = SYMB_LCOLNUM(itercblk)+1;
      }
    }
  }
  fprintf(stream, "col : %ld - %ld\n",
          (long)itercola, (long)(itercolc-1));

#ifdef CSC_LOG
  fprintf(stdout, "<- z_dump4 \n");
#endif
}

/*
 * Function: z_dump5
 *
 * Writes right-hand-side memeber in a file.
 *
 * Parameters:
 *   datacode - z_SolverMatrix containing right-hand-side member.
 *   stream   - FILE * opened in write mode.
 */
void z_dump5(const z_SolverMatrix *datacode,
           FILE               *stream)
{
  pastix_int_t itercblk;
  pastix_int_t iterupdo;
  pastix_int_t itercolo;
#ifdef MULT_SMX
  pastix_int_t itersmx;
#endif

#ifdef CSC_LOG
  fprintf(stdout, "-> z_dump5 \n");
#endif

#ifdef MULT_SMX
  for (itersmx=0; itersmx<datacode->updovct.sm2xnbr; itersmx++)
  {
#endif
    for (itercblk = 0;
         itercblk < SYMB_CBLKNBR;
         itercblk++)
    {
      iterupdo = datacode->updovct.cblktab[itercblk].sm2xind;

      for (itercolo = SYMB_FCOLNUM(itercblk);
           itercolo < SYMB_LCOLNUM(itercblk)+1;
           itercolo++)
      {
#ifdef TYPE_COMPLEX
        fprintf(stream, "%ld (%.13e,%.13e)\n", (long)itercolo,
                creal(datacode->updovct.sm2xtab[iterupdo]), cimag(datacode->updovct.sm2xtab[iterupdo]));
#else
#ifdef MULT_SMX
        fprintf(stream, "%ld %.13e\n", (long)itercolo,
                datacode->updovct.sm2xtab[itersmx*datacode->updovct.sm2xsze+iterupdo]);
#else
        fprintf(stream, "%ld %.13e\n", (long)itercolo,
                datacode->updovct.sm2xtab[iterupdo]);
#endif
#endif
        iterupdo++;
      }
    }
#ifdef MULT_SMX
  }
#endif
#ifdef CSC_LOG
  fprintf(stdout, "<- z_dump5 \n");
#endif
}

/*
 * Function: z_dump6
 *
 * Prints diagonal blocks in the folowing format :
 * > ** block diag <cblknbr> **
 * > <line1> [<value1> <value2> ... ]
 * > <line2> [...                   ]
 *
 * Prints one file dor L and one for U.
 *
 * Parameters:
 *   datacode - z_SolverMatrix.
 *   streamL  - FILE * into which L diagonal blocs will be writen.
 *   streamU  - FILE * into which U diagonal blocs will be writen.
 */
void z_dump6(const z_SolverMatrix *datacode,
           FILE               *streamL,
           FILE               *streamU)
{
  pastix_int_t itercblk;
  pastix_int_t iterbloc;
  pastix_int_t itercolo;
  pastix_int_t iterline;
  pastix_int_t coefindx;
  pastix_int_t stride;
  pastix_int_t i,j;

#ifdef CSC_LOG
  fprintf(stdout, "-> z_dump6 \n");
#endif

  for (itercblk = 0;
       itercblk < SYMB_CBLKNBR;
       itercblk++)
  {

    /* bloc diag */
    iterbloc = SYMB_BLOKNUM(itercblk);
    coefindx = SOLV_COEFIND(iterbloc);
    stride = SOLV_STRIDE(itercblk);

    fprintf(streamL, " ** block diag %ld **\n", (long)itercblk);
    fprintf(streamU, " ** block diag %ld **\n", (long)itercblk);

    for (iterline = SYMB_FROWNUM(iterbloc);
         iterline < SYMB_LROWNUM(iterbloc)+1;
         iterline++)
    {

      fprintf(streamL, "%ld [  ", (long)iterline);
      fprintf(streamU, "%ld [  ", (long)iterline);

      for (itercolo = SYMB_FCOLNUM(itercblk);
           itercolo < SYMB_LCOLNUM(itercblk)+1;
           itercolo++)
      {
        i=iterline - SYMB_FROWNUM(iterbloc);
        j=itercolo - SYMB_FCOLNUM(itercblk);


#ifdef TYPE_COMPLEX
        fprintf(streamL, "(%.4g,%.4g) ",
                creal(SOLV_COEFTAB(itercblk)[coefindx + (j * stride) + i ]),
                cimag(SOLV_COEFTAB(itercblk)[coefindx + (j * stride) + i ]));
        fprintf(streamU, "(%.4g,%.4g) ",
                creal(SOLV_UCOEFTAB(itercblk)[coefindx + (j * stride) + i ]),
                cimag(SOLV_UCOEFTAB(itercblk)[coefindx + (j * stride) + i ]));
#else
        fprintf(streamL, "%.4g ", SOLV_COEFTAB(itercblk)[coefindx + (j * stride) + i ]);
        fprintf(streamU, "%.4g ", SOLV_UCOEFTAB(itercblk)[coefindx + (j * stride) + i ]);
#endif
      }
      fprintf(streamL, "]\n");
      fprintf(streamU, "]\n");

    }
  }

#ifdef CSC_LOG
  fprintf(stdout, "<- z_dump6 \n");
#endif
}

/*
 * Function: z_dump7
 *
 * Writes a vector in the folowing format :
 * > <line1> <value[line1]>
 * > <line2> <value[line1]>
 *
 * Parameters:
 *   v      - vector to write.
 *   stream - FILE * opened in write mode.
 *   nbr    - Size of the vector v.
 */
void z_dump7(pastix_complex64_t *v,
           FILE  *stream,
           pastix_int_t    colnbr)
{
  pastix_int_t iter;

#ifdef CSC_LOG
  fprintf(stdout, "-> z_dump7 \n");
#endif

  for (iter=0; iter<colnbr; iter++)
  {
    fprintf(stream, "%ld %lf\n", (long) iter, (double)v[iter]);
  }

#ifdef CSC_LOG
  fprintf(stdout, "<- z_dump7 \n");
#endif
}
