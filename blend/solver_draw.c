/* Copyright INRIA 2004
**
** This file is part of the Scotch distribution.
**
** The Scotch distribution is libre/free software; you can
** redistribute it and/or modify it under the terms of the
** GNU Lesser General Public License as published by the
** Free Software Foundation; either version 2.1 of the
** License, or (at your option) any later version.
**
** The Scotch distribution is distributed in the hope that
** it will be useful, but WITHOUT ANY WARRANTY; without even
** the implied warranty of MERCHANTABILITY or FITNESS FOR A
** PARTICULAR PURPOSE. See the GNU Lesser General Public
** License for more details.
**
** You should have received a copy of the GNU Lesser General
** Public License along with the Scotch distribution; if not,
** write to the Free Software Foundation, Inc.,
** 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
**
** $Id: symbol_draw.c 285 2005-03-10 10:25:31Z pelegrin $
*/
/************************************************************/
/**                                                        **/
/**   NAME       : symbol_draw.c                           **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                This module draws symbolic matrices in  **/
/**                PostScript (tm) format .                **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 29 sep 1998     **/
/**                                 to     29 sep 1998     **/
/**                # Version 1.0  : from : 26 jun 2002     **/
/**                                 to     26 jun 2002     **/
/**                # Version 1.3  : from : 10 apr 2003     **/
/**                                 to     10 jun 2003     **/
/**                # Version 3.0  : from : 29 feb 2004     **/
/**                                 to     29 feb 2004     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/
#include <time.h>
#include "common.h"
#include "symbol.h"
#include "solver.h"

/*+ Generic PostScript (tm) output definitions. +*/

#define SYMBOL_PSDPI      72    /*+ PostScript dots-per-inch            +*/
#define SYMBOL_PSPICTSIZE 6.6   /*+ PostScript picture size (in inches) +*/

int
solverDrawFunc (
const SolverMatrix * const  solvptr,
FILE * const                stream,
int                         verbose )
{
  pastix_int_t cblknum;                    /* Number of current column block */
  time_t       picttime;                   /* Creation time                  */
  double       pictsize;                   /* Number of distinct coordinates */
  int          o;

  time (&picttime);                               /* Get current time */
  pictsize = (double) (solvptr->nodenbr + 1);     /* Get matrix size  */

  fprintf (stream, "%%!PS-Adobe-2.0 EPSF-2.0\n"); /* Write header */
  fprintf (stream, "%%%%Title: symbolmatrix (%ld,%ld,%ld)\n",
           (long) solvptr->cblknbr, (long) solvptr->bloknbr, (long)solvptr->nodenbr);
  fprintf (stream, "%%%%Creator: symbolDraw (LaBRI, Universite Bordeaux I)\n");
  fprintf (stream, "%%%%CreationDate: %s", ctime (&picttime));
  fprintf (stream, "%%%%BoundingBox: 0 0 %ld %ld\n",
           (long) (SYMBOL_PSPICTSIZE * SYMBOL_PSDPI),
           (long) (SYMBOL_PSPICTSIZE * SYMBOL_PSDPI));
  fprintf (stream, "%%%%Pages: 0\n");
  fprintf (stream, "%%%%EndComments\n");          /* Write shortcuts */
  fprintf (stream, "/c { 4 2 roll pop pop newpath 2 copy 2 copy moveto dup lineto dup lineto closepath fill } bind def\n");
  fprintf (stream, "/d { 4 2 roll pop pop newpath 2 copy 2 copy moveto dup lineto dup lineto closepath } bind def\n");
  fprintf (stream, "/b { 4 copy 2 index exch moveto lineto dup 3 index lineto exch lineto closepath fill pop } bind def\n");
  fprintf (stream, "/a { 4 copy 2 index exch moveto lineto dup 3 index lineto exch lineto closepath pop } bind def\n");
  fprintf (stream, "/r { setrgbcolor } bind def\n");
  fprintf (stream, "/g { setgray } bind def\n");

  fprintf (stream, "0 setlinecap\n");             /* Use miter caps       */
  fprintf (stream, "%f dup scale\n",              /* Print scaling factor */
           (double) SYMBOL_PSDPI * SYMBOL_PSPICTSIZE / pictsize);
  fprintf (stream, "/Times-Roman 70 selectfont\n"); /* activate text in eps file */
  fprintf (stream, "[ 1 0 0 -1 0 %d ] concat\n",  /* Reverse Y coordinate */
           (int) (solvptr->nodenbr + 1));

  fprintf (stream, "0 0\n");                      /* Output fake column block */
  for (cblknum = 0; cblknum < solvptr->cblknbr; cblknum ++) {
    float               coloval[3];               /* Color of diagonal block and previous color */

    coloval[0] = 0.5;
    coloval[1] = 0.5;
    coloval[2] = 0.5;
    if ((coloval[0] == coloval[1]) &&
        (coloval[1] == coloval[2]))
      fprintf (stream, "%.2g g ",
               (float) coloval[0]);
    else
      fprintf (stream, "%.2g %.2g %.2g r \n",
               (float) coloval[0], (float) coloval[1], (float) coloval[2]);

    SolverCblk *cblk   = &solvptr->cblktab[cblknum];

    fprintf (stream, "%ld\t%ld\tc\n",             /* Begin new column block */
             (long) (cblk->fcolnum - solvptr->baseval),
             (long) (cblk->lcolnum - solvptr->baseval + 1));

    pastix_int_t ncols = cblk_colnbr( cblk );
    SolverBlok *blok   = cblk[0].fblokptr+1;
    SolverBlok *lblok  = cblk[1].fblokptr;

    for (; blok<lblok; blok++)
    {
      float               colbval[3];             /* Color of off-diagonal block */
      coloval[0] = colbval[0];                /* Save new color data */
      coloval[1] = colbval[1];
      coloval[2] = colbval[2];


      if ( cblk->cblktype & CBLK_DENSE ) {
          fprintf (stream, "%.2g %.2g %.2g r \n",
                   0.5, 0.5, 0.5);
      }
      else{
          pastix_int_t nrows = blok_rownbr( blok );

          pastix_int_t conso_dense = 2*nrows*ncols;
          pastix_int_t conso_LR    = 0;
          if (blok->LRblock[0].rk != -1){
              conso_LR += (((nrows+ncols) * blok->LRblock[0].rk));
          }
          else{
              conso_LR += nrows*ncols;
          }
          if (blok->LRblock[1].rk != -1){
              conso_LR += (((nrows+ncols) * blok->LRblock[1].rk));
          }
          else{
              conso_LR += nrows*ncols;
          }

          double gain = 1.0 * conso_dense / conso_LR;
          /* printf("Conso LR %ld Dense %ld Gain %f Blok %p\n", conso_LR, conso_dense, gain, blok); */

          /* There is no compression */
          if (gain == 1.){
              fprintf(stream, "%.2g %.2g %.2g r \n",
                      0., 0., 0.);
          }
          /* Small amount of compression: red */
          else if (gain < 5.){
              fprintf(stream, "%.2g %.2g %.2g r \n",
                      gain / 5., 0., 0.);
          }
          /* Huge amount of compression */
          else{
              float color = 0.5 + (gain-5) / 10.;
              if (color > 1)
                  color = 1.;
              fprintf(stream, "%.2g %.2g %.2g r \n",
                       0., color, 0.);
          }
      }

      fprintf (stream, "%ld\t%ld\tb\n",         /* Write block in column block */
               (long) (blok->frownum - solvptr->baseval),
               (long) (blok->lrownum - solvptr->baseval + 1));
    }
  }

  /* Plot numbers */
  if (verbose > 4){
    int nb_bloks = 0;
    int nb_cblks = 0;
    FILE *fd1 = fopen( "contribblok.txt", "r" );
    FILE *fd2 = fopen( "contribcblk.txt", "r" );
    int original_cblk = 1;
    double color = 0.2;

    fprintf (stream, "0 0\n");                      /* Output fake column block */
    for (cblknum = 0; cblknum < solvptr->cblknbr; cblknum ++) {
      SolverCblk *cblk   = &solvptr->cblktab[cblknum];
      int unused, nb_contrib;
      fscanf(fd2, "%d %d %d\n", &unused, &nb_contrib, &original_cblk);

      fprintf (stream, "%.2g g %ld\t%ld\tc\n",             /* Begin new column block */
               color,
               (long) (cblk->fcolnum - solvptr->baseval),
               (long) (cblk->lcolnum - solvptr->baseval + 1));
      if ( ! (cblk->cblktype & CBLK_DENSE) ) {
          fprintf (stream, "%ld\t%ld\t4 copy 3 index exch moveto [ 1 0 0 -1 0 0 ] concat 0.0 0.0 0.0 setrgbcolor (%d) show [ 1 0 0 -1 0 0 ] concat pop\n",
                   (long) (cblk->fcolnum - solvptr->baseval),
                   (long) (cblk->lcolnum - solvptr->baseval + 1),
                   nb_contrib);
      }

      SolverBlok *blok   = cblk[0].fblokptr+1;
      SolverBlok *lblok  = cblk[1].fblokptr;

      for (; blok<lblok; blok++)
      {
        int unused, nb_contrib;
        fscanf(fd1, "%d %d\n", &unused, &nb_contrib);
        fprintf (stream, "%ld\t%ld\ta\n",         /* Write block in column block */
                 (long) (blok->frownum - solvptr->baseval),
                 (long) (blok->lrownum - solvptr->baseval + 1));
        if ( ! (cblk->cblktype & CBLK_DENSE) ) {
            fprintf (stream, "%ld\t%ld\t4 copy 3 index exch moveto [ 1 0 0 -1 0 0 ] concat 1.0 1.0 1.0 setrgbcolor (%d) show [ 1 0 0 -1 0 0 ] concat pop\n",
                     (long) (blok->frownum - solvptr->baseval),
                     (long) (blok->lrownum - solvptr->baseval + 1),
                     nb_contrib);
        }
        nb_bloks++;
      }

      if (original_cblk == 0){
          if (color < 0.3)
              color = 0.8;
          else
              color = 0.2;
      }
      nb_cblks++;
    }
    fclose(fd1);
    fclose(fd2);
  }

  fprintf (stream, "pop pop\n");                  /* Purge last column block indices */
  o = fprintf (stream, "showpage\n");   /* Restore context                 */


  return ((o != EOF) ? 0 : 1);
}

/*+ This routine writes to the given stream
*** a PostScript (tm) picture of the symbolic
*** block matrix, with diagonal blocks in
*** black and off-diagonal blocks in dark gray.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
solverDraw (
const SolverMatrix * const  solvptr,
FILE * const                stream,
int verbose )
{
    return (solverDrawFunc (solvptr, stream, verbose));
}
