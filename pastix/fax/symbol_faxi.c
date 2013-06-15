/*
   File: symbol_faxi_nomerge.c

   Part of a parallel direct block solver.
   This is the incomplete block symbolic
   factorization routine.

   This code is based on the one of
   symbol_faxi.c . In fact, this code
   is exactly the one of version 1.1 of
   fax.

   This file can be included by faxi_graph.c.

   Authors:
     - Francois Pellegrini

   Dates:
     Version 3.0 - from 07 dec 2004 to 15 dec 2004
*/

/*
  Section: The defines and macros.
*/
/*
  Macros:

  SYMBOL_FAXI_INCLUDED - Has to be defined if the code his included
                         in an external file function (see <faxi_graph.c>)
  SYMBOL_FAXI          - Defined if SYMBOL_FAXI_INCLUDED is not...
                         But does not seem to be used...

*/
#ifndef SYMBOL_FAXI_INCLUDED                      /* Can be included by faxi_graph.c */
#define SYMBOL_FAXI

#include "common.h"
#include "symbol.h"
#include "order.h"
#include "fax.h"
#include "symbol_faxi.h"
#endif /* SYMBOL_FAXI_INCLUDED */

/*
  Macro: SYMBOL_FAXI_ITERATOR

  Loop for all adjacent edges, used in <symbolFaxi>.
  Must be defined in including file if SYMBOL_FAXI_INCLUDED is defined.

  Parameters:
    ngbdptr - Neighbour pointer.
    vertnum - Vertex index.
    vertend - Iterator.
*/
#ifndef SYMBOL_FAXI_INCLUDED
#define SYMBOL_FAXI_ITERATOR(ngbdptr, vertnum, vertend)			\
  for (vertend  = ngbfrst (ngbdptr, vertnum);				\
       vertend != baseval - 1;						\
       vertend  = ngbnext (ngbdptr)) {
/*
  Macro: SYMBOL_FAXI_VERTEX_DEGREE

  Computes the number of adjacent edges to a vertex.

  Parameters:
    ngbdptr - Neighbour pointer.
    vertnum - Vertex index.
*/
#define SYMBOL_FAXI_VERTEX_DEGREE(ngbdptr, vertnum)			\
  (ngbdegr ((ngbdptr), (vertnum)))

/*
  Section: Functions
*/

/*
  Function: symbolFaxi

  Incomplete symbolic factorization routine
  with limitation of level-of-fill value.

  This routine computes the block symbolic
  factorization of the given matrix graph
  according to the given vertex ordering.

  The algorithm is similar to the one of
  complete symbolic factorization, except
  that a level of fill value is recorded
  for each block, which serves to compute
  the level of fill for each filled block.

  Parameters:
    symbptr - Symbolic block matrix [based]
    vertnbr - Number of vertices
    edgenbr - Number of edges
    baseval - Base value
    ngbdptr - Neighbor bookkeeping area
    ngbfrst - First neighbor function
    ngbnext - Next neighbor function
    ngbdegr - Vertex degree function (upper bound)
    ordeptr - Matrix ordering
    levfmax - Inclusive maximum level of fill for blocks

  Returns:
    0  - on success.
    !0 - on error.

*/
int
symbolFaxi (SymbolMatrix * const        symbptr,
            const pastix_int_t                   vertnbr,
            const pastix_int_t                   edgenbr,
            const pastix_int_t                   baseval,
            void * const                ngbdptr,
            pastix_int_t                         ngbfrst (void * const, const pastix_int_t),
            pastix_int_t                         ngbnext (void * const),
            pastix_int_t                         ngbdegr (void * const, const pastix_int_t),
            const Order * const         ordeptr,
            const pastix_int_t                   levfmax)
#endif /* SYMBOL_FAXI_INCLUDED */
{
  pastix_int_t                       vertnum;              /* Vertex number of current column                   */
  pastix_int_t                       vertend;              /* Current end vertex number                         */
  const pastix_int_t * restrict      permtax;              /* Based access to direct permutation array          */
  const pastix_int_t * restrict      peritax;              /* Based access to inverse permutation array         */
  const pastix_int_t * restrict      rangtax;              /* Based access to column block range array          */
  pastix_int_t * restrict            ctrbtax;              /* Based access to array of contribution chains      */
  SymbolCblk * restrict     cblktax;              /* Based access to column block array                */
  pastix_int_t                       cblknum;              /* Based number of current column block              */
  pastix_int_t                       cblkctr;              /* Based number of current contributing column block */
  SymbolBlok * restrict     bloktax;              /* Based access to incomplete block array            */
  pastix_int_t                       bloknum;              /* Based number of current first free block slot     */
  pastix_int_t                       blokmax;              /* Maximum number of blocks in array                 */
  SymbolFaxiTlok * restrict tloktab;              /* Beginning of array of temporary blocks            */
  pastix_int_t                       ctrbsum;              /* Number of contributing blocks for column block    */
  pastix_int_t * restrict            sorttab;              /* Beginning of sort area                            */
  pastix_int_t                       sortnbr;              /* Number of vertices in sort area and hash table    */
  pastix_int_t * restrict            hashtab;              /* Hash vertex table                                 */
  pastix_int_t                       hashmsk;              /* Mask for access to hash table                     */
  pastix_int_t                       colend;               /* Column number of vertex neighbor                  */

  permtax = ordeptr->permtab - baseval;           /* Compute array bases */
  peritax = ordeptr->peritab - baseval;
  rangtax = ordeptr->rangtab - baseval;

  blokmax  = ordeptr->cblknbr * (2 + edgenbr / vertnbr) + 2; /* Estimate size of initial block array */

  {     /* Allocate arrays for factoring   */
    pastix_int_t        * ctrbtab = NULL; /* Array for contribution chaining */
    SymbolCblk * cblktab = NULL; /* Column block array              */
    SymbolBlok * bloktab = NULL; /* Incomplete block array          */

    MALLOC_INTERN(ctrbtab, ordeptr->cblknbr,     pastix_int_t);
    MALLOC_INTERN(cblktab, ordeptr->cblknbr + 1, SymbolCblk);
    MALLOC_INTERN(bloktab, blokmax,              SymbolBlok);

    cblktax = cblktab - baseval;                  /* Set based accesses */
    bloktax = bloktab - baseval;
    ctrbtax = ctrbtab - baseval;

    memset (ctrbtab, ~0, ordeptr->cblknbr * sizeof (pastix_int_t)); /* Initialize column block contributions link array */
  }

  bloknum = baseval;
  for (cblknum = baseval; cblknum < baseval + ordeptr->cblknbr; cblknum ++) { /* For all column blocks */
    pastix_int_t                 colnum;                   /* Number of current column [based]                  */
    pastix_int_t                 colmax;                   /* Maximum column index for current column block     */

    {                                             /* Compute offsets and check for array size */
      pastix_int_t                 degrsum;
      pastix_int_t                 hashsiz;
      pastix_int_t                 hashmax;
      pastix_int_t                 ctrbtmp;
      pastix_int_t                 sortoft;                /* Offset of sort array                   */
      pastix_int_t                 tlokoft;                /* Offset of temporary block array        */
      pastix_int_t                 tlndoft;                /* Offset of end of temporary block array */
      pastix_int_t                 tlokmax;

      colnum = rangtax[cblknum];
      colmax = rangtax[cblknum + 1];              /* Get maximum column value */

      cblktax[cblknum].fcolnum = colnum;          /* Set column block data */
      cblktax[cblknum].lcolnum = colmax - 1;
      cblktax[cblknum].bloknum = bloknum;

      degrsum = 1;                                /* One for diagonal block                    */
      for ( ; colnum < colmax; colnum ++)         /* For all columns                           */
        degrsum += SYMBOL_FAXI_VERTEX_DEGREE (ngbdptr, peritax[colnum]); /* Add column degrees */

      for (hashmax = 256; hashmax < degrsum; hashmax *= 2) ; /* Get upper bound on hash table size */
      hashsiz = hashmax << 2;                     /* Fill hash table at 1/4 of capacity            */
      hashmsk = hashsiz - 1;

      for (ctrbsum = 2, ctrbtmp = ctrbtax[cblknum]; /* Follow chain of contributing column blocks */
           ctrbtmp != ~0; ctrbtmp = ctrbtax[ctrbtmp])
        ctrbsum += cblktax[ctrbtmp + 1].bloknum - cblktax[ctrbtmp].bloknum - 2; /* Sum contributing column blocks */
      ctrbsum *= 2;                               /* Each contribution can split a block into 2 other pieces      */

      tlokmax = degrsum + ctrbsum;
      sortoft = tlokmax * sizeof (SymbolBlok);
      if ((hashsiz * (pastix_int_t)sizeof(pastix_int_t)) > sortoft)     /* Compute offset of sort area */
        sortoft = (hashsiz * sizeof (pastix_int_t));
      tlokoft = sortoft + degrsum * sizeof (pastix_int_t); /* Compute offset of temporary block area */
      tlndoft = tlokoft + tlokmax * sizeof (SymbolFaxiTlok); /* Compute end of area         */

      if (((unsigned char *) (bloktax + bloknum) + tlndoft) > /* If not enough room */
          ((unsigned char *) (bloktax + blokmax))) {
        SymbolBlok *        bloktmp;              /* Temporary pointer for array resizing */

        do
          blokmax = blokmax + (blokmax >> 2) + 4; /* Increase block array size by 25% as long as it does not fit */
        while (((unsigned char *) (bloktax + bloknum) + tlndoft) > ((unsigned char *) (bloktax + blokmax)));

        if ((bloktmp = (SymbolBlok *) memRealloc (bloktax + baseval, (blokmax * sizeof (SymbolBlok)))) == NULL) {
          errorPrint ("symbolFaxi: out of memory (2)");
          memFree    (bloktax + baseval);
          memFree    (cblktax + baseval);
          memFree    (ctrbtax + baseval);
          return     (1);
        }
        bloktax = bloktmp - baseval;
      }

      hashtab = (pastix_int_t *)            (bloktax + bloknum);
      sorttab = (pastix_int_t *)            ((unsigned char *) hashtab + sortoft);
      tloktab = (SymbolFaxiTlok *) ((unsigned char *) hashtab + tlokoft);

      memset (hashtab, ~0, hashsiz * sizeof (pastix_int_t)); /* Initialize hash table */
    }

    sortnbr = 0;                                  /* No vertices yet                 */
    for (colnum = rangtax[cblknum]; colnum < colmax; colnum ++) { /* For all columns */
      pastix_int_t                 hashnum;

      vertnum = peritax[colnum];                  /* Get associated vertex       */
      SYMBOL_FAXI_ITERATOR (ngbdptr, vertnum, vertend) /* For all adjacent edges */

        colend = permtax[vertend];                /* Get end column number */

        if (colend < colmax)                      /* If end vertex number in left columns */
          continue;                               /* Skip to next neighbor                */

        for (hashnum = (colend * SYMBOL_FAXI_HASHPRIME) & hashmsk; ; /* Search end column in hash table */
             hashnum = (hashnum + 1) & hashmsk) {
          pastix_int_t *               hashptr;

          hashptr = hashtab + hashnum;            /* Point to hash slot           */
          if (*hashptr == colend)                 /* If end column in hash table  */
            break;                                /* Skip to next end column      */
          if (*hashptr == ~0) {                   /* If slot is empty             */
            *hashptr = colend;                    /* Set column in hash table     */
            sorttab[sortnbr ++] = colend;         /* Add end column to sort array */
            break;
          }
        }
      }                                           /* End of loop on neighbors */
    }                                             /* End of loop on columns   */

    intSort1asc1 (sorttab, sortnbr);              /* Sort neighbor array */

    cblkctr = cblknum;
    if (ctrbtax[cblknum] == ~0) {                 /* If column is not to be updated */
      pastix_int_t                 sortnum;

      bloktax[bloknum].frownum = cblktax[cblknum].fcolnum; /* Build diagonal block */
      bloktax[bloknum].lrownum = cblktax[cblknum].lcolnum;
      bloktax[bloknum].cblknum = cblknum;
      bloktax[bloknum].levfval = 0;               /* Original block */
      bloknum ++;

      for (sortnum = 0; sortnum < sortnbr; ) {    /* For all entries in sorted array */

        colend = sorttab[sortnum];
        if (colend >= rangtax[cblkctr + 1]) {     /* If column block number to be found */
          pastix_int_t                 cblktmm;            /* Median value                       */
          pastix_int_t                 cblktmx;            /* Maximum value                      */

          for (cblkctr ++,                        /* Find new column block by dichotomy */
               cblktmx = ordeptr->cblknbr + baseval;
               cblktmx - cblkctr > 1; ) {
            cblktmm = (cblktmx + cblkctr) >> 1;
            if (rangtax[cblktmm] <= colend)
              cblkctr = cblktmm;
            else
              cblktmx = cblktmm;
          }
        }

        bloktax[bloknum].frownum = colend;        /* Set beginning of new block */
        while ((++ sortnum < sortnbr) &&          /* Scan extent of block       */
               (sorttab[sortnum] - 1 == sorttab[sortnum - 1]) &&
               (sorttab[sortnum] < rangtax[cblkctr + 1])) ;
        bloktax[bloknum].lrownum = sorttab[sortnum - 1]; /* Set end of block */
        bloktax[bloknum].cblknum = cblkctr;
        bloktax[bloknum].levfval = 0;             /* Original block */
        bloknum ++;                               /* One more block */
      }
    }
    else {                                        /* Column will be updated           */
      pastix_int_t                 sortnum;                /* Current index in sort array      */
      pastix_int_t                 tloknum;                /* Current index on temporary block */
      pastix_int_t                 tlokfre;                /* Index of first free block        */

      tloktab[0].nextnum = 1;                     /* Dummy temporary block for insertion */
      tloknum            = 1;

      for (sortnum = 0; sortnum < sortnbr; ) {    /* For all entries in sorted array */

        colend = sorttab[sortnum];
        if (colend >= rangtax[cblkctr + 1]) {     /* If column block number to be found */
          pastix_int_t                 cblktmm;            /* Median value                       */
          pastix_int_t                 cblktmx;            /* Maximum value                      */

          for (cblkctr ++,                        /* Find new column block by dichotomy */
               cblktmx = ordeptr->cblknbr + baseval;
               cblktmx - cblkctr > 1; ) {
            cblktmm = (cblktmx + cblkctr) >> 1;
            if (rangtax[cblktmm] <= colend)
              cblkctr = cblktmm;
            else
              cblktmx = cblktmm;
          }
        }
        tloktab[tloknum].frownum = colend;        /* Set beginning of new block */
        while ((++ sortnum < sortnbr) &&          /* Scan extent of block       */
               (sorttab[sortnum] - 1 == sorttab[sortnum - 1]) &&
               (sorttab[sortnum] < rangtax[cblkctr + 1])) ;
        tloktab[tloknum].lrownum = sorttab[sortnum - 1]; /* Set end of block */
        tloktab[tloknum].cblknum = cblkctr;
        tloktab[tloknum].levfval = 0;             /* Original block */
        tloktab[tloknum].nextnum = tloknum + 1;   /* Chain block    */
        tloknum = tloknum + 1;
      }
      tloktab[tloknum].frownum =                  /* Build trailing block */
      tloktab[tloknum].lrownum = vertnbr + baseval;
      tloktab[tloknum].cblknum = ordeptr->cblknbr + baseval;
      tloktab[tloknum].nextnum = ~0;              /* Set end of chain */

      tlokfre = ++ tloknum;                       /* Build free chain for possible contributing blocks */
      for ( ; tloknum < tlokfre + ctrbsum; tloknum = tloknum + 1)
        tloktab[tloknum].nextnum = tloknum + 1;
      tloktab[tloknum].nextnum = ~0;              /* Set end of free chain */

      for (cblkctr = ctrbtax[cblknum]; cblkctr != ~0; cblkctr = ctrbtax[cblkctr]) { /* Follow chain */
        pastix_int_t                 blokctr;              /* Current index of contributing column block     */
        pastix_int_t                 levffac;              /* Minimum level of fill over facing block(s)     */
        pastix_int_t                 tloklst;              /* Index of previous temporary block              */

        blokctr = cblktax[cblkctr].bloknum + 1;   /* Get index of first extra-diagonal block                 */
        levffac = bloktax[blokctr].levfval;       /* Get level of fill of first extra-diagonal block         */
        for (blokctr ++; (blokctr < cblktax[cblkctr + 1].bloknum) && /* For all blocks facing diagonal block */
             (bloktax[blokctr].cblknum == cblknum); blokctr ++) {
          if (bloktax[blokctr].levfval < levffac) /* If facing block has smaller level of fill */
            levffac = bloktax[blokctr].levfval;   /* Keep smallest level of fill               */
        }

        tloklst = 0;                              /* Previous is diagonal block            */
        tloknum = tloktab[0].nextnum;             /* Current is first extra-diagonal block */

        for ( ; blokctr < cblktax[cblkctr + 1].bloknum; blokctr ++) { /* For all blocks in contributing column block   */
          pastix_int_t                 frownum;            /* Current first and last rows of contributing block(s) to be merged */
          pastix_int_t                 lrownum;
          pastix_int_t                 levfval;            /* Level of fill of block(s) to be created */

          frownum = bloktax[blokctr].frownum;     /* Get extents of block to merge */
          lrownum = bloktax[blokctr].lrownum;
          levfval = 1 + ((bloktax[blokctr].levfval > levffac) ? bloktax[blokctr].levfval : levffac); /* Compute level of fill */

          if (levfval > levfmax)                  /* If above maximum level of fill allowed */
            continue;                             /* Forget about the contribution          */

          while (tloktab[tloknum].lrownum < frownum) { /* Skip unmatched chained blocks */
            tloklst = tloknum;
            tloknum = tloktab[tloknum].nextnum;
          }

          for ( ; ; ) {                           /* Process chained blocks   */
            if (lrownum < tloktab[tloknum].frownum) { /* If block has no mate */
#ifdef FAX_DEBUG
              if (tlokfre == ~0) {
                errorPrint ("symbolFaxi: internal error (1)");
                memFree    (bloktax + baseval);
                memFree    (cblktax + baseval);
                memFree    (ctrbtax + baseval);
                return     (1);
              }
#endif /* FAX_DEBUG */
              tloktab[tlokfre].frownum = frownum; /* Copy block data */
              tloktab[tlokfre].lrownum = lrownum;
              tloktab[tlokfre].cblknum = bloktax[blokctr].cblknum;
              tloktab[tlokfre].levfval = levfval;
              tloktab[tloklst].nextnum = tlokfre; /* Chain new block */
              tloklst                  = tlokfre;
              tlokfre                  = tloktab[tloklst].nextnum;
              tloktab[tloklst].nextnum = tloknum;
              break;                              /* Skip to next contribution block */
            }

            if (tloktab[tloknum].levfval < levfval) { /* If chained block has lower level */
              if (frownum < tloktab[tloknum].frownum) {
#ifdef FAX_DEBUG
                if (tlokfre == ~0) {
                  errorPrint ("symbolFaxi: internal error (2)");
                  memFree    (bloktax + baseval);
                  memFree    (cblktax + baseval);
                  memFree    (ctrbtax + baseval);
                  return     (1);
                }
#endif /* FAX_DEBUG */
                tloktab[tlokfre].frownum = frownum; /* Copy block data */
                tloktab[tlokfre].lrownum = tloktab[tloknum].frownum - 1;
                tloktab[tlokfre].cblknum = tloktab[tloknum].cblknum; /* Same as bloktab[blokctr].cblknum */
                tloktab[tlokfre].levfval = levfval;
                tloktab[tloklst].nextnum = tlokfre; /* Chain new block */
                tloklst                  = tlokfre;
                tlokfre                  = tloktab[tloklst].nextnum;
                tloktab[tloklst].nextnum = tloknum;
              }
              frownum = tloktab[tloknum].lrownum + 1; /* Re-start block after chained block */

              if (frownum > lrownum)              /* If end of contribution block reached */
                break;                            /* Skip to next contribution block      */

              tloklst = tloknum;                  /* Else process next chained block with remains of same contribution block */
              tloknum = tloktab[tloknum].nextnum;
            }
            else {                                /* If contribution block has lower or same level                 */
              if (tloktab[tloknum].lrownum > lrownum) { /* If chained block spans after contribution block         */
                if (tloktab[tloknum].frownum < frownum) { /*  If chained block spans before contribution block too */
#ifdef FAX_DEBUG
                  if (tlokfre == ~0) {
                    errorPrint ("symbolFaxi: internal error (3)");
                    memFree    (bloktax + baseval);
                    memFree    (cblktax + baseval);
                    memFree    (ctrbtax + baseval);
                    return     (1);
                  }
#endif /* FAX_DEBUG */
                  tloktab[tlokfre].frownum = lrownum + 1; /* Build new block for second remaining part and resize */
                  tloktab[tlokfre].lrownum = tloktab[tloknum].lrownum;
                  tloktab[tlokfre].cblknum = tloktab[tloknum].cblknum;
                  tloktab[tlokfre].levfval = tloktab[tloknum].levfval;
                  tloktab[tloknum].lrownum = frownum - 1;
                  tloklst                  = tloknum;
                  tloknum                  = tlokfre;
                  tlokfre                  = tloktab[tloknum].nextnum;
                  tloktab[tloknum].nextnum = tloktab[tloklst].nextnum;
                  tloktab[tloklst].nextnum = tloknum;
                }
                else                              /* Chained block does not span before contribution block */
                  tloktab[tloknum].frownum = lrownum + 1; /* Resize chained block to avoid collision       */
              }
              else {                              /* Chained block does not span after contribution block     */
                if (tloktab[tloknum].frownum < frownum) { /* If chained block spans before contribution block */
                  tloktab[tloknum].lrownum = frownum - 1; /* Resize chained block to avoid collision          */
                  tloklst                  = tloknum; /* Skip to next block as this one is now before         */
                  tloknum                  = tloktab[tloklst].nextnum;
                }
                else {                            /* Chained block completely shadowed */
                  tloktab[tloklst].nextnum = tloktab[tloknum].nextnum; /* Remove it    */
                  tloktab[tloknum].nextnum = tlokfre;
                  tlokfre                  = tloknum;
                  tloknum                  = tloktab[tloklst].nextnum;
                }
              }
            }
          }
        }
      }

      bloktax[bloknum].frownum = cblktax[cblknum].fcolnum; /* Build diagonal block */
      bloktax[bloknum].lrownum = cblktax[cblknum].lcolnum;
      bloktax[bloknum].cblknum = cblknum;
      bloktax[bloknum].levfval = 0;               /* Original block */
      bloknum ++;

      for (tloknum = tloktab[0].nextnum;          /* For all chained blocks         */
           tloktab[tloknum].nextnum != ~0;        /* Until trailer block is reached */
           tloknum = tloktab[tloknum].nextnum) {  /* Copy block data to block array */
        if ((tloktab[tloknum].cblknum == bloktax[bloknum - 1].cblknum) &&
            (tloktab[tloknum].levfval == bloktax[bloknum - 1].levfval) &&
            (tloktab[tloknum].frownum == bloktax[bloknum - 1].lrownum + 1))
          bloktax[bloknum - 1].lrownum = tloktab[tloknum].lrownum; /* Merge blocks with same level of fill */
        else {
          bloktax[bloknum].frownum = tloktab[tloknum].frownum;
          bloktax[bloknum].lrownum = tloktab[tloknum].lrownum;
          bloktax[bloknum].cblknum = tloktab[tloknum].cblknum;
          bloktax[bloknum].levfval = tloktab[tloknum].levfval;
          bloknum ++;
        }
      }
    }
    if ((bloknum - cblktax[cblknum].bloknum) > 2) { /* If more than one extra-diagonal blocks exist                 */
      ctrbtax[cblknum] = ctrbtax[bloktax[cblktax[cblknum].bloknum + 1].cblknum]; /* Link contributing column blocks */
      ctrbtax[bloktax[cblktax[cblknum].bloknum + 1].cblknum] = cblknum;
    }
  }
  cblktax[cblknum].fcolnum =                      /* Set last column block data */
  cblktax[cblknum].lcolnum = vertnbr + baseval;
  cblktax[cblknum].bloknum = bloknum;

  memFree (ctrbtax + baseval);                    /* Free contribution link array */

  symbptr->baseval = baseval;                     /* Fill in matrix fields */
  symbptr->cblknbr = ordeptr->cblknbr;
  symbptr->bloknbr = bloknum - baseval;
  symbptr->cblktab = cblktax + baseval;
  symbptr->bloktab = (SymbolBlok *) memRealloc (bloktax + baseval, (bloknum - baseval) * sizeof (SymbolBlok)); /* Set array to its exact size */
  symbptr->nodenbr = vertnbr;

#ifdef FAX_DEBUG
  if (symbolCheck (symbptr) != 0) {
    errorPrint ("symbolFaxi: internal error (4)");
    symbolExit (symbptr);
    return     (1);
  }
#endif /* FAX_DEBUG */

  return PASTIX_SUCCESS;
}
