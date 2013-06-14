#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "common.h"
#include "symbol.h"
#include "symbolrand.h"

PASTIX_INT hazard(PASTIX_INT a, PASTIX_INT b)
{
  float r;
  PASTIX_INT h;

#ifdef DEBUG_BLEND
  ASSERT(b>=a,MOD_BLEND);
#endif

  r = ((float)rand())/ RAND_MAX;
#ifdef DEBUG_M
  ASSERT(r>= 0,MOD_BLEND);
  ASSERT(r<=1,MOD_BLEND);
#endif

  h = (PASTIX_INT) (a + floor( r*(b-a)+0.5) );

#ifdef DEBUG_BLEND
  ASSERT(h >= a,MOD_BLEND);
  ASSERT(h<=b,MOD_BLEND);
#endif
  return h;
}


void symbolRand(SymbolMatrix *symbmtx, PASTIX_INT h1, PASTIX_INT h2)
{
  
  PASTIX_INT i, j;
  PASTIX_INT bloknbr;
  PASTIX_INT odb;
  PASTIX_INT blok2erase;
  PASTIX_INT bloknum;
  PASTIX_INT erasednbr;
  PASTIX_INT rabotednbr;


  SymbolCblk *cblktab = NULL;
  SymbolBlok *bloktab = NULL;

  bloknbr = symbmtx->bloknbr;
  
  MALLOC_INTERN(cblktab, symbmtx->cblknbr+1, SymbolCblk);
  memcpy(cblktab, symbmtx->cblktab, sizeof(SymbolCblk) * (symbmtx->cblknbr+1));

  MALLOC_INTERN(bloktab, bloknbr, SymbolBlok);

  bloknum = 0;
  erasednbr = 0;
  rabotednbr = 0;
  for(i=0;i<symbmtx->cblknbr;i++)
    {
      odb = symbmtx->cblktab[i+1].bloknum - symbmtx->cblktab[i].bloknum - 1;
      if(odb>2)
	{
	  blok2erase = hazard(0, odb-1);
#ifdef DEBUG_BLEND
	  ASSERT(blok2erase >= 0,MOD_BLEND);
	  ASSERT(blok2erase < odb,MOD_BLEND);
#endif
	}
      else
	blok2erase = 0;


      cblktab[i].bloknum = bloknum;
      memcpy(bloktab + bloknum, symbmtx->bloktab + symbmtx->cblktab[i].bloknum, sizeof(SymbolBlok));
#ifdef DEBUG_BLEND
      ASSERT(bloktab[bloknum].frownum == cblktab[i].fcolnum,MOD_BLEND);
      ASSERT(bloktab[bloknum].lrownum == cblktab[i].lcolnum,MOD_BLEND);
#endif
      bloknum++;

      for(j=symbmtx->cblktab[i].bloknum+1;j<symbmtx->cblktab[i+1].bloknum;j++)
	{
	  if(symbmtx->bloktab[j].levfval != 0) /** DO NOT ERASE A BLOCK FROM A **/
	    if(blok2erase > 0 && hazard(0,h1) == 1 && j>symbmtx->cblktab[i].bloknum+1) /** Never drop the FIRST ODB !!!! (elimination tree not good)**/
	      {
		blok2erase--;
		bloknbr--;
		erasednbr++;
		continue;
	      }
	  memcpy(bloktab + bloknum, symbmtx->bloktab + j, sizeof(SymbolBlok));

	  /** Reduce dimension of this block **/
	  if(symbmtx->bloktab[j].levfval != 0)
	    if(hazard(0,h2) == 1)
	      {
		bloktab[bloknum].frownum = hazard(symbmtx->bloktab[j].frownum, symbmtx->bloktab[j].lrownum);
		bloktab[bloknum].lrownum = hazard(bloktab[bloknum].frownum, symbmtx->bloktab[j].lrownum);
		rabotednbr++;
#ifdef DEBUG_BLEND
		ASSERT(bloktab[bloknum].lrownum >=  bloktab[bloknum].frownum,MOD_BLEND);
#endif
		
	      }


	  bloknum++;
	}
    }
  
  cblktab[i].bloknum = bloknum;
  fprintf(stdout, "Initial CBLKNBR %ld BLOCKNBR %ld \n", (long)symbmtx->cblknbr, (long)symbmtx->bloknbr);
  fprintf(stdout, "SYMBOL_RABOT (%ld %ld): %ld blocks erased %ld blocks raboted \n", (long)h1, (long)h2, (long)erasednbr, (long)rabotednbr);

#ifdef DEBUG_BLEND
  ASSERT(bloknum == bloknbr,MOD_BLEND);
#endif
  
  symbmtx->bloknbr = bloknbr;
  
  memFree_null(symbmtx->cblktab);
  symbmtx->cblktab = cblktab;
  
  memFree_null(symbmtx->bloktab);
  symbmtx->bloktab = bloktab;
  
#ifdef DEBUG_BLEND
  for(i=0;i<symbmtx->cblknbr-1;i++)
    ASSERT(symbmtx->cblktab[i+1].bloknum - symbmtx->cblktab[i].bloknum > 1,MOD_BLEND);
  for(i=0;i<symbmtx->bloknbr;i++)
    {
      PASTIX_INT cblknum;
      cblknum = symbmtx->bloktab[i].cblknum;
      ASSERT(symbmtx->cblktab[cblknum].fcolnum <=  symbmtx->bloktab[i].frownum,MOD_BLEND);
      ASSERT(symbmtx->cblktab[cblknum].lcolnum >=  symbmtx->bloktab[i].lrownum,MOD_BLEND);
    }
#endif


  
}
