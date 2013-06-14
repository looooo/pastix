#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "common.h"
#include "elimin.h"


pastix_int_t egraphInit(EliminGraph *egraph)
{
  egraph->baseval = 0;
  egraph->vertnbr = 0;
  egraph->verttab = NULL;
  egraph->inbltab = NULL;
  egraph->ownetab = NULL;
  return 1;
}

void egraphExit(EliminGraph *egraph)
{
  memFree_null(egraph->verttab);
  memFree_null(egraph->inbltab);
  memFree_null(egraph->ownetab);
  memFree_null(egraph);
}


pastix_int_t treeInit(EliminTree *etree)
{
  etree->baseval = 0;
  etree->nodenbr = 0;
  etree->nodetab = NULL;
  etree->sonstab = NULL;
  return 1;
}

void treeExit(EliminTree *etree)
{
    memFree_null(etree->nodetab);
    memFree_null(etree->sonstab);
    memFree_null(etree);
}


void treePlot(EliminTree *etree, FILE *out)
{
  pastix_int_t i;

  fprintf(out,
	  "digraph G {\n"
	  "\tcolor=white\n"
	  "rankdir=BT;\n");

  for (i=0;  i < etree->nodenbr; i++)
    {
      if ((etree->nodetab[i]).fathnum == -1)
	continue;
      fprintf(out, "\t\"%ld\"->\"%ld\"\n", (long)i, (long)((etree->nodetab[i]).fathnum));
    }

  fprintf(out, "}\n");
}
	    
