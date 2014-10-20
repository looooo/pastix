#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <assert.h>

#include "common.h"
#include "ftgt.h"
#include "symbol.h"
#include "queue.h"
#include "bulles.h"
#include "updown.h"
#include "solver.h"
#include "fanboth2.h"

#undef DEBUG_BLEND

/** reduction is the percentage of the total AUB size that we do not want to violate **/
pastix_int_t Malt2(SolverMatrix *solvmtx, double reduction)
{
  const pastix_int_t fanFixSize = sizeof(FanInTarget);
  pastix_int_t p, procdst;
  pastix_int_t i, j, cursor=0;
  pastix_int_t accessnbr=0;
  pastix_int_t *ftgtsizetab   = NULL;
  pastix_int_t *ftgtaccesstab = NULL;
  pastix_int_t *extraftgtcurtab   = NULL;
  pastix_int_t *extraftgtfirtab   = NULL;
  pastix_int_t extracursor    = 0;
  pastix_int_t totalftgtsize = 0;
  pastix_int_t *ftgtsendnbr  = NULL;
  pastix_int_t *ftgtcursor   = NULL;
  pastix_int_t **ftgtsendtab = NULL;
  pastix_int_t allocmem; /* the allocated memory */
  pastix_int_t maxmem; /* the max memory       */
  pastix_int_t maxalloc = 0;
  pastix_int_t overmem = 0;  /* the excedent memmory that can't be overlap by fanboth */
  pastix_int_t extraftgtnbr = 0; /* the number of extra fan in target */
  ExtraFtgt *extraftgttab = NULL;
  FanInTarget *newftgttab = NULL;
  pastix_int_t *newnumtab          = NULL;
#ifdef DEBUG_BLEND
  pastix_int_t debugaccess = 0;
  pastix_int_t *debugnew2old = NULL;
  pastix_int_t *debugtab = NULL;
#endif

  Queue partialFtgtInd; /** The index of the partial ftgt in the ftgttab **/
  Queue partialFtgtExtraInd;  /** The index from which ftgt is issue the partial ftgt **/
  Queue taskQueue; /* the task ordered by priority */
  Queue latestQueue; /* the ftgt ind ordered by the latest updatest order **/
  Queue toSendQueue; /* the ftgt fully aggregated but not sent ordered by decraesing size */ 
  pastix_int_t *senttab =  NULL; /** The table to know is an ftgt is already sent **/

  
  
  MALLOC_INTERN(ftgtsendnbr, solvmtx->clustnbr, pastix_int_t);
  MALLOC_INTERN(ftgtcursor,  solvmtx->clustnbr, pastix_int_t);
  MALLOC_INTERN(ftgtsendtab, solvmtx->clustnbr, pastix_int_t *);
#ifdef DEBUG_BLEND
  MALLOC_INTERN(debugtab,    solvmtx->ftgtnbr,  pastix_int_t);
#endif
  MALLOC_INTERN(senttab,     solvmtx->ftgtnbr,  pastix_int_t); 

  bzero(senttab, sizeof(pastix_int_t)*solvmtx->ftgtnbr);
  queueInit(&partialFtgtInd, solvmtx->ftgtnbr/2 + 1);
  queueInit(&partialFtgtExtraInd , solvmtx->ftgtnbr/2 + 1);
  queueInit(&taskQueue , solvmtx->tasknbr);
  queueInit(&latestQueue , solvmtx->ftgtnbr/2 + 1);
  queueInit(&toSendQueue , solvmtx->ftgtnbr/2 + 1);

  

  /** Compute the size in octets of fan in target **/
  /** And compute the send queues                 **/
  /** Be carefull fan in target are already expanded **/
  MALLOC_INTERN(ftgtsizetab, solvmtx->ftgtnbr, pastix_int_t);
  bzero(ftgtsendnbr,  solvmtx->clustnbr*sizeof(pastix_int_t));
 
  for(i=0;i<solvmtx->ftgtnbr;i++)
    {
      ftgtsizetab[i] = fanFixSize 
	+ (solvmtx->ftgttab[i].infotab[FTGT_LCOLNUM] -  solvmtx->ftgttab[i].infotab[FTGT_FCOLNUM] + 1)*
	  (solvmtx->ftgttab[i].infotab[FTGT_LROWNUM] -  solvmtx->ftgttab[i].infotab[FTGT_FROWNUM] + 1)*sizeof(double);
      totalftgtsize += ftgtsizetab[i];
      ftgtsendnbr[solvmtx->ftgttab[i].infotab[FTGT_PROCDST]]++;
    }
  for(p=0;p<solvmtx->clustnbr;p++)
    if(ftgtsendnbr[p] > 0)
      {
	MALLOC_INTERN(ftgtsendtab[p], ftgtsendnbr[p], pastix_int_t);
      }
  /** Compute the maxmimum memory allowed to AUB allocations **/
  maxmem = (pastix_int_t)(((float)reduction*totalftgtsize)/100);
  fprintf(stdout, "Total AUB allocations %ld Coef Allocation %ld pourcentage of AUB/Coef %g \n", 
	  (long)totalftgtsize, (long)(solvmtx->coefnbr*sizeof(double)), ((float)totalftgtsize)/(solvmtx->coefnbr*sizeof(double))*100 );
  
  /** Count the number of access to fan in target **/
  for(i=0;i<solvmtx->tasknbr;i++)
    {
      if( (solvmtx->tasktab[i].taskid == E2) && (solvmtx->indtab[solvmtx->tasktab[i].indnum] >= 0) )  
	accessnbr++;
      if(solvmtx->tasktab[i].taskid == COMP_1D)
	{
	  pastix_int_t odb_nbr;
	  odb_nbr = (solvmtx->cblktab[solvmtx->tasktab[i].cblknum+1].fblokptr - solvmtx->bloktab) - solvmtx->tasktab[i].bloknum -1;
	  for(j=solvmtx->tasktab[i].indnum;j<solvmtx->tasktab[i].indnum + (odb_nbr*(odb_nbr+1))/2;j++)
	    if(solvmtx->indtab[j] >= 0) 
	      accessnbr++;
	}
    }
#ifdef DEBUG_BLEND
  {
    pastix_int_t debugaccess = 0;
    for(i=0;i<solvmtx->ftgtnbr;i++)
      debugaccess += solvmtx->ftgttab[i].infotab[FTGT_CTRBNBR];

    /*if(debugaccess != accessnbr)*/
    /*fprintf(stderr, "Debugaccess %ld accesnbr %ld \n", (long)debugaccess, (long)accessnbr);*/
    ASSERT(debugaccess == accessnbr,MOD_BLEND);
  }
#endif

  /** Fill the ftgt access tab **/
  MALLOC_INTERN(ftgtaccesstab, accessnbr+1, pastix_int_t);
  queueClear(&taskQueue);
  for(i=0;i<solvmtx->tasknbr;i++)
    queueAdd(&taskQueue, i, (double)(solvmtx->tasktab[i].prionum));
  
  while(queueSize(&taskQueue)>0)
    {
      i = queueGet(&taskQueue);
      if( (solvmtx->tasktab[i].taskid == E2) && (solvmtx->indtab[solvmtx->tasktab[i].indnum] >= 0) )
	{
	  ftgtaccesstab[cursor] = solvmtx->indtab[solvmtx->tasktab[i].indnum];
	  cursor++;
	}
      if(solvmtx->tasktab[i].taskid == COMP_1D)
	{
	  pastix_int_t odb_nbr;
	  odb_nbr = (solvmtx->cblktab[solvmtx->tasktab[i].cblknum+1].fblokptr - solvmtx->bloktab)- solvmtx->tasktab[i].bloknum -1;
	  for(j=solvmtx->tasktab[i].indnum;j<solvmtx->tasktab[i].indnum + (odb_nbr*(odb_nbr+1))/2;j++)
	    if(solvmtx->indtab[j] >= 0) 
	      {
		ftgtaccesstab[cursor] = solvmtx->indtab[j];
		cursor++;
	      }
	}
    }
#ifdef DEBUG_BLEND
  ASSERT(cursor==accessnbr,MOD_BLEND);
#endif
  
  /** Fill the ftgt send queues **/
  bzero(ftgtsendnbr,  solvmtx->clustnbr*sizeof(pastix_int_t));
  for(i=0;i<solvmtx->ftgtnbr;i++)
    {
      ftgtsendtab[solvmtx->ftgttab[i].infotab[FTGT_PROCDST]][ftgtsendnbr[solvmtx->ftgttab[i].infotab[FTGT_PROCDST]]] = i;
      ftgtsendnbr[solvmtx->ftgttab[i].infotab[FTGT_PROCDST]]++;
    }

  /** Walk through the ftgt **/
  MALLOC_INTERN(extraftgttab, accessnbr, ExtraFtgt);

  MALLOC_INTERN(extraftgtcurtab, solvmtx->ftgtnbr, pastix_int_t);
  MALLOC_INTERN(extraftgtfirtab, solvmtx->ftgtnbr, pastix_int_t);
  memset(extraftgtcurtab, -1, solvmtx->ftgtnbr*sizeof(pastix_int_t));
  memset(extraftgtfirtab, -1, solvmtx->ftgtnbr*sizeof(pastix_int_t));

  bzero(ftgtcursor, solvmtx->clustnbr*sizeof(pastix_int_t));
  allocmem = 0;
  extraftgtnbr = 0;
  cursor = 0;
  for(i=0;i<accessnbr;i++)
    {
      /*fprintf(stdout, "%ld     %ld \n", (long)i, (long)allocmem);*/
      if(solvmtx->ftgttab[ftgtaccesstab[i]].infotab[FTGT_CTRBCNT] == solvmtx->ftgttab[ftgtaccesstab[i]].infotab[FTGT_CTRBNBR])
	{
	  pastix_int_t nextaccess;
	  allocmem += ftgtsizetab[ftgtaccesstab[i]];


	  /*extraftgtcurtab[ftgtaccesstab[i]]= -1;*/

	  while( (allocmem > maxmem) && (queueSize(&latestQueue)>0) )
	    {
	      cursor = getFtgtInd2(solvmtx, senttab, &toSendQueue, &latestQueue);
	      /*printf("TOTO access %ld cursor %ld %allocmem %ld maxmem %ld \n", (long)i, (long)cursor, (long)allocmem, (long)maxmem);*/
	      if(cursor>=0)
		{
		  procdst = solvmtx->ftgttab[cursor].infotab[FTGT_PROCDST];
		  if(solvmtx->ftgttab[cursor].infotab[FTGT_CTRBCNT] == 0)
		    {
		      solvmtx->ftgttab[cursor].infotab[FTGT_PRIONUM] =   solvmtx->ftgttab[ftgtsendtab[procdst][ftgtcursor[procdst]]].infotab[FTGT_PRIONUM];
		    }
		  else
		      {
		  if(extraftgtfirtab[cursor] < 0)
		    extraftgtfirtab[cursor] = extracursor;
		  extraftgtnbr++;
#ifdef DEBUG_BLEND
		  /*ASSERT(solvmtx->ftgttab[cursor].infotab[FTGT_CTRBCNT] > 0,MOD_BLEND);*/
#endif
		  allocmem -= ftgtsizetab[cursor];
		 

		  queueAdd2(&partialFtgtInd, ftgtsendtab[procdst][ftgtcursor[procdst]], 
			    (double)(ftgtsendtab[procdst][ftgtcursor[procdst]]), i);
		  queueAdd2(&partialFtgtExtraInd,  extracursor, 
			    (double)(ftgtsendtab[procdst][ftgtcursor[procdst]]), i);
		  extraftgttab[extracursor].ctrbnbr = solvmtx->ftgttab[cursor].infotab[FTGT_CTRBNBR] 
				- solvmtx->ftgttab[cursor].infotab[FTGT_CTRBCNT];
		  extraftgttab[extracursor].ctrbcnt = extraftgttab[extracursor].ctrbnbr;
		  extraftgttab[extracursor].prionum = solvmtx->ftgttab[ftgtsendtab[procdst][ftgtcursor[procdst]]].infotab[FTGT_PRIONUM];
		  extraftgttab[extracursor].indnum  = ftgtsendtab[procdst][ftgtcursor[procdst]];
		  extraftgttab[extracursor].ftgtnum = cursor;
		  extraftgttab[extracursor].next    = -1;
		  if(extraftgtcurtab[cursor] >= 0)
		    extraftgttab[extraftgtcurtab[cursor]].next = extracursor;
		  extraftgtcurtab[cursor] = extracursor;
#ifdef DEBUG_BLEND
		  ASSERT(ftgtsendtab[procdst][ftgtcursor[procdst]] <= cursor,MOD_BLEND);
		  ASSERT(solvmtx->ftgttab[ftgtsendtab[procdst][ftgtcursor[procdst]]].infotab[FTGT_CTRBCNT] > 0,MOD_BLEND);
#endif
		  /*fprintf(stdout, "Add Partial Ftgt %ld at %ld ctrbnbr %ld issue from %ld  TO task %ld on proc %ld \n", (long)extracursor, (long)ftgtsendtab[procdst][ftgtcursor[procdst]], (long)extraftgttab[extracursor].ctrbnbr, (long)cursor, (long)solvmtx->ftgttab[cursor].infotab[FTGT_TASKDST], (long)procdst);*/
		  solvmtx->ftgttab[cursor].infotab[FTGT_CTRBNBR] = solvmtx->ftgttab[cursor].infotab[FTGT_CTRBCNT];
		  extracursor++;
		      }
		}
	      else
		{
		  if(allocmem - maxmem > overmem)
		    overmem = allocmem - maxmem;
		}
	    }
	  /** Put the new allocated Ftgt in the queue **/
	  nextaccess = getFtgtNextAccess(i, accessnbr, ftgtaccesstab);
	  queueAdd(&latestQueue, ftgtaccesstab[i], -(float)nextaccess);
	}
      if(allocmem > maxalloc)
	maxalloc = allocmem;


      solvmtx->ftgttab[ftgtaccesstab[i]].infotab[FTGT_CTRBCNT]--;
      if(solvmtx->ftgttab[ftgtaccesstab[i]].infotab[FTGT_CTRBCNT] == 0)
	{
	  queueAdd(&toSendQueue, ftgtaccesstab[i], 
		   (double)(-ftgtsizetab[ftgtaccesstab[i]]));
	  /*fprintf(stdout, "Ftgt %ld acess %ld \n", (long)ftgtaccesstab[i], (long)i);*/
	  procdst = solvmtx->ftgttab[ftgtaccesstab[i]].infotab[FTGT_PROCDST];
	  while((ftgtcursor[procdst] < ftgtsendnbr[procdst]) && (solvmtx->ftgttab[ftgtsendtab[procdst][ftgtcursor[procdst]]].infotab[FTGT_CTRBCNT] == 0))
	     {
#ifdef DEBUG_BLEND
	       debugtab[ftgtsendtab[procdst][ftgtcursor[procdst]]] = i;
#endif
	       senttab[ftgtsendtab[procdst][ftgtcursor[procdst]]] = 1;
	       /*if(extraftgtcurtab[ftgtsendtab[procdst][ftgtcursor[procdst]]] < 0)*/
		 allocmem -= ftgtsizetab[ftgtsendtab[procdst][ftgtcursor[procdst]]];
		 ftgtcursor[procdst]++;
	     }
	}
    }
  
#ifdef DEBUG_BLEND
  /*for(i=0;i<solvmtx->ftgtnbr;i++)
    fprintf(stdout, "ftgtctrb %ld %ld \n", (long)solvmtx->ftgttab[i].infotab[FTGT_CTRBCNT], (long)solvmtx->ftgttab[i].infotab[FTGT_CTRBNBR]);*/
  /*ASSERT(solvmtx->ftgttab[i].infotab[FTGT_CTRBCNT]==0,MOD_BLEND);*/
  /*if(!(solvmtx->ftgttab[i].infotab[FTGT_CTRBCNT]==0))
    fprintf(stdout, "ftgtnbr %ld ftgtctrb[%ld] %ld %ld \n", (long)solvmtx->ftgtnbr, (long)i, (long)solvmtx->ftgttab[i].infotab[FTGT_CTRBCNT], (long)solvmtx->ftgttab[i].infotab[FTGT_CTRBNBR]);*/
  
  ASSERT(allocmem == 0,MOD_BLEND);
  fprintf(stdout, "Constructing newftgttab \n"); 
#endif
  fprintf(stdout, "ftgtnbr %ld Limit %ld Maxalloc %ld overmem %ld extraftgtnbr %ld \n", (long)solvmtx->ftgtnbr, (long)maxmem, (long)maxalloc, (long)overmem, (long)extraftgtnbr);
  fprintf(stdout, "REAL REDUCTION %g percents of the total AUB \n", (((float)maxalloc)/totalftgtsize)*100);
  fprintf(stdout, "Now AUB/COEF in dynamic execution %g \n", (((float)maxalloc)/(solvmtx->coefnbr*sizeof(double)))*100);
  /** Generation of the new fan-in-target **/
  /* For avoid side effect */
  queueAdd(&partialFtgtInd, solvmtx->ftgtnbr, (double)(solvmtx->ftgtnbr));
  MALLOC_INTERN(newftgttab, solvmtx->ftgtnbr + extracursor, FanInTarget);
  MALLOC_INTERN(newnumtab,  solvmtx->ftgtnbr, pastix_int_t);
  cursor = 0;
  j = 0;
  for(i=0;i<solvmtx->ftgtnbr;i++)
    {
      while(queueSize(&partialFtgtInd) && i == queueRead(&partialFtgtInd))
	{
	  queueGet(&partialFtgtInd);
	  j = queueGet(&partialFtgtExtraInd);
#ifdef DEBUG_BLEND
	  ASSERT( extraftgttab[j].indnum == i,MOD_BLEND);
#endif
	  /*fprintf(stdout, "Add Partial %ld at %ld before Ftgt %ld \n", (long)j, (long)cursor, (long)i);*/
	  /*memcpy(&(newftgttab[cursor]), &(solvmtx->ftgttab[extraftgttab[j].ftgtnum]), sizeof(FanInTarget));*/
	  newftgttab[cursor].infotab[FTGT_PROCDST] = solvmtx->ftgttab[extraftgttab[j].ftgtnum].infotab[FTGT_PROCDST];
	  newftgttab[cursor].infotab[FTGT_TASKDST] = solvmtx->ftgttab[extraftgttab[j].ftgtnum].infotab[FTGT_TASKDST];
	  newftgttab[cursor].infotab[FTGT_BLOKDST] = solvmtx->ftgttab[extraftgttab[j].ftgtnum].infotab[FTGT_BLOKDST];
	  newftgttab[cursor].infotab[FTGT_FCOLNUM] = solvmtx->ftgttab[extraftgttab[j].ftgtnum].infotab[FTGT_FCOLNUM];
	  newftgttab[cursor].infotab[FTGT_LCOLNUM] = solvmtx->ftgttab[extraftgttab[j].ftgtnum].infotab[FTGT_LCOLNUM];
	  newftgttab[cursor].infotab[FTGT_FROWNUM] = solvmtx->ftgttab[extraftgttab[j].ftgtnum].infotab[FTGT_FROWNUM];
	  newftgttab[cursor].infotab[FTGT_LROWNUM] = solvmtx->ftgttab[extraftgttab[j].ftgtnum].infotab[FTGT_LROWNUM];
	  newftgttab[cursor].infotab[FTGT_PRIONUM] = extraftgttab[j].prionum;
	  newftgttab[cursor].infotab[FTGT_CTRBNBR] = extraftgttab[j].ctrbnbr;
#ifdef DEBUG_BLEND
	  ASSERT( extraftgttab[j].ctrbnbr > 0,MOD_BLEND);
#endif
	  extraftgttab[j].ftgtnewnum = cursor;
	  cursor++;
	}
      /*memcpy(&(newftgttab[cursor]), &(solvmtx->ftgttab[i]), sizeof(FanInTarget));*/
      
      newftgttab[cursor].infotab[FTGT_PROCDST] = solvmtx->ftgttab[i].infotab[FTGT_PROCDST];
      newftgttab[cursor].infotab[FTGT_TASKDST] = solvmtx->ftgttab[i].infotab[FTGT_TASKDST];
      newftgttab[cursor].infotab[FTGT_BLOKDST] = solvmtx->ftgttab[i].infotab[FTGT_BLOKDST];
      newftgttab[cursor].infotab[FTGT_FCOLNUM] = solvmtx->ftgttab[i].infotab[FTGT_FCOLNUM];
      newftgttab[cursor].infotab[FTGT_LCOLNUM] = solvmtx->ftgttab[i].infotab[FTGT_LCOLNUM];
      newftgttab[cursor].infotab[FTGT_FROWNUM] = solvmtx->ftgttab[i].infotab[FTGT_FROWNUM];
      newftgttab[cursor].infotab[FTGT_LROWNUM] = solvmtx->ftgttab[i].infotab[FTGT_LROWNUM];
      newftgttab[cursor].infotab[FTGT_PRIONUM] = solvmtx->ftgttab[i].infotab[FTGT_PRIONUM];
      newftgttab[cursor].infotab[FTGT_CTRBNBR] = solvmtx->ftgttab[i].infotab[FTGT_CTRBNBR];
      newnumtab[i] = cursor;
      cursor++;
    }
  memFree_null(solvmtx->ftgttab);
  solvmtx->ftgttab = newftgttab;
#ifdef DEBUG_BLEND
  /*fprintf(stdout, "Newftgttab constructed \n" );*/
  ASSERT(queueSize(&partialFtgtInd) == 1,MOD_BLEND);
  MALLOC_INTERN(debugnew2old, solvmtx->ftgtnbr+extracursor, pastix_int_t);
  /*memset(debugnew2old, -1, sizeof(pastix_int_t)*(solvmtx->ftgtnbr+extracursor));*/
  for(i=0;i<solvmtx->ftgtnbr+extracursor;i++)
    debugnew2old[i] = -1;
  for(i=0;i<solvmtx->ftgtnbr;i++)
    debugnew2old[newnumtab[i]] = i;
  bzero(ftgtsendnbr,  solvmtx->clustnbr*sizeof(pastix_int_t));
#endif
  for(i=0;i<solvmtx->ftgtnbr;i++)
      extraftgtcurtab[i] = extraftgtfirtab[i];
  solvmtx->ftgtnbr += extracursor;
#ifdef DEBUG_BLEND
  for(p=0;p<solvmtx->clustnbr;p++)
    {
      ftgtsendnbr[p] = 0;
      ftgtcursor[p] = 0;
    }
  
  for(i=0;i<solvmtx->ftgtnbr;i++)
    {
      ftgtsendnbr[solvmtx->ftgttab[i].infotab[FTGT_PROCDST]]++;
    }
  for(p=0;p<solvmtx->clustnbr;p++)
    if(ftgtsendnbr[p] > 0)
      {
	MALLOC_INTERN(ftgtsendtab[p], ftgtsendnbr[p], pastix_int_t);
      }
  bzero(ftgtsendnbr,  solvmtx->clustnbr*sizeof(pastix_int_t));
  for(i=0;i<solvmtx->ftgtnbr;i++)
    {
      ftgtsendtab[solvmtx->ftgttab[i].infotab[FTGT_PROCDST]][ftgtsendnbr[solvmtx->ftgttab[i].infotab[FTGT_PROCDST]]] = i;
      ftgtsendnbr[solvmtx->ftgttab[i].infotab[FTGT_PROCDST]]++;
    }
#endif
  
 
  
  /** Update the index of ftgt in the tasks **/
#ifdef DEBUG_BLEND
  for(i=0;i<solvmtx->ftgtnbr;i++)
    solvmtx->ftgttab[i].infotab[FTGT_CTRBCNT] = solvmtx->ftgttab[i].infotab[FTGT_CTRBNBR];
  for(i=0;i<extracursor;i++)
    ASSERT(newnumtab[extraftgttab[i].indnum] > extraftgttab[i].ftgtnewnum,MOD_BLEND);

  debugaccess = 0;
#endif

  queueClear(&taskQueue);
  for(i=0;i<solvmtx->tasknbr;i++)
    queueAdd(&taskQueue, i, (double)(solvmtx->tasktab[i].prionum));
  while(queueSize(&taskQueue)>0)
    {
      pastix_int_t oldindnum;
      i = queueGet(&taskQueue);
#ifdef DEBUG_BLEND
      fprintf(stdout, "S: task %ld prionum %ld \n", (long)i, (long)solvmtx->tasktab[i].prionum);
#endif
      if( (solvmtx->tasktab[i].taskid == E2) && (solvmtx->indtab[solvmtx->tasktab[i].indnum] >= 0) )
	{
	  
	  oldindnum = solvmtx->indtab[solvmtx->tasktab[i].indnum];
	  if(extraftgtcurtab[oldindnum] < 0)
	    {
	      solvmtx->indtab[solvmtx->tasktab[i].indnum] = newnumtab[oldindnum];
#ifdef DEBUG_BLEND
	      ASSERT(solvmtx->ftgttab[newnumtab[oldindnum]].infotab[FTGT_CTRBCNT] > 0,MOD_BLEND);
	      solvmtx->ftgttab[newnumtab[oldindnum]].infotab[FTGT_CTRBCNT]--;
	      if(solvmtx->ftgttab[newnumtab[oldindnum]].infotab[FTGT_CTRBCNT] == 0)
		{
		  fprintf(stdout, "Ftgt %ld access %ld \n", (long)oldindnum, (long)debugaccess);
		  procdst = solvmtx->ftgttab[newnumtab[oldindnum]].infotab[FTGT_PROCDST];
		  while((ftgtcursor[procdst] < ftgtsendnbr[procdst]) && (solvmtx->ftgttab[ftgtsendtab[procdst][ftgtcursor[procdst]]].infotab[FTGT_CTRBCNT] == 0))
		    {
		      fprintf(stdout, "S: ftgt %ld \n", (long)ftgtsendtab[procdst][ftgtcursor[procdst]]);
		      if(debugnew2old[ftgtsendtab[procdst][ftgtcursor[procdst]]] >= 0)
			{
			  fprintf(stdout, "Old Ftgt %ld send at access %ld \n", (long)debugnew2old[ftgtsendtab[procdst][ftgtcursor[procdst]]], (long)debugaccess);
			  ASSERT(debugtab[debugnew2old[ftgtsendtab[procdst][ftgtcursor[procdst]]]] == debugaccess,MOD_BLEND);
			}
		      else
			fprintf(stdout, "ExtraFtgt %ld send at access %ld \n", (long)ftgtsendtab[procdst][ftgtcursor[procdst]], (long)debugaccess);
		      ftgtcursor[procdst]++;
		    }
		}
#endif
	    }
	  else
	    {
	      solvmtx->indtab[solvmtx->tasktab[i].indnum] = extraftgttab[extraftgtcurtab[oldindnum]].ftgtnewnum;
	      extraftgttab[extraftgtcurtab[oldindnum]].ctrbcnt--;
#ifdef DEBUG_BLEND
	      solvmtx->ftgttab[extraftgttab[extraftgtcurtab[oldindnum]].ftgtnewnum].infotab[FTGT_CTRBCNT]--;
#endif
	      if(extraftgttab[extraftgtcurtab[oldindnum]].ctrbcnt == 0)
		{
		  /*fprintf(stdout, "Finish Partial ftgt %ld issue from %ld \n", (long)extraftgtcurtab[oldindnum], (long)oldindnum);*/
		  /*fprintf(stdout, "Next %ld \n", (long)extraftgttab[extraftgtcurtab[oldindnum]].next);*/
#ifdef DEBUG_BLEND
		  ASSERT(solvmtx->ftgttab[newnumtab[extraftgttab[extraftgtcurtab[oldindnum]].indnum]].infotab[FTGT_CTRBCNT] > 0,MOD_BLEND);
		  procdst = solvmtx->ftgttab[newnumtab[extraftgttab[extraftgtcurtab[oldindnum]].indnum]].infotab[FTGT_PROCDST];
		  while((ftgtcursor[procdst] < ftgtsendnbr[procdst]) && (solvmtx->ftgttab[ftgtsendtab[procdst][ftgtcursor[procdst]]].infotab[FTGT_CTRBCNT] == 0))
		    {
		      fprintf(stdout, "S: ftgt %ld \n", (long)ftgtsendtab[procdst][ftgtcursor[procdst]]);
		      if(debugnew2old[ftgtsendtab[procdst][ftgtcursor[procdst]]] >= 0)
			{
			  fprintf(stdout, "Old Ftgt %ld send at access %ld \n", (long)debugnew2old[ftgtsendtab[procdst][ftgtcursor[procdst]]], (long)debugaccess);
			  ASSERT(debugtab[debugnew2old[ftgtsendtab[procdst][ftgtcursor[procdst]]]] == debugaccess,MOD_BLEND);
			}
		      else
			  fprintf(stdout, "ExtraFtgt %ld send at access %ld \n", (long)ftgtsendtab[procdst][ftgtcursor[procdst]], (long)debugaccess);
		      ftgtcursor[procdst]++;
		    }
#endif
		  extraftgtcurtab[oldindnum] = extraftgttab[extraftgtcurtab[oldindnum]].next;
		}
	    }
#ifdef DEBUG_BLEND
	  debugaccess++;
#endif
	}
      if(solvmtx->tasktab[i].taskid == COMP_1D)
	{
	  pastix_int_t odb_nbr;
	  odb_nbr = (solvmtx->cblktab[solvmtx->tasktab[i].cblknum+1].fblokptr - solvmtx->bloktab) - solvmtx->tasktab[i].bloknum -1;
	  for(j=solvmtx->tasktab[i].indnum;j<solvmtx->tasktab[i].indnum + (odb_nbr*(odb_nbr+1))/2;j++)
	    if(solvmtx->indtab[j] >= 0) 
	      {
		oldindnum = solvmtx->indtab[j];
		if(extraftgtcurtab[oldindnum] < 0)
		  {
		    solvmtx->indtab[j] = newnumtab[oldindnum];
#ifdef DEBUG_BLEND
		    ASSERT(solvmtx->ftgttab[newnumtab[oldindnum]].infotab[FTGT_CTRBNBR] > 0,MOD_BLEND);
		    ASSERT(solvmtx->ftgttab[newnumtab[oldindnum]].infotab[FTGT_CTRBCNT] > 0,MOD_BLEND);
		    solvmtx->ftgttab[newnumtab[oldindnum]].infotab[FTGT_CTRBCNT]--;
		    if(solvmtx->ftgttab[newnumtab[oldindnum]].infotab[FTGT_CTRBCNT] == 0)
		      {
			fprintf(stdout, "Ftgt %ld access %ld \n", (long)oldindnum, (long)debugaccess);

			procdst = solvmtx->ftgttab[newnumtab[oldindnum]].infotab[FTGT_PROCDST];
			while((ftgtcursor[procdst] < ftgtsendnbr[procdst]) && (solvmtx->ftgttab[ftgtsendtab[procdst][ftgtcursor[procdst]]].infotab[FTGT_CTRBCNT] == 0))
			  {
			    fprintf(stdout, "S: ftgt %ld \n", (long)ftgtsendtab[procdst][ftgtcursor[procdst]]);
			    if(debugnew2old[ftgtsendtab[procdst][ftgtcursor[procdst]]] >= 0)
			      {
				fprintf(stdout, "Old Ftgt %ld send at access %ld \n", (long)debugnew2old[ftgtsendtab[procdst][ftgtcursor[procdst]]], (long)debugaccess);
				ASSERT(debugtab[debugnew2old[ftgtsendtab[procdst][ftgtcursor[procdst]]]] == debugaccess,MOD_BLEND);
			      }
			    else
			      fprintf(stdout, "ExtraFtgt %ld send at access %ld \n", (long)ftgtsendtab[procdst][ftgtcursor[procdst]], (long)debugaccess);
			    ftgtcursor[procdst]++;
			  }
		      }
#endif
		  }
		else
		  {
		    solvmtx->indtab[j] = extraftgttab[extraftgtcurtab[oldindnum]].ftgtnewnum;
		    extraftgttab[extraftgtcurtab[oldindnum]].ctrbcnt--;
#ifdef DEBUG_BLEND
		    solvmtx->ftgttab[extraftgttab[extraftgtcurtab[oldindnum]].ftgtnewnum].infotab[FTGT_CTRBCNT]--;
#endif
		    if(extraftgttab[extraftgtcurtab[oldindnum]].ctrbcnt == 0)
		      {
			/*fprintf(stdout, "Finish Partial ftgt %ld issue from %ld inserted before %ld \n", (long)extraftgtcurtab[oldindnum], (long)oldindnum, (long)extraftgttab[extraftgtcurtab[oldindnum]].indnum);*/

			/*fprintf(stdout, "Next %ld \n", (long)extraftgttab[extraftgtcurtab[oldindnum]].next);*/
			
#ifdef DEBUG_BLEND
			ASSERT(solvmtx->ftgttab[newnumtab[extraftgttab[extraftgtcurtab[oldindnum]].indnum]].infotab[FTGT_CTRBCNT] > 0,MOD_BLEND);
			procdst = solvmtx->ftgttab[newnumtab[extraftgttab[extraftgtcurtab[oldindnum]].indnum]].infotab[FTGT_PROCDST];
			while((ftgtcursor[procdst] < ftgtsendnbr[procdst]) 
			      && (solvmtx->ftgttab[ftgtsendtab[procdst][ftgtcursor[procdst]]].infotab[FTGT_CTRBCNT] == 0))
			  {
			    fprintf(stdout, "S: ftgt %ld \n", (long)ftgtsendtab[procdst][ftgtcursor[procdst]]);
			    if(debugnew2old[ftgtsendtab[procdst][ftgtcursor[procdst]]] >= 0)
			      {
				fprintf(stdout, "Ftgt %ld send at access %ld \n", (long)debugnew2old[ftgtsendtab[procdst][ftgtcursor[procdst]]], (long)debugaccess);
				ASSERT(debugtab[debugnew2old[ftgtsendtab[procdst][ftgtcursor[procdst]]]] == debugaccess,MOD_BLEND);
			      }
			    else
			      fprintf(stdout, "ExtraFtgt %ld send at access %ld \n", (long)ftgtsendtab[procdst][ftgtcursor[procdst]], (long)debugaccess);
			    ftgtcursor[procdst]++;
			  }
#endif
			extraftgtcurtab[oldindnum] = extraftgttab[extraftgtcurtab[oldindnum]].next;
		      }
		  }
#ifdef DEBUG_BLEND
		debugaccess++;
#endif
	      }
	}
    }
#ifdef DEBUG_BLEND
  for(i=0;i<solvmtx->clustnbr;i++)
    {
      fprintf(stdout, "i %ld ftgtcursor%ld ftgtsendnbr %ld \n", (long)i, (long)ftgtcursor[i], (long)ftgtsendnbr[i]);
      ASSERT(ftgtcursor[i] == ftgtsendnbr[i],MOD_BLEND);
    }
  for(i=0;i<solvmtx->ftgtnbr;i++)
    solvmtx->ftgttab[i].infotab[FTGT_CTRBCNT] == 0;
#endif
  
  /** Reset the ctrcnt of the ftgt **/
  for(i=0;i<solvmtx->ftgtnbr;i++)
    solvmtx->ftgttab[i].infotab[FTGT_CTRBCNT] = solvmtx->ftgttab[i].infotab[FTGT_CTRBNBR];

#ifdef DEBUG_BLEND
  {
    pastix_int_t debugaccess = 0;
    for(i=0;i<solvmtx->ftgtnbr;i++)
      debugaccess += solvmtx->ftgttab[i].infotab[FTGT_CTRBNBR];

    if(debugaccess != accessnbr)
      fprintf(stdout, "Debugaccess %ld accesnbr %ld \n", (long)debugaccess, (long)accessnbr);
    ASSERT(debugaccess == accessnbr,MOD_BLEND);
  }
#endif

    
#ifdef DEBUG_BLEND

  fprintf(stdout, "Cursor %ld total %ld \n", (long)cursor, (long)solvmtx->ftgtnbr);
#endif
 
  /** Desallocations **/
  memFree_null(newnumtab);
  if(extracursor > 0)
    memFree_null(extraftgttab);
  queueExit(&taskQueue);
  queueExit(&partialFtgtInd);
  queueExit(&partialFtgtExtraInd);
  queueExit(&latestQueue);
  queueExit(&toSendQueue);
#ifdef DEBUG_BLEND
  memFree_null(debugtab);
#endif
  memFree_null(senttab);
  memFree_null(ftgtsizetab);
  memFree_null(ftgtaccesstab);
  for(p=0;p<solvmtx->clustnbr;p++)
    if(ftgtsendnbr[p]>0)
      memFree_null(ftgtsendtab[p]);
  memFree_null(ftgtsendtab);
  memFree_null(ftgtsendnbr);
  memFree_null(ftgtcursor);
  memFree_null(extraftgtcurtab);
  memFree_null(extraftgtfirtab);
  /*fprintf(stdout, " accessnbr %ld maxalloc %ld \n", (long)accessnbr, (long)maxalloc);*/
  fprintf(stdout, "Percentage of AUB dynamique / totalcoef %g \n", ((float)maxalloc)/(solvmtx->coefnbr*sizeof(double)));
  return maxalloc;

}

pastix_int_t getFtgtInd2(SolverMatrix *solvmtx, pastix_int_t *senttab, Queue *toSendQueue, Queue *latestQueue)
{
  pastix_int_t ftgtnum = -1;
  pastix_int_t ind;
  while(queueSize(toSendQueue)>0)
    {
      ind = queueGet(toSendQueue);
#ifdef DEBUG_BLEND
      ASSERT(solvmtx->ftgttab[ind].infotab[FTGT_CTRBCNT] == 0,MOD_BLEND);
#endif
      if(senttab[ind] == 0)
	{
	    /*fprintf(stdout, "Se"); */
	  ftgtnum = ind;
	  break;
	}
    }
  if(ftgtnum <0)
    while(queueSize(latestQueue)>0)
      {
	ind = queueGet(latestQueue);
	/*ASSERT(solvmtx->ftgttab[ind].infotab[FTGT_CTRBCNT] < solvmtx->ftgttab[ind].infotab[FTGT_CTRBNBR],MOD_BLEND);*/
	if(solvmtx->ftgttab[ind].infotab[FTGT_CTRBCNT] > 0) /* OIMBE >= 0 aussi ?? **/
	  /*&&
	 (solvmtx->ftgttab[ind].infotab[FTGT_CTRBCNT] < solvmtx->ftgttab[ind].infotab[FTGT_CTRBNBR]))*/
	  {
	    ftgtnum = ind;
	    break;
	  }
      }
  
#ifdef DEBUG_BLEND
  ASSERT(ftgtnum >=0,MOD_BLEND);
#endif
  return ftgtnum;
}

pastix_int_t getFtgtNextAccess(pastix_int_t ind, pastix_int_t accessnbr, pastix_int_t *ftgtaccesstab)
{
  pastix_int_t i;
  pastix_int_t maxind = -1;
  for(i=ind;i<accessnbr;i++)
    {
      if(ftgtaccesstab[i] == ftgtaccesstab[ind])
	maxind = i;
    }
  return maxind;
}


