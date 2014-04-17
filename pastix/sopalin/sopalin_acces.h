#define SYMB_CBLKNBR      datacode->cblknbr
#define SYMB_BLOKNBR      datacode->bloknbr
#define SYMB_NODENBR      datacode->nodenbr   /* Is it useful ?*/
#define SYMB_BLOKNUM(x)   (datacode->cblktab[x].fblokptr - datacode->bloktab)
#define SYMB_FCOLNUM(x)   datacode->cblktab[x].fcolnum
#define SYMB_LCOLNUM(x)   datacode->cblktab[x].lcolnum
#define SYMB_FROWNUM(x)   datacode->bloktab[x].frownum
#define SYMB_LROWNUM(x)   datacode->bloktab[x].lrownum
#define SYMB_CBLKNUM(x)   datacode->bloktab[x].cblknum /*<0 if remote*/

#define CBLK_BLOKNBR(x)   (SYMB_BLOKNUM(x+1) - SYMB_BLOKNUM(x))
#define CBLK_COLNBR(x)    (SYMB_LCOLNUM(x) - SYMB_FCOLNUM(x) + 1)
#define BLOK_ROWNBR(x)    (SYMB_LROWNUM(x) - SYMB_FROWNUM(x) + 1)
#define HBLOCK_FROWNUM(x)  datacode->hbloktab[x].frownum
#define HBLOCK_LROWNUM(x)  datacode->hbloktab[x].lrownum
#define HBLOCK_ROWNBR(x)   (HBLOCK_LROWNUM(x) - HBLOCK_FROWNUM(x) + 1)
#define HBLOCK_COEFIND(x)  datacode->hbloktab[x].coefind
#ifdef NAPA_SOPALIN /* ILU(k) */
#define HBLOCK_ISFACING(j,b)                            \
  (((HBLOCK_FROWNUM(j)>=SYMB_FROWNUM(b)) &&             \
    (HBLOCK_LROWNUM(j)<=SYMB_LROWNUM(b))) ||            \
   ((HBLOCK_FROWNUM(j)<=SYMB_FROWNUM(b)) &&             \
    (HBLOCK_LROWNUM(j)>=SYMB_LROWNUM(b))) ||            \
   ((HBLOCK_FROWNUM(j)<=SYMB_FROWNUM(b)) &&             \
    (HBLOCK_LROWNUM(j)>=SYMB_FROWNUM(b))) ||            \
   ((HBLOCK_FROWNUM(j)<=SYMB_LROWNUM(b)) &&             \
    (HBLOCK_LROWNUM(j)>=SYMB_LROWNUM(b))))
#else
#define HBLOCK_ISFACING(j,b)                            \
  ((HBLOCK_FROWNUM(j)>=SYMB_FROWNUM(b)) &&              \
   (HBLOCK_LROWNUM(j)<=SYMB_LROWNUM(b)))
#endif
#ifdef NAPA_SOPALIN /* ILU(k) */
#define BLOCK_ISFACING(j,b)                           \
  (((SYMB_FROWNUM(j)>=SYMB_FROWNUM(b)) &&             \
    (SYMB_LROWNUM(j)<=SYMB_LROWNUM(b))) ||            \
   ((SYMB_FROWNUM(j)<=SYMB_FROWNUM(b)) &&             \
    (SYMB_LROWNUM(j)>=SYMB_LROWNUM(b))) ||            \
   ((SYMB_FROWNUM(j)<=SYMB_FROWNUM(b)) &&             \
    (SYMB_LROWNUM(j)>=SYMB_FROWNUM(b))) ||            \
   ((SYMB_FROWNUM(j)<=SYMB_LROWNUM(b)) &&             \
    (SYMB_LROWNUM(j)>=SYMB_LROWNUM(b))))
#else
#define BLOCK_ISFACING(j,b)                           \
  ((SYMB_FROWNUM(j)>=SYMB_FROWNUM(b)) &&              \
   (SYMB_LROWNUM(j)<=SYMB_LROWNUM(b)))
#endif

#define SOLV_COEFMAX      datacode->coefmax /* max elements for greatest cblk*/
#define SOLV_COEFNBR      datacode->coefnbr
#define SOLV_FTGTCNT      datacode->ftgtcnt
#define SOLV_FTGTNBR      datacode->ftgtnbr
#define SOLV_INDNBR       datacode->indnbr
#define SOLV_INDTAB       datacode->indtab
#define SOLV_PROCNBR      datacode->clustnbr
#define SOLV_PROCNUM      datacode->clustnum
#define SOLV_TASKNBR      datacode->tasknbr
#define SOLV_TASKTAB      datacode->tasktab

#define SOLV_COEFIND(x)   datacode->bloktab[x].coefind

#define SOLV_COEFTAB(x)   datacode->cblktab[x].coeftab
#define SOLV_PROCDIAG(x)  datacode->cblktab[x].procdiag
#define SOLV_STRIDE(x)    datacode->cblktab[x].stride
#define SOLV_UCOEFTAB(x)  datacode->cblktab[x].ucoeftab

#define SOLV_TTSKTAB(x)    datacode->ttsktab[me][x]
#define SOLV_TTSKNBR       datacode->ttsknbr[me]
#define SOLV_PROC2CLUST(x) datacode->proc2clust[x]
#ifdef SMP_SOPALIN
#define SOLV_THRDNBR       datacode->thrdnbr
#define SOLV_BUBLNBR       datacode->bublnbr
#else
#define SOLV_THRDNBR       1
#define SOLV_BUBLNBR       1
#endif /* SMP_SOPALIN */


#define FANIN_COEFTAB(x)  datacode->ftgttab[x].coeftab
#define FANIN_INFOTAB(x)  datacode->ftgttab[x].infotab

#define FANIN_BLOKDST(x)  datacode->ftgttab[x].infotab[FTGT_BLOKDST]
#define FANIN_CTRBCNT(x)  datacode->ftgttab[x].infotab[FTGT_CTRBCNT]
#define FANIN_CTRBNBR(x)  datacode->ftgttab[x].infotab[FTGT_CTRBNBR]
#define FANIN_FCOLNUM(x)  datacode->ftgttab[x].infotab[FTGT_FCOLNUM]
#define FANIN_FROWNUM(x)  datacode->ftgttab[x].infotab[FTGT_FROWNUM]
#define FANIN_LCOLNUM(x)  datacode->ftgttab[x].infotab[FTGT_LCOLNUM]
#define FANIN_LROWNUM(x)  datacode->ftgttab[x].infotab[FTGT_LROWNUM]
#define FANIN_PRIONUM(x)  datacode->ftgttab[x].infotab[FTGT_PRIONUM]
#define FANIN_TASKDST(x)  datacode->ftgttab[x].infotab[FTGT_TASKDST]
#define FANIN_GCBKDST(x)  datacode->ftgttab[x].infotab[FTGT_GCBKDST]  /* Only in OOC     */
#define FANIN_IDTRACE(x)  datacode->ftgttab[x].infotab[FTGT_IDTRACE]  /* Only with trace */

#define FANIN_CBLKDST(x)  FANIN_TASKDST(x)

#ifdef SMP_SOPALIN
#define FANIN_PROCDST(x)  SOLV_PROC2CLUST(datacode->ftgttab[x].infotab[FTGT_PROCDST])
#else /* SMP_SOPALIN */
#define FANIN_PROCDST(x)  datacode->ftgttab[x].infotab[FTGT_PROCDST]
#endif /* SMP_SOPALIN */

#define TASK_BLOKNUM(x)   datacode->tasktab[x].bloknum
#define TASK_CBLKNUM(x)   datacode->tasktab[x].cblknum
#define TASK_CTRBCNT(x)   datacode->tasktab[x].ctrbcnt
#define TASK_FTGTCNT(x)   datacode->tasktab[x].ftgtcnt
#define TASK_INDNUM(x)    datacode->tasktab[x].indnum
#define TASK_PRIONUM(x)   datacode->tasktab[x].prionum
#define TASK_PROC(x)      SOLV_PROCDIAG(TASK_CBLKNUM(x))
#define TASK_TASKID(x)    datacode->tasktab[x].taskid
#define TASK_THREADID(x)  datacode->tasktab[x].threadid

#define UPDOWN_SM2XTAB          datacode->updovct.sm2xtab
#define UPDOWN_SM2XMAX          datacode->updovct.sm2xmax
#define UPDOWN_SM2XIND(x)       datacode->updovct.cblktab[x].sm2xind
#define UPDOWN_CTRBCNT(x)       datacode->updovct.cblktab[x].ctrbcnt
#define UPDOWN_CTRBNBR(x)       datacode->updovct.cblktab[x].ctrbnbr
#define UPDOWN_MSGNBR(x)        datacode->updovct.cblktab[x].msgnbr
#define UPDOWN_MSGCNT(x)        datacode->updovct.cblktab[x].msgcnt
#define UPDOWN_BROWPROCTAB(x)   datacode->updovct.cblktab[x].browproctab
#define UPDOWN_BROWCBLKTAB(x)   datacode->updovct.cblktab[x].browcblktab
#define UPDOWN_SM2XSZE          datacode->updovct.sm2xsze
#define UPDOWN_SM2XNBR          datacode->updovct.sm2xnbr
#define UPDOWN_BROWPROCNBR(x)   datacode->updovct.cblktab[x].browprocnbr
#define UPDOWN_GCBLK2LIST(x)    datacode->updovct.gcblk2list[x]
#define UPDOWN_GCBLK2LISTNBR    datacode->updovct.gcblk2listnbr
#define UPDOWN_LISTPTR(x)       datacode->updovct.listptr[x]
#define UPDOWN_LISTPTRNBR       datacode->updovct.listptrnbr
#define UPDOWN_LISTCBLK(x)      datacode->updovct.listcblk[x]
#define UPDOWN_LISTBLOK(x)      datacode->updovct.listblok[x]
#define UPDOWN_LISTNBR          datacode->updovct.listnbr
#define SOLV_GCBLKNBR           datacode->gcblknbr
#define UPDOWN_GCBLK2GLIST(x)   datacode->updovct.gcblk2glist[x]
#define UPDOWN_GCBLK2GLISTNBR   datacode->updovct.gcblk2glistnbr
#define UPDOWN_GLISTPTR(x)      datacode->updovct.glistptr[x]
#define UPDOWN_GLISTPTRNBR      datacode->updovct.glistptrnbr
#define UPDOWN_GLISTCBLK(x)     datacode->updovct.glistcblk[x]
#define UPDOWN_GLISTBLOK(x)     datacode->updovct.glistblok[x]
#define UPDOWN_GLISTPROC(x)     datacode->updovct.glistproc[x]
#define UPDOWN_GLISTNBR         datacode->updovct.glistnbr
#define UPDOWN_LOC2GLOB(x)      datacode->updovct.loc2glob[x]
#define UPDOWN_LOC2GLOBNBR      datacode->updovct.loc2globnbr
#define UPDOWN_LBLK2GCBLK(x)    datacode->updovct.lblk2gcblk[x]
#define UPDOWN_GCBLKNBR         datacode->updovct.gcblknbr
#define UPDOWN_GNODENBR         datacode->updovct.gnodenbr
#define UPDOWN_UPMSGNBR         datacode->updovct.upmsgnbr
#define UPDOWN_DOWNMSGNBR       datacode->updovct.downmsgnbr

/* next info is now stored in task struct */ 
#define SOLV_FTGTIND(x)   (-SYMB_CBLKNUM(x))
/* next should be good for 1D ??? */

#define PACKMAX  datacode->nbftmax
#define PACKAREA datacode->arftmax

#define PASTIX_COMM sopalin_data->sopar->pastix_comm

#ifdef TEST_ISEND
#define MAX_S_REQUESTS                                                \
  ( ( THREAD_FUNNELED_ON )?(8096):                                      \
    (MAX(1, (8096/(SOLV_THRDNBR)))))
#else
#define MAX_S_REQUESTS 1
#endif

#ifdef TEST_IRECV
#define MAX_R_REQUESTS (SOLV_PROCNBR-1) /*+ Inutile de mettre plus de requette que de proc +*/
#else
#define MAX_R_REQUESTS 1
#endif

#define THREAD_FUNNELED_ON (                                    \
    sopalin_data->sopar->iparm[IPARM_THREAD_COMM_MODE] &        \
    API_THREAD_FUNNELED)
#define THREAD_FUNNELED_OFF (!THREAD_FUNNELED_ON)

#define THREAD_COMM_ON  (                                       \
    sopalin_data->sopar->iparm[IPARM_THREAD_COMM_MODE] &        \
    ( API_THREAD_FUNNELED|API_THREAD_COMM_ONE|                   \
      API_THREAD_COMM_DEFINED|API_THREAD_COMM_NBPROC ) )
#define THREAD_COMM_OFF (!THREAD_COMM_ON)

#define TASK_TASK2ESP( __i ) ( -((__i) + 2) )
#define TASK_ESP2TASK( __i ) ( -((__i) + 2) )

#define SOLV_HCBLKNBR     datacode->hcblknbr
#define SOLV_HBLOKNBR     (datacode->hcblktab[SOLV_HCBLKNBR].fblokptr - datacode->hbloktab)
#define HCBLK_BLOKNUM(x)   (datacode->hcblktab[x].fblokptr - datacode->hbloktab)
#define HCBLK_FCOLNUM(x)   datacode->hcblktab[x].fcolnum
#define HCBLK_LCOLNUM(x)   datacode->hcblktab[x].lcolnum
#define HCBLK_COLNBR(x)    (HCBLK_LCOLNUM(x) - HCBLK_FCOLNUM(x) + 1)
#define HCBLK_STRIDE(x)    datacode->hcblktab[x].stride
#define HCBLK_OWNER(x)     datacode->hcblktab[x].procdiag
#define HCBLK_GCBLK(x)     datacode->hcblktab[x].gcblknum
#define SOLV_GCBLK2HALO(x) -(datacode->gcblk2halo[x]+1)
#define SOLV_GCBLK2LOC(x)    datacode->gcblk2halo[x]-1
#define SOLV_GCBLKISHALO(x) (datacode->gcblk2halo[x] < 0)
