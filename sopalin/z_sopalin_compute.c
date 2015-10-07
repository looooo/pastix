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
 File: sopalin_compute.c

 Computation functions.

 Pierre Ramet : fev 2003

 */

void API_CALL(z_factor_diag)  (z_Sopalin_Data_t *sopalin_data, pastix_int_t me, pastix_int_t c);
void API_CALL(z_factor_trsm1d)(z_Sopalin_Data_t *sopalin_data, pastix_int_t me, pastix_int_t c);
void API_CALL(z_compute_contrib_compact)(z_Sopalin_Data_t *sopalin_data, pastix_int_t me, pastix_int_t c, pastix_int_t b1, pastix_int_t b2, pastix_int_t usediag);
void API_CALL(z_add_contrib_local)      (z_Sopalin_Data_t *sopalin_data, pastix_int_t me, pastix_int_t b1,pastix_int_t b2,pastix_int_t c,pastix_int_t b3,pastix_int_t cbl);
void API_CALL(z_add_contrib_target)     (z_Sopalin_Data_t *sopalin_data, pastix_int_t me, pastix_int_t b1,pastix_int_t b2,pastix_int_t task,pastix_int_t t);

/*
 * Compute tasks
 */
void API_CALL(z_compute_1d)    (z_Sopalin_Data_t *sopalin_data, pastix_int_t me, pastix_int_t task);
void API_CALL(z_compute_1dgemm)(z_Sopalin_Data_t *sopalin_data, pastix_int_t me, pastix_int_t task, pastix_int_t i, pastix_int_t b2);

#include "z_compute_gemdm.c"
#include "z_compute_diag.c"
#include "z_compute_trsm.c"

/****************************************************************************/
/* COMPUTE TASK 1D                                                          */
/****************************************************************************/

/*
 * Compute the update in one big block for all the left block column
 * multiply by the one in regards with the diagonal block
 * The result is stored in a temporary buffer to be added part
 * by part to the target cblk
 */
void API_CALL(z_compute_contrib_compact)(z_Sopalin_Data_t *sopalin_data, pastix_int_t me, pastix_int_t c, pastix_int_t b1, pastix_int_t b2, pastix_int_t usediag)
{
    z_SolverMatrix  *datacode    = sopalin_data->datacode;
    z_Thread_Data_t *thread_data = sopalin_data->thread_data[me];
    pastix_complex64_t         *gaik, *gb, *gc;
#ifdef CHOL_SOPALIN
    pastix_complex64_t         *gajk;
#endif
    pastix_int_t            dima, dimi, dimj, stride, k;
    (void)b2; (void)gb; (void)usediag;

    gb = thread_data->maxbloktab1; /* C in U for LU, B for LDLt */
    gc = thread_data->maxbloktab2; /* C in L */
    stride = SOLV_STRIDE(c);

    /* Matrix A = Aik */
    gaik = &(SOLV_COEFTAB(c)[SOLV_COEFIND(b1)]);

    /* Compute M */
    dimi = 0;
    for (k=b1; k<SYMB_BLOKNUM(c+1); k++)
        dimi += SYMB_LROWNUM(k) - SYMB_FROWNUM(k) + 1;

    /* Compute N */
    dimj = SYMB_LROWNUM(b1) - SYMB_FROWNUM(b1) + 1;/* ATTENTION stride-dima;*/
#ifdef PASTIX_ESC
    for (k=b2; k>b1; k--)
        dimj += SYMB_LROWNUM(k) - SYMB_FROWNUM(k) + 1;
#endif

    /* Compute K */
    dima = SYMB_LCOLNUM(c) - SYMB_FCOLNUM(c) + 1;

    thread_data->firstbloktab  = b1;
    thread_data->stridebloktab = dimi;

    ASSERTDBG(dimj*thread_data->stridebloktab <= SOLV_COEFMAX, MOD_SOPALIN);

    /* multAikB(gc,gaik,gb,stride,dimj,dima,dimi) */
#ifdef CHOL_SOPALIN
#ifdef SOPALIN_LU
    /* L update */
    gajk = &(SOLV_UCOEFTAB(c)[SOLV_COEFIND(b1)]);
    SOPALIN_GEMM("N", "T", dimi, dimj, dima,
                 fun,   gaik, stride,
                 gajk, stride,
                 fzero, gc,   dimi);

    /* U update */
    gaik = &(SOLV_UCOEFTAB(c)[SOLV_COEFIND(b1)]);
    gajk = &(SOLV_COEFTAB(c)[ SOLV_COEFIND(b1)]);
    SOPALIN_GEMM("N", "T", dimi, dimj, dima,
                 fun,   gaik, stride,
                 gajk, stride,
                 fzero, gb,   dimi);
#else /* SOPALIN_LU */
    gajk = &(SOLV_COEFTAB(c)[ SOLV_COEFIND(b1)]);
    SOPALIN_GEMM("N",
                 "C",
                 dimi, dimj, dima,
                 fun,   gaik, stride,
                 gajk, stride,
                 fzero, gc,   dimi);
#endif /* SOPALIN_LU */
#else /* CHOL_SOPALIN */
    if ( usediag == 1 )
    {
        int ldw = SOLV_COEFMAX;
        pastix_complex64_t *D = &(SOLV_UCOEFTAB(c)[SOLV_COEFIND(SYMB_BLOKNUM(c))]);
        gb = &(SOLV_COEFTAB(c)[SOLV_COEFIND(b1)]);

        API_CALL(z_CORE_gemdm)(PastixNoTrans,
#ifdef HERMITIAN
                             PastixConjTrans,
#else
                             PastixTrans,
#endif
                             dimi, dimj, dima,
                             fun,   gaik, stride,
                             gb,   stride,
                             fzero, gc,   dimi,
                             D, iun,
                             thread_data->maxbloktab1, ldw );
    }
    else
    {
        gb = thread_data->maxbloktab1+(stride-dima-dimi); /* Correction pour le pere Goudin */
#ifdef HERMITIAN
        SOPALIN_GEMM("N",
                     "C",
                     dimi, dimj, dima,
                     fun,   gaik, stride,
                     gb,   stride-dima,
                     fzero, gc,   dimi);
#else
        SOPALIN_GEMM("N",
                     "T",
                     dimi, dimj, dima,
                     fun,   gaik, stride,
                     gb,   stride-dima,
                     fzero, gc,   dimi);
#endif
    }
#endif /* CHOL_SOPALIN */
}

/*
 Function: API_CALL(z_add_contrib_local)

 Add contribution from b2 on b3.

 Parameters:
 sopalin_data - Global sopalin structure
 me           - Computing thread number
 b1           - extra digonal block of c
 b2           - extra diagonal block after b1 on c
 c            - Column block
 b3           - Extra digonal block on cbl facing b2
 cbl          - Column block facing b1

 */
void API_CALL(z_add_contrib_local)(z_Sopalin_Data_t *sopalin_data, pastix_int_t me, pastix_int_t b1, pastix_int_t b2, pastix_int_t c, pastix_int_t b3, pastix_int_t cbl)
{
    pastix_int_t frownum, lrownum, ofrownum, olrownum;
    pastix_int_t dimi, dimj, stridea, strideb, step, k;
    pastix_complex64_t *ga,*gb;
#ifdef SOPALIN_LU
    pastix_complex64_t *ga2, *gb2;
#endif /* SOPALIN_LU */
    z_SolverMatrix  *datacode    = sopalin_data->datacode;
    z_Thread_Data_t *thread_data = sopalin_data->thread_data[me];
#ifdef DEBUG_SOPALIN_NAPA
    int flag;
#endif
    (void)c;

    /* column block in front of b1 */
    /* cbl=SYMB_CBLKNUM(b1); ??? */
    stridea = SOLV_STRIDE(cbl);

    print_debug(DBG_SOPALIN_COMP1D,
                "FR(b2) %ld FR(b3) %ld\nLR(b3) %ld LR(b2) %ld\nFR(b1) %ld FC(cb) %ld\nLC(cb) %ld LR(b1) %ld\n",
                (long)SYMB_FROWNUM(b2), (long)SYMB_FROWNUM(b3),
                (long)SYMB_LROWNUM(b3), (long)SYMB_LROWNUM(b2),
                (long)SYMB_FROWNUM(b1), (long)SYMB_FCOLNUM(cbl),
                (long)SYMB_LCOLNUM(cbl),(long)SYMB_LROWNUM(b1));

#ifdef NAPA_SOPALIN
    ASSERTDBG((SYMB_FROWNUM(b1)  >= SYMB_FCOLNUM(cbl)) &&
              (SYMB_LCOLNUM(cbl) >= SYMB_LROWNUM(b1) ), MOD_SOPALIN);
#else
    ASSERTDBG((SYMB_FROWNUM(b2)  >= SYMB_FROWNUM(b3) ) &&
              (SYMB_LROWNUM(b3)  >= SYMB_LROWNUM(b2) ) &&
              (SYMB_FROWNUM(b1)  >= SYMB_FCOLNUM(cbl)) &&
              (SYMB_LCOLNUM(cbl) >= SYMB_LROWNUM(b1) ), MOD_SOPALIN);
#endif

    ga = &(SOLV_COEFTAB(cbl)[ SOLV_COEFIND(b3)+
                              (SYMB_FROWNUM(b1)-SYMB_FCOLNUM(cbl))*stridea+
                              (SYMB_FROWNUM(b2)-SYMB_FROWNUM(b3))]);
#ifdef SOPALIN_LU
    if ((b3 == SYMB_BLOKNUM(cbl)) && (b1 != b2)) {
        ga2 = &(SOLV_COEFTAB(cbl)[SOLV_COEFIND(b3)+
                                  (SYMB_FROWNUM(b1)-SYMB_FCOLNUM(cbl))+
                                  (SYMB_FROWNUM(b2)-SYMB_FROWNUM(b3))*stridea]);
    }
    else {
        ga2 = &(SOLV_UCOEFTAB(cbl)[SOLV_COEFIND(b3)+
                                   (SYMB_FROWNUM(b1)-SYMB_FCOLNUM(cbl))*stridea+
                                   (SYMB_FROWNUM(b2)-SYMB_FROWNUM(b3))]);
    }
#endif /* SOPALIN_LU */

    /* vertical dimension */
    dimj = SYMB_LROWNUM(b1) - SYMB_FROWNUM(b1) + 1;
    /* vertical dimension */
    dimi = SYMB_LROWNUM(b2) - SYMB_FROWNUM(b2) + 1;

    strideb = thread_data->stridebloktab;

    step = b2 - b1;
    for (k=b1; k<b2; k++)
        step += SYMB_LROWNUM(k) - SYMB_FROWNUM(k);

    ASSERTDBG(step+dimi    <= strideb,     MOD_SOPALIN);
    ASSERTDBG(strideb*dimj <= SOLV_COEFMAX,MOD_SOPALIN);

#ifdef PASTIX_ESC
    /* Décalage du bloc colonne de k + decalage du bloc k inutile */
    for (k=thread_data->firstbloktab; k<b1; k++)
        step += (SYMB_LROWNUM(k) - SYMB_FROWNUM(k) + 1) * (strideb + 1);
#endif

    gb  = &(thread_data->maxbloktab2[step]);
#ifdef SOPALIN_LU
    gb2 = &(thread_data->maxbloktab1[step]);
#endif /* SOPALIN_LU */

#ifdef NAPA_SOPALIN
    ofrownum=SYMB_FROWNUM(b2);
    olrownum=SYMB_LROWNUM(b2);
    b3--;
#ifdef DEBUG_SOPALIN_NAPA
    flag = 1;
#endif
    do {
#ifdef DEBUG_SOPALIN_NAPA
        pastix_int_t trace = 0;
        /* il peut y avoir plusieurs cibles partielles */
        if (!flag)
        {
            print_debug(DBG_SOPALIN_NAPA, "ILU: plusieurs cibles locales\n");
        }
#endif
        frownum=ofrownum;
        lrownum=olrownum;
        b3++;
#ifdef DEBUG_SOPALIN_NAPA
        if ((!flag) || (SYMB_FROWNUM(b3)>frownum) ||
            (SYMB_LROWNUM(b3)<lrownum))
        {
            trace = 1;
            if (flag)
            {
                print_debug(DBG_SOPALIN_NAPA, "\nILU: debug local SF=%ld F=%ld SL=%ld L=%ld gb=%ld ga=%ld\n",
                            (long)SYMB_FROWNUM(b3), (long)frownum,
                            (long)SYMB_LROWNUM(b3), (long)lrownum,
                            (long)0, (long)(SOLV_COEFIND(b3)+
                                            (SYMB_FROWNUM(b1)-SYMB_FCOLNUM(cbl))*stridea+
                                            (SYMB_FROWNUM(b2)-SYMB_FROWNUM(b3))));
            }
        }
#endif
        if (SYMB_FROWNUM(b3)>frownum)
        {
            frownum=SYMB_FROWNUM(b3);
            print_debug(DBG_SOPALIN_NAPA, "ILU: tronque frownum\n");
        }
        if (SYMB_LROWNUM(b3)<lrownum)
        {
            lrownum=SYMB_LROWNUM(b3);
            print_debug(DBG_SOPALIN_NAPA, "ILU: tronque lrownum\n");
        }
        dimi=lrownum-frownum+1;
        /*
         gb=&(maxbloktab2[me][step])+frownum-ofrownum;
         */
        gb=&(thread_data->maxbloktab2[step+frownum-ofrownum]);
        ga=&(SOLV_COEFTAB(cbl)[SOLV_COEFIND(b3)+
                               (SYMB_FROWNUM(b1)-SYMB_FCOLNUM(cbl))*stridea+
                               (frownum-SYMB_FROWNUM(b3))]);
#ifdef SOPALIN_LU
        if ( b3!=SYMB_BLOKNUM(cbl) )
        {
            ga2=&(SOLV_UCOEFTAB(cbl)[SOLV_COEFIND(b3)+
                                     (SYMB_FROWNUM(b1)-SYMB_FCOLNUM(cbl))*stridea+
                                     (frownum-SYMB_FROWNUM(b3))]);
        }
        else if ( b1 != b2 ) /* Store U and L directly in L for factorization */
        {
            ga2=&(SOLV_COEFTAB(cbl)[SOLV_COEFIND(b3)+
                                    (SYMB_FROWNUM(b1)-SYMB_FCOLNUM(cbl))+
                                    (frownum-SYMB_FROWNUM(b3))*stridea]);
        }
        else {
            /* We are on diagonal block => we don't do anything */
        }
        gb2=&(thread_data->maxbloktab1[step+frownum-ofrownum]);
#endif /* SOPALIN_LU */
#ifdef DEBUG_SOPALIN_NAPA
        if (trace)
        {
            print_debug(DBG_SOPALIN_NAPA,
                        "ILU: debug local SF=%ld F=%ld SL=%ld L=%ld gb=%ld ga=%ld\n",
                        (long)SYMB_FROWNUM(b3), (long)frownum,
                        (long)SYMB_LROWNUM(b3), (long)lrownum,
                        (long)(frownum-ofrownum),
                        (long)(SOLV_COEFIND(b3) +
                               (SYMB_FROWNUM(b1)-SYMB_FCOLNUM(cbl))*stridea+
                               (frownum-SYMB_FROWNUM(b3))));
        }
#endif

        ASSERTDBG((SYMB_FROWNUM(b3)<=frownum) &&
                  (SYMB_LROWNUM(b3)>=lrownum),MOD_SOPALIN);

#else /* NAPA_SOPALIN */

        ASSERTDBG((SYMB_FROWNUM(b3)<=SYMB_FROWNUM(b2)) &&
                  (SYMB_LROWNUM(b3)>=SYMB_LROWNUM(b2)),MOD_SOPALIN);

#endif /* NAPA_SOPALIN */

        MUTEX_LOCK(&(sopalin_data->mutex_blok[b3]));
        SOPALIN_GESM( "N", "N", dimi, dimj,
                      fun, gb, strideb,
                      ga, stridea);
#ifdef SOPALIN_LU
        if ( b3!=SYMB_BLOKNUM(cbl) )
        {
            SOPALIN_GESM("N","N",dimi,dimj,fun,gb2,strideb,ga2,stridea);
        }
        else if ( b1 != b2 ) /* Store U and L directly in L for factorization */
        {
            SOPALIN_GESM("T","N",dimj,dimi,fun,gb2,strideb,ga2,stridea);
        }
        else {
            /* We are on diagonal block => we don't do anything */
        }
#endif
        MUTEX_UNLOCK(&(sopalin_data->mutex_blok[b3]));

#ifdef NAPA_SOPALIN
#ifdef DEBUG_SOPALIN_NAPA
        flag = 0;
#endif
    } while  ((b3+1<SYMB_BLOKNUM(cbl+1)) &&
              (
                  ((SYMB_FROWNUM(b3+1)<=SYMB_FROWNUM(b2)) &&
                   (SYMB_LROWNUM(b3+1)>=SYMB_FROWNUM(b2))) ||
                  ((SYMB_LROWNUM(b3+1)>=SYMB_LROWNUM(b2)) &&
                   (SYMB_FROWNUM(b3+1)<=SYMB_LROWNUM(b2))) ||
                  ((SYMB_FROWNUM(b3+1)>=SYMB_FROWNUM(b2)) &&
                   (SYMB_LROWNUM(b3+1)<=SYMB_LROWNUM(b2))) ||
                  ((SYMB_FROWNUM(b3+1)<=SYMB_FROWNUM(b2)) &&
                   (SYMB_LROWNUM(b3+1)>=SYMB_LROWNUM(b2)))
               ));
#endif
}

void API_CALL(z_add_contrib_target)(z_Sopalin_Data_t *sopalin_data, pastix_int_t me, pastix_int_t b1,pastix_int_t b2, pastix_int_t task, pastix_int_t t)
{
    pastix_int_t dimi,dimj,stridea,strideb,step,k;
    pastix_complex64_t *ga,*gb;
    z_SolverMatrix  *datacode    = sopalin_data->datacode;
    z_Thread_Data_t *thread_data = sopalin_data->thread_data[me];
    (void)task;

    stridea = FANIN_LROWNUM(t)-FANIN_FROWNUM(t)+1;

    /* Allocation du buffer pour la contribution a envoyer */
#ifdef ALLOC_FTGT
    MUTEX_LOCK(&(sopalin_data->mutex_fanin[t]));
#ifdef OOC_FTGT
    print_debug(DBG_OOC_FTGT, "WAIT %4d %4d\n", (int)t, (int) task);
    z_ooc_wait_for_ftgt(sopalin_data, t, me);

    ASSERTDBG(((unsigned long)(*(((double*)FANIN_COEFTAB(t))-1))) ==
              sizeof(pastix_complex64_t) * ((FANIN_LROWNUM(t)-FANIN_FROWNUM(t)+1) * (FANIN_LCOLNUM(t)-FANIN_FCOLNUM(t)+1)   *
                                        ((sopalin_data->sopar->factotype == API_FACT_LU)?2:1))
              , MOD_SOPALIN);
#else
    if (FANIN_COEFTAB(t)==NULL)
    {
        pastix_int_t j;
        pastix_int_t ftgtsize = (FANIN_LROWNUM(t)-FANIN_FROWNUM(t)+1)
            *(FANIN_LCOLNUM(t)-FANIN_FCOLNUM(t)+1);
#ifdef SOPALIN_LU
        ftgtsize *= 2;
#endif

        MALLOC_INTERN(FANIN_COEFTAB(t), ftgtsize, pastix_complex64_t);
        for (j=0;j<ftgtsize;j++)
            FANIN_COEFTAB(t)[j] = fzero;

        print_debug(DBG_SOPALIN_ALLOC, "alloc fanin coeff %x\n",(unsigned int)(intptr_t)FANIN_COEFTAB(t));

        STATS_ADD( ftgtsize );
    }
#endif /* OOC_FTGT */
    MUTEX_UNLOCK(&(sopalin_data->mutex_fanin[t]));
#endif /* ALLOC_FTGT */

    print_debug(DBG_SOPALIN_COMP1D,
                "FR(b2) %ld FR(tg) %ld\nLR(tg) %ld LR(b2) %ld\nFR(b1) %ld FC(tg) %ld\nLC(tg) %ld LR(b1) %ld\n",
                (long)SYMB_FROWNUM(b2), (long)FANIN_FROWNUM(t),
                (long)FANIN_LROWNUM(t), (long)SYMB_LROWNUM(b2),
                (long)SYMB_FROWNUM(b1), (long)FANIN_FCOLNUM(t),
                (long)FANIN_LCOLNUM(t), (long)SYMB_LROWNUM(b1));

    ASSERTDBG((SYMB_FROWNUM(b2)>=FANIN_FROWNUM(t))&&
              (FANIN_LROWNUM(t)>=SYMB_LROWNUM(b2))&&
              (SYMB_FROWNUM(b1)>=FANIN_FCOLNUM(t))&&
              (FANIN_LCOLNUM(t)>=SYMB_LROWNUM(b1)),MOD_SOPALIN);

    ga=&(FANIN_COEFTAB(t)[(SYMB_FROWNUM(b1)-FANIN_FCOLNUM(t))*stridea+
                          (SYMB_FROWNUM(b2)-FANIN_FROWNUM(t))]);

    /* vertical dimension */
    dimj=SYMB_LROWNUM(b1)-SYMB_FROWNUM(b1)+1;
    /* vertical dimension */
    dimi=SYMB_LROWNUM(b2)-SYMB_FROWNUM(b2)+1;

    strideb=thread_data->stridebloktab;

    /* Indice de debut de bloc dans maxbloktab */
    step = b2 - b1;
    for (k=b1; k<b2; k++)
        step += SYMB_LROWNUM(k) - SYMB_FROWNUM(k);

    ASSERTDBG(step+dimi<=strideb,MOD_SOPALIN);
    ASSERTDBG(strideb*dimj<=SOLV_COEFMAX,MOD_SOPALIN);

#ifdef PASTIX_ESC
    /* Décalage du bloc colonne de k + decalage du bloc k inutile */
    for (k=thread_data->firstbloktab; k<b1; k++)
        step += (SYMB_LROWNUM(k) - SYMB_FROWNUM(k) + 1) * (strideb + 1);
#endif

    gb=&(thread_data->maxbloktab2[step]);

    MUTEX_LOCK(&(sopalin_data->mutex_fanin[t]));

    SOPALIN_GESM("N","N",dimi,dimj,fun,gb,strideb,ga,stridea);

#ifdef SOPALIN_LU
    if ( (b1!=b2) &&
         (SYMB_FROWNUM(b2)>= FANIN_FCOLNUM(t)) &&
         (SYMB_LROWNUM(b2)<= FANIN_LCOLNUM(t)))
    {
        ga=&(FANIN_COEFTAB(t)[(SYMB_FROWNUM(b1)-FANIN_FCOLNUM(t))+
                              (SYMB_FROWNUM(b2)-FANIN_FROWNUM(t))*stridea]);
        gb=&(thread_data->maxbloktab1[step]);
        SOPALIN_GESM("T","N",dimj,dimi,fun,gb,strideb,ga,stridea);
    }
    else
    {
        /* ga += ftgtsize/2 */
        ga+=(FANIN_LROWNUM(t)-FANIN_FROWNUM(t)+1)
            *(FANIN_LCOLNUM(t)-FANIN_FCOLNUM(t)+1);
        gb=&(thread_data->maxbloktab1[step]);
        SOPALIN_GESM("N","N",dimi,dimj,fun,gb,strideb,ga,stridea);
    }
#endif

    /* WARNING : sauvegarder que si necessaire */
    z_ooc_save_ftgt(sopalin_data, task, t, me);

    FANIN_CTRBCNT(t)--;
    /* Add fanin to the correct ready heap */
#ifdef COMM_REORDER
    if (FANIN_CTRBCNT(t) == 0)
    {
        pastix_int_t dest = FANIN_PROCDST(t);
        MUTEX_UNLOCK(&(sopalin_data->mutex_fanin[t]));
        if (THREAD_FUNNELED_ON)
        {
            MUTEX_LOCK(&(sopalin_data->mutex_comm));
            queueAdd2(sopalin_data->sendqueue, t,
                      (double)(dest+1), (pastix_int_t)FANIN_PRIONUM(t));
            MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
        }
        else
        {
            MUTEX_LOCK(&(sopalin_data->mutex_queue_fanin[dest]));
            queueAdd2(&(sopalin_data->fanintgtsendqueue[dest]), t,
                      ((double)FANIN_PRIONUM(t)), t);
            MUTEX_UNLOCK(&(sopalin_data->mutex_queue_fanin[dest]));
        }
    }
    else
#endif
        MUTEX_UNLOCK(&(sopalin_data->mutex_fanin[t]));
}


/*
 Function: API_CALL(z_compute_1d)

 Factorisation of one block column

 Parameters:
 sopalin_data -
 me           - Thread number
 task         - z_Task number

 */
void API_CALL(z_compute_1d)(z_Sopalin_Data_t *sopalin_data, pastix_int_t me, pastix_int_t task)
{
    z_SolverMatrix  *datacode    = sopalin_data->datacode;
#if defined(TRACE_SOPALIN) || defined(PASTIX_DYNSCHED)
    z_Thread_Data_t *thread_data = sopalin_data->thread_data[me];
#endif
#ifdef PASTIX_DYNSCHED
    int            esp         = sopalin_data->sopar->iparm[IPARM_ESP];
    int            esparam     = sopalin_data->sopar->iparm[IPARM_ESP_THRESHOLD];
    pastix_int_t     t;
#endif
    pastix_int_t            c, fblknum, lblknum;
    pastix_int_t            i, ii, jj, n;
    pastix_int_t            dimb, dimb2;
#ifdef OOC
    pastix_int_t            tooc;
#endif

    c = TASK_CBLKNUM(task);

    if (sopalin_data->sopar->schur == API_YES)
    {
        pastix_int_t lN = sopalin_data->sopar->gN * sopalin_data->sopar->iparm[IPARM_DOF_NBR] - 1;
        if ( SYMB_LCOLNUM(c) == lN )
            return;
    }

    trace_begin_task(thread_data->tracefile,
                     SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                     STATE_L2_COMP1D, task);

    API_CALL(z_factor_diag)(sopalin_data, me, c);

    fblknum = SYMB_BLOKNUM(c);
    lblknum = SYMB_BLOKNUM(c + 1);

    /* Maximal M dimension for the GEMMs */
    dimb = SOLV_STRIDE(c) - ( SYMB_LROWNUM(fblknum) - SYMB_FROWNUM(fblknum) + 1 );
    fblknum++;

    /* if there is an extra-diagonal bloc in column block */
    if ( fblknum < lblknum )
    {
        API_CALL(z_factor_trsm1d)(sopalin_data, me, c);
    }
    {
        z_SolverCblk * cblk =sopalin_data->datacode->cblktab+c;
        char name[256];
        sprintf(name, "cblk_%ld_after_trf_trsm", (long)cblk->gcblknum);
        z_cblk_save(cblk, name, cblk->coeftab);
    }

    /* for all extra-diagonal column blocks */
    n = 0;
    for (i=fblknum; i<lblknum; )
    {
        ii = 1;

        /* N dimension of the GEMM */
        dimb2 = SYMB_LROWNUM(i) - SYMB_FROWNUM(i) + 1;
#ifdef PASTIX_ESC
        while( ( (i+ii) < lblknum )
               && ( SYMB_CBLKNUM(i) == SYMB_CBLKNUM(i+ii) ) ) {
            dimb2 += SYMB_LROWNUM(i+ii) - SYMB_FROWNUM(i+ii) + 1;
            ii++;
        }
#endif

#ifdef PASTIX_DYNSCHED
        t  = SOLV_INDTAB[TASK_INDNUM(task)+n];
#ifdef ESP_A
        if (esp && (dimb2*dimb2 > esparam) )
#else
            if (esp && (dimb2*dimb  > esparam) )
#endif
            {
                pastix_int_t prionum;
#ifdef ESP_WRITE
                pastix_int_t tid = TASK_THREADID(-t);
#else
                pastix_int_t tid = TASK_THREADID(task);
#endif
                if (t > 0) {
                    prionum = (double)(FANIN_PRIONUM(t));
                }
                else {
                    prionum = (double)(TASK_PRIONUM(-t));
                }

                /*
                 * We take the average priority of the current task and the
                 * target task to avoid to delay all the tasks at the same
                 * moment, creating conflicts to acces the mutex
                 */
                prionum = (pastix_int_t)floor( (( 2.* prionum + TASK_PRIONUM(task)) / 3.) - 1. );

                MUTEX_LOCK(&(sopalin_data->tasktab_mutex[tid]));
                queueAdd2(&(sopalin_data->taskqueue[tid]), TASK_TASK2ESP(task), prionum, i);
                sopalin_data->tasktab_indice[tid]--;
                MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[tid]));
                pthread_cond_broadcast(&(sopalin_data->tasktab_cond[tid]));

                thread_data->esp++;

                print_debug(DBG_SOPALIN_COMP1D,
                            "COMP1D %05d Ajout bloc %05d / taskdst %05d / priorite %05d / me %d / n %03d\n",
                            (int)task, (int)i, (int)-t, (int)prionum, (int)tid, (int)n);
            }
            else
#endif /* PASTIX_DYNSCHED */
            {
                API_CALL(z_compute_1dgemm)(sopalin_data, me, task, i, i+ii);
            }

        dimb -= dimb2;
        for (jj=0; jj<ii; jj++,i++)
            n += (lblknum - i);
    }

    trace_end_task(thread_data->tracefile,
                   SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                   STATE_L2_COMP1D, task);
}

void API_CALL(z_compute_1dgemm)(z_Sopalin_Data_t *sopalin_data, pastix_int_t me, pastix_int_t task, pastix_int_t i, pastix_int_t b2)
{
    z_SolverMatrix  *datacode    = sopalin_data->datacode;
#ifdef TRACE_SOPALIN
    z_Thread_Data_t *thread_data = sopalin_data->thread_data[me];
#endif
    pastix_int_t            j, t, fblknum, lblknum, n, usediag = 0;
    pastix_int_t            c           = TASK_CBLKNUM(task);
#ifdef OOC
    pastix_int_t            tooc;
#endif

    fblknum = SYMB_BLOKNUM(c); /* Be careful, not the same fblknum than in comp1d[plus] */
    lblknum = SYMB_BLOKNUM(c + 1);

    n = i - fblknum;
    n = (n * (n - 1)) / 2;
    n = (i - fblknum - 1) * (lblknum - fblknum) - n;

    trace_begin_task(thread_data->tracefile,
                     SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                     STATE_L2_COMP1DGEMM, task);

    /* Si on arrive par une tache qui n'a pas l'info on la recalcule */
    if (b2 == -1)
    {
        b2 = i+1;
#ifdef PASTIX_ESC
        while( (b2 < lblknum)
               && (SYMB_CBLKNUM(i) == SYMB_CBLKNUM(b2)))
            b2++;
#endif
        usediag = 1;
    }

    /* Compute contributions (GEMM) */
    API_CALL(z_compute_contrib_compact)(sopalin_data, me, c, i, b2-1, usediag);

#ifdef OOC
    if ((tooc = SOLV_INDTAB[TASK_INDNUM(task)+(n)]) < 0)
    {
        z_ooc_wait_for_cblk(sopalin_data, TASK_CBLKNUM(-tooc), me);
    }
#endif /* OOC */

    /* for all following blocks in block column */
    for (; i<b2; i++)
        for (j=i; j<lblknum; j++)
        {
            t = SOLV_INDTAB[TASK_INDNUM(task)+(n++)];
            if (t < 0)
            {
                pastix_int_t b;
                /* if the contribution is local */
#ifdef NAPA_SOPALIN
                /* if the contribution really exist */
                ASSERT(-t<SOLV_TASKNBR,MOD_SOPALIN); /* pas possible dans blend */
#endif
                b = TASK_BLOKNUM(-t);
                if (TASK_TASKID(-t)==COMP_1D)
                {
                    b = SYMB_BLOKNUM(TASK_CBLKNUM(-t));
                    /* if TASK 1D -> look for bloknum */
#ifdef NAPA_SOPALIN
                    while (!(((SYMB_FROWNUM(j)>=SYMB_FROWNUM(b)) &&
                              (SYMB_LROWNUM(j)<=SYMB_LROWNUM(b))) ||
                             ((SYMB_FROWNUM(j)<=SYMB_FROWNUM(b)) &&
                              (SYMB_LROWNUM(j)>=SYMB_LROWNUM(b))) ||
                             ((SYMB_FROWNUM(j)<=SYMB_FROWNUM(b)) &&
                              (SYMB_LROWNUM(j)>=SYMB_FROWNUM(b))) ||
                             ((SYMB_FROWNUM(j)<=SYMB_LROWNUM(b)) &&
                              (SYMB_LROWNUM(j)>=SYMB_LROWNUM(b)))))
#else
                        while (!((SYMB_FROWNUM(j)>=SYMB_FROWNUM(b)) &&
                                 (SYMB_LROWNUM(j)<=SYMB_LROWNUM(b))))
#endif
                        {
                            b++;
                            ASSERTDBG(b<SYMB_BLOKNUM(TASK_CBLKNUM(-t)+1),MOD_SOPALIN);
                        }
                }

                print_debug(DBG_SOPALIN_COMPUTE,
                            "%ld add local contrib %ld %ld %ld %ld (%ld %ld) %ld\n",
                            (long)me, (long)i, (long)j, (long)c, (long)b, (long)-t,
                            (long)TASK_TASKID(-t), (long)TASK_CTRBCNT(-t));

                z_ooc_wait_for_cblk(sopalin_data, TASK_CBLKNUM(-t), me);

                API_CALL(z_add_contrib_local)(sopalin_data, me, i, j, c, b, TASK_CBLKNUM(-t));

                z_ooc_save_coef(sopalin_data, task, TASK_CBLKNUM(-t), me);

                MUTEX_LOCK(&(sopalin_data->mutex_task[-t]));
                TASK_CTRBCNT(-t)--;
                ASSERTDBG((TASK_CTRBCNT(-t) >= 0), MOD_SOPALIN);
#ifdef PASTIX_DYNSCHED
                if ((!TASK_CTRBCNT(-t)) && (sopalin_data->taskmark[-t] == -1))
                {
                    pastix_int_t iter;

#if (DBG_PASTIX_DYNSCHED > 0)
                    ASSERTDBG(sopalin_data->taskmark[-t] == -1, MOD_SOPALIN);
#endif
                    sopalin_data->taskmark[-t]++;
                    MUTEX_UNLOCK(&(sopalin_data->mutex_task[-t]));

                    iter = TASK_THREADID(-t);
                    MUTEX_LOCK(&(sopalin_data->tasktab_mutex[iter]));
                    queueAdd(&(sopalin_data->taskqueue[iter]),
                             -t, (double)(TASK_PRIONUM(-t)));
                    MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[iter]));
                    pthread_cond_broadcast(&(sopalin_data->tasktab_cond[iter]));
                }
                else
                    MUTEX_UNLOCK(&(sopalin_data->mutex_task[-t]));
#else
                MUTEX_UNLOCK(&(sopalin_data->mutex_task[-t]));
                pthread_cond_broadcast(&(sopalin_data->cond_task[-t]));
#endif
            }
            else
            {
                /* if the contribution is not local */
#ifdef NAPA_SOPALIN
                /* if the contribution really exist */
                if (t<SOLV_FTGTNBR) {
#endif

                    print_debug(DBG_SOPALIN_COMPUTE,
                                "%ld add fanin contrib %ld %ld %ld %ld %ld\n",
                                (long)me, (long)i, (long)j, (long)c, (long)t,
                                (long)FANIN_CTRBCNT(t));

                    API_CALL(z_add_contrib_target)(sopalin_data, me, i, j, task, t);

#ifdef NAPA_SOPALIN
                }
                else
                {
                    print_debug(DBG_SOPALIN_NAPA, "ILU: drop (c=%ld,b=%ld)\n", (long)i, (long)j);
                }
#endif
            }
        }

#ifdef PASTIX_DUMP_CBLK
    {
        n = i - fblknum;
        n = (n * (n - 1)) / 2;
        n = (i - fblknum - 1) * (lblknum - fblknum) - n;
        if ((t = SOLV_INDTAB[TASK_INDNUM(task)+(n)]) < 0) {
            z_SolverCblk *cblk  = sopalin_data->datacode->cblktab+c;
            z_SolverCblk *fcblk = sopalin_data->datacode->cblktab+TASK_CBLKNUM(-t)-1;
            z_SolverBlok *blok  = sopalin_data->datacode->bloktab+i-1;
            char name[256];
            sprintf(name, "cblk_%z_after_gemm_%z_%z_%z_on_%d", fcblk->gcblknum,
                    cblk->gcblknum, blok - cblk->fblokptr, fcblk->gcblknum,
                    sopalin_data->datacode->clustnum);
            z_cblk_save(fcblk, name, fcblk->coeftab);
        }
    }
#endif

#ifdef OOC
    if (tooc < 0)
    {
        z_ooc_save_coef(sopalin_data, task, TASK_CBLKNUM(-tooc), me);
    }
#endif

    trace_end_task(thread_data->tracefile,
                   SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                   STATE_L2_COMP1DGEMM, task);

    if (THREAD_FUNNELED_OFF)
    {
        pastix_int_t dest;
        for (dest=0;dest<SOLV_PROCNBR;dest++)
        {
            if (dest == SOLV_PROCNUM) continue;
            API_CALL(send_all_fanin)(sopalin_data, me, dest);
        }
    }
}

#ifdef PASTIX_WITH_STARPU
#include "starpu_zkernels.c"
#include "starpu_zupdo_kernels.c"
#endif
