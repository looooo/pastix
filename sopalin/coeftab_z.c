/**
 *
 * @file coeftab_z.c
 *
 *  PaStiX factorization routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 * @precisions normal z -> s d c
 *
 **/
#define _GNU_SOURCE
#include "common.h"
#include "bcsc.h"

void
coeftab_zdumpcblk( const SolverCblk *cblk,
                   void *array,
                   FILE *stream );

/* Section: Functions */

void
coeftab_zffbcsc( const SolverMatrix  *solvmtx,
                 const pastix_bcsc_t *bcsc,
                 pastix_int_t         itercblk )
{
    const bcsc_format_t *csccblk = bcsc->cscftab + itercblk;
    SolverCblk *solvcblk = solvmtx->cblktab + itercblk;
    SolverBlok *solvblok;
    SolverBlok *solvblok2 = (solvcblk+1)->fblokptr;
    pastix_complex64_t *lcoeftab = solvcblk->lcoeftab;
    pastix_complex64_t *ucoeftab = solvcblk->ucoeftab;
    pastix_complex64_t *Lvalues = bcsc->Lvalues;
    pastix_complex64_t *Uvalues = bcsc->Uvalues;
    pastix_int_t itercoltab, iterval, coefindx;

    for (itercoltab=0; itercoltab<csccblk->colnbr; itercoltab++)
    {
        pastix_int_t frow = csccblk->coltab[itercoltab];
        pastix_int_t lrow = csccblk->coltab[itercoltab+1];
        solvblok = solvcblk->fblokptr;

        for (iterval=frow; iterval<lrow; iterval++)
        {
            pastix_int_t rownum = bcsc->rowtab[iterval];

            /* If values in the lower part of the matrix */
            if (rownum >= (solvcblk->fcolnum+itercoltab))
            {
                while ((solvblok < solvblok2) &&
                       ((solvblok->lrownum < rownum) ||
                        (solvblok->frownum > rownum)))
                {
                    solvblok++;
                }

                if ( solvblok < solvblok2 )
                {
                    coefindx  = solvblok->coefind;
                    coefindx += rownum - solvblok->frownum;
                    coefindx += solvcblk->stride * itercoltab;

                    lcoeftab[coefindx] = Lvalues[iterval];

                    if ( (ucoeftab != NULL) &&
                         (rownum > (solvcblk->fcolnum + itercoltab)) )
                    {
#if defined(PRECISION_z) || defined(PRECISION_c)
                        if (bcsc->mtxtype == PastixHermitian)
                            ucoeftab[coefindx] = conj(Uvalues[iterval]);
                        else
#endif
                            ucoeftab[coefindx] = Uvalues[iterval];
                    }
                }
                else {
                    /* printf("ILU: csc2solv drop coeff from CSC c=%ld(%ld) l=%ld(%ld) cblk=%ld fcol=%ld lcol=%ld\n", */
                    /*        (long)datacode->cblktab[itercblk].fcolnum+ */
                    /*        (long)itercoltab,(long)itercoltab, */
                    /*        (long)CSC_ROW(cscmtx,iterval),(long)iterval, */
                    /*        (long)itercblk, */
                    /*        (long)datacode->cblktab[itercblk].fcolnum, */
                    /*        (long)datacode->cblktab[itercblk].lcolnum); */
                }
            }
        }
    }
}


/*
 * Function: z_CoefMatrix_Allocate
 *
 * Allocate matrix coefficients in coeftab and ucoeftab.
 *
 * Should be first called with me = -1 to allocated coeftab.
 * Then, should be called with me set to thread ID
 * to allocate column blocks coefficients arrays.
 *
 * Parameters
 *
 *    datacode  - solverMatrix
 *    factotype - factorization type (LU, LLT ou LDLT)
 *    me        - thread number. (-1 for first call,
 *                from main thread. >=0 to allocate column blocks
 *     assigned to each thread.)
 */
void
coeftab_zinitcblk( const SolverMatrix  *solvmtx,
                   const pastix_bcsc_t *bcsc,
                   pastix_int_t itercblk,
                   int fakefillin, int factoLU )
{
    SolverCblk *cblk = solvmtx->cblktab + itercblk;
    pastix_int_t coefnbr = cblk->stride * cblk_colnbr( cblk );
    pastix_int_t j;

    /* If not NULL, allocated to store the shur complement for exemple */
    assert( cblk->lcoeftab == NULL );

    MALLOC_INTERN( cblk->lcoeftab, coefnbr, pastix_complex64_t );
    memset( cblk->lcoeftab, 0, coefnbr * sizeof(pastix_complex64_t) );

    if ( factoLU ) {
        MALLOC_INTERN( cblk->ucoeftab, coefnbr, pastix_complex64_t );
        memset( cblk->ucoeftab, 0, coefnbr * sizeof(pastix_complex64_t) );
    }
    else {
        cblk->ucoeftab = NULL;
    }

    /**
     * Fake initializatuion of the bloc column such that:
     *       - L is filled with 1.
     *       - U is filled with 2.
     *       - The diagonal is made dominant
     */
    if ( fakefillin ) {
        pastix_complex64_t *L = cblk->lcoeftab;
        pastix_complex64_t *U = cblk->ucoeftab;

        for (j=0; j<coefnbr; j++, L++)
        {
            *L = (pastix_complex64_t)1.;
        }

        if ( factoLU ) {
            for (j=0; j<coefnbr; j++, U++)
            {
                *U = (pastix_complex64_t)2.;
            }
        }

        /* for (j=0; j<size; itercol++) */
        /* { */
        /*     /\* On s'assure que la matrice est diagonale dominante *\/ */
        /*     SOLV_COEFTAB(itercblk)[index+itercol*stride+itercol] = (pastix_complex64_t) (UPDOWN_GNODENBR*UPDOWN_GNODENBR); */
        /* } */
    }
    else
    {
        coeftab_zffbcsc( solvmtx, bcsc, itercblk );
    }

#if defined(PASTIX_DUMP_COEFTAB)
    {
        FILE *f;
        char *filename;

        asprintf( &filename, "Lcblk%05ld.txt", itercblk );
        f = fopen( filename, "w" );
        coeftab_zdumpcblk( cblk, cblk->lcoeftab, f );
        fclose( f );
        free( filename );

        if ( cblk->ucoeftab ) {
            asprintf( &filename, "Ucblk%05ld.txt", itercblk );
            f = fopen( filename, "w" );
            coeftab_zdumpcblk( cblk, cblk->ucoeftab, f );
            fclose( f );
            free( filename );
        }
    }
#endif /* defined(PASTIX_DUMP_COEFTAB) */
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
void
coeftab_zdumpcblk( const SolverCblk *cblk,
                   void *array,
                   FILE *stream )
{
    pastix_complex64_t *coeftab = (pastix_complex64_t*)array;
    SolverBlok *blok;
    pastix_int_t itercol;
    pastix_int_t iterrow;
    pastix_int_t coefindx;

    for (itercol  = cblk->fcolnum;
         itercol <= cblk->lcolnum;
         itercol++)
    {
        /* Diagonal Block */
        blok     = cblk->fblokptr;
        coefindx = blok->coefind;
        coefindx += (itercol - cblk->fcolnum) * cblk->stride;

        for (iterrow  = blok->frownum;
             iterrow <= blok->lrownum;
             iterrow++, coefindx++)
        {
            if ((cabs( coeftab[coefindx] ) > 0.) &&
                (itercol <= iterrow))
            {
#if defined(PRECISION_z) || defined(PRECISION_c)
                fprintf(stream, "%ld %ld (%13e,%13e)\n",
                        (long)itercol, (long)iterrow,
                        creal(coeftab[coefindx]), cimag(coeftab[coefindx]));
#else
                fprintf(stream, "%ld %ld %13e\n",
                        (long)itercol, (long)iterrow,
                        coeftab[coefindx]);
#endif
            }
        }

        /* Off diagonal blocks */
        blok++;
        while( blok < (cblk+1)->fblokptr )
        {
            coefindx  = blok->coefind;
            coefindx += (itercol - cblk->fcolnum) * cblk->stride;

            for (iterrow  = blok->frownum;
                 iterrow <= blok->lrownum;
                 iterrow++, coefindx++)
            {
                if (cabs( coeftab[coefindx]) > 0.)
                {
#if defined(PRECISION_z) || defined(PRECISION_c)
                    fprintf(stream, "%ld %ld (%13e,%13e)\n",
                            (long)itercol, (long)iterrow,
                            creal(coeftab[coefindx]), cimag(coeftab[coefindx]));
#else
                    fprintf(stream, "%ld %ld %13e\n",
                            (long)itercol, (long)iterrow,
                            coeftab[coefindx]);
#endif
                }
            }
            blok++;
        }
    }
}

void
coeftab_zdump( const SolverMatrix *solvmtx,
               const char   *filename )
{
    SolverCblk *cblk = solvmtx->cblktab;
    pastix_int_t itercblk;
    FILE *stream = fopen( filename, "w" );

    for (itercblk=0; itercblk<solvmtx->cblknbr; itercblk++, cblk++)
    {
        coeftab_zdumpcblk( cblk, cblk->lcoeftab, stream );
        if ( NULL != cblk->ucoeftab )
            coeftab_zdumpcblk( cblk, cblk->lcoeftab, stream );
    }

    fclose( stream );
}


/* /\* void *\/ */
/* /\* coeftab_zinit( SopalinParam    *sopar, *\/ */
/* /\*                SolverMatrix    *datacode, *\/ */
/* /\*                pthread_mutex_t *mutex, *\/ */
/* /\*                pastix_int_t     factotype, *\/ */
/* /\*                pastix_int_t     me ) *\/ */
/* /\* { *\/ */
/* /\*   pastix_int_t i; *\/ */
/* /\*   pastix_int_t itercblk, coefnbr; *\/ */

/* /\* #ifdef PASTIX_WITH_STARPU *\/ */
/* /\*   /\\* For CUDA devices we have no allocation (yet?) *\\/ *\/ */
/* /\*   if ( sopar->iparm[IPARM_STARPU] == API_YES && me >= SOLV_THRDNBR) *\/ */
/* /\*     return; *\/ */
/* /\* #endif /\\* WITH_STARPU *\\/ *\/ */
/* /\* #ifndef OOC *\/ */
/* /\*   { *\/ */
/* /\*     /\\* On ne passe pas ici en OOC *\\/ *\/ */
/* /\*     pastix_int_t bubnum  = me; *\/ */
/* /\*     pastix_int_t task; *\/ */

/* /\* #  ifdef PASTIX_DYNSCHED *\/ */
/* /\*     while (bubnum != -1) *\/ */
/* /\*       { *\/ */
/* /\*         pastix_int_t fcandnum = datacode->btree->nodetab[bubnum].fcandnum; *\/ */
/* /\*         pastix_int_t lcandnum = datacode->btree->nodetab[bubnum].lcandnum; *\/ */
/* /\*         for (i=(me-fcandnum);i < datacode->ttsknbr[bubnum]; i+=(lcandnum-fcandnum+1)) *\/ */
/* /\* #  else *\/ */
/* /\*     for (i=0; i < datacode->ttsknbr[bubnum]; i++) *\/ */
/* /\* #  endif /\\* PASTIX_DYNSCHED *\\/ *\/ */

/* /\*       { *\/ */
/* /\*         task = datacode->ttsktab[bubnum][i]; *\/ */
/* /\*         itercblk = TASK_CBLKNUM(task); *\/ */
/* /\*         coefnbr  = SOLV_STRIDE(itercblk) * (SYMB_LCOLNUM(itercblk) - *\/ */
/* /\*                                             SYMB_FCOLNUM(itercblk) + 1); *\/ */
/* /\*         if ((TASK_TASKID(task) == COMP_1D) *\/ */
/* /\*             || (TASK_TASKID(task) == DIAG)) *\/ */
/* /\*           { *\/ */
/* /\*             datacode->cblktab[itercblk].procdiag = me; *\/ */
/* /\*             if (SOLV_COEFTAB(itercblk) == NULL) *\/ */
/* /\*               { /\\* If not NULL it should be the schur *\\/ *\/ */
/* /\*                 MALLOC_INTERN(SOLV_COEFTAB(itercblk), coefnbr, pastix_complex64_t); *\/ */
/* /\*               } *\/ */
/* /\*           } *\/ */
/* /\*         else if ( SOLV_COEFIND(TASK_BLOKNUM(task)) == 0 ) *\/ */
/* /\*           { *\/ */
/* /\*             MUTEX_LOCK(mutex); *\/ */
/* /\*             datacode->cblktab[itercblk].procdiag = me; *\/ */
/* /\*             if (SOLV_COEFTAB(itercblk) == NULL) *\/ */
/* /\*               { *\/ */
/* /\*                 MALLOC_INTERN(SOLV_COEFTAB(itercblk), coefnbr, pastix_complex64_t); *\/ */
/* /\*               } *\/ */
/* /\*             MUTEX_UNLOCK(mutex); *\/ */
/* /\*           } *\/ */
/* /\*       } *\/ */

/* /\* #  ifdef PASTIX_DYNSCHED *\/ */
/* /\*         bubnum = BFATHER(datacode->btree, bubnum); *\/ */
/* /\*       } *\/ */
/* /\* #  endif /\\* PASTIX_DYNSCHED *\\/ *\/ */
/* /\*   } *\/ */
/* /\* #endif /\\* OOC *\\/ *\/ */

/* /\*   /\\* *\/ */
/* /\*    * Allocate LU coefficient arrays *\/ */
/* /\*    * We also use it to store the diagonal in LDLt factorization using esp *\/ */
/* /\*    *\\/ *\/ */
/* /\*   if ( (factotype == API_FACT_LU) /\\* LU *\\/ *\/ */
/* /\*        || ( (factotype == API_FACT_LDLT) && sopar->iparm[IPARM_ESP] ) ) *\/ */
/* /\*     { *\/ */
/* /\* #ifndef OOC *\/ */
/* /\*       { *\/ */
/* /\*         /\\* On ne passe pas ici en OOC *\\/ *\/ */
/* /\*         pastix_int_t bubnum  = me; *\/ */
/* /\*         pastix_int_t task; *\/ */

/* /\* #  ifdef PASTIX_DYNSCHED *\/ */
/* /\*         while (bubnum != -1) *\/ */
/* /\*           { *\/ */
/* /\*             pastix_int_t fcandnum = datacode->btree->nodetab[bubnum].fcandnum; *\/ */
/* /\*             pastix_int_t lcandnum = datacode->btree->nodetab[bubnum].lcandnum; *\/ */
/* /\*             for (i=(me-fcandnum); i < datacode->ttsknbr[bubnum]; i+=(lcandnum-fcandnum+1)) *\/ */
/* /\* #  else *\/ */
/* /\*         for (i=0; i < datacode->ttsknbr[bubnum]; i++) *\/ */
/* /\* #  endif /\\* PASTIX_DYNSCHED *\\/ *\/ */

/* /\*           { *\/ */
/* /\*             task = datacode->ttsktab[bubnum][i]; *\/ */
/* /\*             itercblk = TASK_CBLKNUM(task); *\/ */
/* /\*             if ( (me != datacode->cblktab[itercblk].procdiag) *\/ */
/* /\*                  || (SOLV_UCOEFTAB(itercblk) != NULL) ) *\/ */
/* /\*               { *\/ */
/* /\*                 continue; *\/ */
/* /\*               } *\/ */

/* /\*             if ( (factotype == API_FACT_LDLT) && sopar->iparm[IPARM_ESP] ) *\/ */
/* /\*               { *\/ */
/* /\*                 coefnbr  = SYMB_LCOLNUM(itercblk) - SYMB_FCOLNUM(itercblk) + 1; *\/ */
/* /\*               } *\/ */
/* /\*             else *\/ */
/* /\*               { *\/ */
/* /\*                 coefnbr  = SOLV_STRIDE(itercblk) * (SYMB_LCOLNUM(itercblk) - *\/ */
/* /\*                                                     SYMB_FCOLNUM(itercblk) + 1); *\/ */
/* /\*               } *\/ */

/* /\*             MALLOC_INTERN(SOLV_UCOEFTAB(itercblk), coefnbr, pastix_complex64_t); *\/ */
/* /\*           } *\/ */

/* /\* #  ifdef PASTIX_DYNSCHED *\/ */
/* /\*             bubnum = BFATHER(datacode->btree, bubnum); *\/ */
/* /\*           } *\/ */
/* /\* #  endif /\\* PASTIX_DYNSCHED *\\/ *\/ */
/* /\*       } *\/ */
/* /\* #endif /\\* OOC *\\/ *\/ */
/* /\*     } *\/ */
/* /\* } *\/ */

/* /\* /\\* *\/ */
/* /\*  * Function: z_CoefMatrix_Init *\/ */
/* /\*  * *\/ */
/* /\*  * Init coeftab and ucoeftab coefficients. *\/ */
/* /\*  * *\/ */
/* /\*  * Parameters: *\/ */
/* /\*  *    datacode     - solverMatrix *\/ */
/* /\*  *    barrier      - Barrier used for thread synchronisation. *\/ */
/* /\*  *    me           - Thread ID *\/ */
/* /\*  *    iparm        - Integer parameters array. *\/ */
/* /\*  *    transcsc     - vecteur transcsc *\/ */
/* /\*  *    sopalin_data - <z_Sopalin_Data_t> structure. *\/ */
/* /\*  *\\/ *\/ */
/* /\* void z_CoefMatrix_Init(z_SolverMatrix         *datacode, *\/ */
/* /\*                        sopthread_barrier_t  *barrier, *\/ */
/* /\*                        pastix_int_t                   me, *\/ */
/* /\*                        pastix_int_t                  *iparm, *\/ */
/* /\*                        pastix_complex64_t               **transcsc, *\/ */
/* /\*                        z_Sopalin_Data_t       *sopalin_data) *\/ */
/* /\* { *\/ */

/* /\*     pastix_int_t j, itercblk; *\/ */
/* /\*     pastix_int_t i, coefnbr; *\/ */

/* /\* #ifdef PASTIX_WITH_STARPU *\/ */
/* /\*     /\\* For CUDA devices we have no allocation (yet?) *\\/ *\/ */
/* /\*     if ( iparm[IPARM_STARPU] == API_YES && me >= SOLV_THRDNBR) *\/ */
/* /\*         return; *\/ */
/* /\* #endif /\\* WITH_STARPU *\\/ *\/ */

/* /\*     /\\* Remplissage de la matrice *\\/ *\/ */
/* /\*     if (iparm[IPARM_FILL_MATRIX] == API_NO) *\/ */
/* /\*     { *\/ */
/* /\*         /\\* Remplissage par bloc *\\/ *\/ */
/* /\*         pastix_int_t bubnum  = me; *\/ */
/* /\* #ifdef PASTIX_DYNSCHED *\/ */
/* /\*         while (bubnum != -1) *\/ */
/* /\*         { *\/ */
/* /\*             pastix_int_t fcandnum = datacode->btree->nodetab[bubnum].fcandnum; *\/ */
/* /\*             pastix_int_t lcandnum = datacode->btree->nodetab[bubnum].lcandnum; *\/ */
/* /\*             for (i=(me-fcandnum);i < datacode->ttsknbr[bubnum]; i+=(lcandnum-fcandnum+1)) *\/ */
/* /\* #else *\/ */
/* /\*                 for (i=0; i < datacode->ttsknbr[bubnum]; i++) *\/ */
/* /\* #endif /\\* PASTIX_DYNSCHED *\\/ *\/ */

/* /\*                 { *\/ */
/* /\*                     pastix_int_t task; *\/ */
/* /\*                     pastix_int_t k = i; *\/ */
/* /\* #ifdef OOC *\/ */
/* /\*                     /\\* En OOC, on inverse la boucle pour conserver les premiers blocs en mémoire *\\/ *\/ */
/* /\*                     k = datacode->ttsknbr[bubnum]-i-1; *\/ */
/* /\* #endif *\/ */
/* /\*                     task = datacode->ttsktab[bubnum][k]; *\/ */
/* /\*                     itercblk = TASK_CBLKNUM(task); *\/ */

/* /\*                     if ( me != datacode->cblktab[itercblk].procdiag ) *\/ */
/* /\*                         continue; *\/ */

/* /\*                     coefnbr  = SOLV_STRIDE(itercblk) * (SYMB_LCOLNUM(itercblk) - SYMB_FCOLNUM(itercblk) + 1); *\/ */

/* /\*                     z_ooc_wait_for_cblk(sopalin_data, itercblk, me); *\/ */

/* /\*                     /\\* initialisation du bloc colonne *\\/ *\/ */
/* /\*                     for (j=0 ; j < coefnbr ; j++) *\/ */
/* /\*                     { *\/ */
/* /\*                         SOLV_COEFTAB(itercblk)[j] = ZERO; *\/ */
/* /\*                         if (iparm[IPARM_FACTORIZATION] == API_FACT_LU) *\/ */
/* /\*                             SOLV_UCOEFTAB(itercblk)[j] = ZERO; *\/ */
/* /\*                     } *\/ */

/* /\*                     /\\* remplissage *\\/ *\/ */
/* /\*                     z_Csc2solv_cblk(sopalin_data->sopar->cscmtx, datacode, *transcsc, itercblk); *\/ */

/* /\*                     z_ooc_save_coef(sopalin_data, task, itercblk, me); *\/ */
/* /\*                 } *\/ */

/* /\* #ifdef PASTIX_DYNSCHED *\/ */
/* /\*             bubnum = BFATHER(datacode->btree, bubnum); *\/ */
/* /\*         } *\/ */
/* /\* #endif /\\* PASTIX_DYNSCHED *\\/ *\/ */

/* /\* #ifdef DEBUG_COEFINIT *\/ */
/* /\*         if (me == 0) *\/ */
/* /\*         { *\/ */
/* /\*             FILE *transfile; *\/ */
/* /\*             char transfilename[10]; *\/ */
/* /\*             sprintf(transfilename, "trans%ld.%ld",(long) me,(long) SOLV_PROCNUM); *\/ */
/* /\*             transfile = fopen(transfilename, "w"); *\/ */
/* /\*             z_dump7(*transcsc, transfile); *\/ */
/* /\*         } *\/ */
/* /\* #endif *\/ */
/* /\*         /\\* Libération de mémoire *\\/ *\/ */
/* /\* #ifdef STARPU_INIT_SMP *\/ */
/* /\*         if (sopalin_data->sopar->iparm[IPARM_STARPU] == API_NO) *\/ */
/* /\* #endif /\\* STARPU_INIT_SMP *\\/ *\/ */
/* /\*             SYNCHRO_X_THREAD(SOLV_THRDNBR, *barrier); *\/ */
/* /\*         if (me == 0) *\/ */
/* /\*         { *\/ */
/* /\*             if (iparm[IPARM_FACTORIZATION] != API_FACT_LU) *\/ */
/* /\*             { *\/ */
/* /\*                 if (*transcsc != NULL) *\/ */
/* /\*                     memFree_null(*transcsc); *\/ */
/* /\*             } *\/ */
/* /\*             else *\/ */
/* /\*             { *\/ */
/* /\*                 if (iparm[IPARM_SYM] == API_SYM_YES || iparm[IPARM_SYM] == API_SYM_HER) /\\* Symmetric *\\/ *\/ */
/* /\*                     *transcsc = NULL; *\/ */
/* /\*                 else /\\* Unsymmetric *\\/ *\/ */
/* /\*                     memFree_null(*transcsc); *\/ */
/* /\*             } *\/ */
/* /\*         } *\/ */
/* /\*     } *\/ */
/* /\*     else  /\\* fake factorisation *\\/ *\/ */
/* /\*     { *\/ */

/* /\*         /\\* deadcode *\\/ *\/ */
/* /\*         pastix_int_t itercol; *\/ */

/* /\*         /\\* Initialisation de la matrice à 0 et 1 ou 2 *\\/ *\/ */
/* /\*         pastix_int_t task; *\/ */
/* /\*         pastix_int_t bubnum  = me; *\/ */

/* /\* #ifdef PASTIX_DYNSCHED *\/ */
/* /\*         while (bubnum != -1) *\/ */
/* /\*         { *\/ */
/* /\*             pastix_int_t fcandnum = datacode->btree->nodetab[bubnum].fcandnum; *\/ */
/* /\*             pastix_int_t lcandnum = datacode->btree->nodetab[bubnum].lcandnum; *\/ */
/* /\*             for (i=(me-fcandnum);i < datacode->ttsknbr[bubnum]; i+=(lcandnum-fcandnum+1)) *\/ */
/* /\* #else *\/ */
/* /\*                 for (i=0; i < datacode->ttsknbr[bubnum]; i++) *\/ */
/* /\* #endif /\\* PASTIX_DYNSCHED *\\/ *\/ */

/* /\*                 { *\/ */
/* /\*                     task = datacode->ttsktab[bubnum][i]; *\/ */
/* /\*                     itercblk = TASK_CBLKNUM(task); *\/ */
/* /\*                     if ( me != datacode->cblktab[itercblk].procdiag ) *\/ */
/* /\*                         continue; *\/ */
/* /\*                     coefnbr  = SOLV_STRIDE(itercblk) * (SYMB_LCOLNUM(itercblk) - SYMB_FCOLNUM(itercblk) + 1); *\/ */

/* /\*                     for (j=0 ; j < coefnbr ; j++) *\/ */
/* /\*                     { *\/ */
/* /\*                         if (iparm[IPARM_FILL_MATRIX] == API_NO) *\/ */
/* /\*                             SOLV_COEFTAB(itercblk)[j] = ZERO; *\/ */
/* /\*                         else *\/ */
/* /\*                             SOLV_COEFTAB(itercblk)[j] = UN; *\/ */
/* /\*                         if (iparm[IPARM_FACTORIZATION] == API_FACT_LU) *\/ */
/* /\*                         { *\/ */
/* /\*                             if (iparm[IPARM_FILL_MATRIX] == API_NO) *\/ */
/* /\*                                 SOLV_UCOEFTAB(itercblk)[j] = ZERO; *\/ */
/* /\*                             else *\/ */
/* /\*                                 SOLV_UCOEFTAB(itercblk)[j] = DEUX; *\/ */
/* /\*                         } *\/ */
/* /\*                     } *\/ */
/* /\* #ifdef _UNUSED_ *\/ */
/* /\*                 } *\/ */
/* /\* #endif *\/ */
/* /\*         } *\/ */

/* /\* #ifdef PASTIX_DYNSCHED *\/ */
/* /\*         bubnum = BFATHER(datacode->btree, bubnum); *\/ */
/* /\*     } *\/ */
/* /\* #endif /\\* PASTIX_DYNSCHED *\\/ *\/ */

/* /\*     /\\* 2 eme phase de l'initialisation de la matrice *\\/ *\/ */
/* /\*     for (i=0; i < SOLV_TTSKNBR; i++) *\/ */
/* /\*     { *\/ */
/* /\*         itercblk = TASK_CBLKNUM(SOLV_TTSKTAB(i)); *\/ */
/* /\*         coefnbr  = SOLV_STRIDE(itercblk) * (SYMB_LCOLNUM(itercblk) - SYMB_FCOLNUM(itercblk) + 1); *\/ */

/* /\*         z_ooc_wait_for_cblk(sopalin_data, itercblk, me); *\/ */

/* /\*         /\\* initialisation du bloc colonne *\\/ *\/ */
/* /\*         for (j=0 ; j < coefnbr ; j++) *\/ */
/* /\*         { *\/ */
/* /\*             SOLV_COEFTAB(itercblk)[j] = UN; *\/ */
/* /\*             if (iparm[IPARM_FACTORIZATION] == API_FACT_LU) *\/ */
/* /\*                 SOLV_UCOEFTAB(itercblk)[j] = DEUX; *\/ */
/* /\*         } *\/ */

/* /\*         /\\* if we are on a diagonal bloc *\\/ *\/ */
/* /\*         if (SYMB_FCOLNUM(itercblk) == SYMB_FROWNUM(SYMB_BLOKNUM(itercblk))) *\/ */
/* /\*         { *\/ */
/* /\*             pastix_int_t index  = SOLV_COEFIND(SYMB_BLOKNUM(itercblk)); *\/ */
/* /\*             pastix_int_t size   = SYMB_LCOLNUM(itercblk) - SYMB_FCOLNUM(itercblk) + 1; *\/ */
/* /\*             pastix_int_t stride = SOLV_STRIDE(itercblk); *\/ */

/* /\*             for (itercol=0; itercol<size; itercol++) *\/ */
/* /\*             { *\/ */
/* /\*                 /\\* On s'assure que la matrice est diagonale dominante *\\/ *\/ */
/* /\*                 SOLV_COEFTAB(itercblk)[index+itercol*stride+itercol] = (pastix_complex64_t) (UPDOWN_GNODENBR*UPDOWN_GNODENBR); *\/ */
/* /\*             } *\/ */
/* /\*             /\\* copie de la partie de block diag de U dans L *\\/ *\/ */
/* /\*             if (iparm[IPARM_FACTORIZATION] == API_FACT_LU) *\/ */
/* /\*             { *\/ */
/* /\*                 pastix_int_t iterrow; *\/ */
/* /\*                 for (itercol=0; itercol<size; itercol++) *\/ */
/* /\*                 { *\/ */
/* /\*                     for (iterrow=itercol+1; iterrow<size; iterrow++) *\/ */
/* /\*                     { *\/ */
/* /\*                         SOLV_COEFTAB(itercblk)[index+iterrow*stride+itercol] = SOLV_UCOEFTAB(itercblk)[index+itercol*stride+iterrow]; *\/ */
/* /\*                     } *\/ */
/* /\*                 } *\/ */
/* /\*             } *\/ */
/* /\*         } *\/ */

/* /\*         z_ooc_save_coef(sopalin_data, SOLV_TTSKTAB(i), itercblk, me); *\/ */
/* /\*     } *\/ */
/* /\* #ifdef _UNUSED_ *\/ */
/* /\* } *\/ */
/* /\* #endif *\/ */
/* /\* printf("fin false fill-in\n"); *\/ */
/* /\* } *\/ */


/* /\* if (iparm[IPARM_FREE_CSCPASTIX] == API_CSC_FREE) *\/ */
/* /\*  { *\/ */
/* /\*    SYNCHRO_X_THREAD(SOLV_THRDNBR, *barrier); *\/ */

/* /\*    /\\* Internal csc is useless if we don't want to do refinement step *\\/ *\/ */
/* /\*    if (me == 0) *\/ */
/* /\*    { *\/ */
/* /\*      if ((iparm[IPARM_END_TASK] < API_TASK_SOLVE) || *\/ */
/* /\*          ((iparm[IPARM_END_TASK] < API_TASK_REFINE) && *\/ */
/* /\*           (iparm[IPARM_RHS_MAKING] == API_RHS_B))) *\/ */
/* /\*      { *\/ */
/* /\*        z_CscExit(sopalin_data->sopar->cscmtx); *\/ */
/* /\*      } *\/ */
/* /\*      else *\/ */
/* /\*      { *\/ */
/* /\*        errorPrintW("The internal CSC can't be freed if you want to use refinement or if you don't give one RHS.\n"); *\/ */
/* /\*      } *\/ */
/* /\*    } *\/ */
/* /\*  } *\/ */
/* /\* } *\/ */

/* /\* /\\* *\/ */
/* /\*  Function: z_CoefMatrix_Free *\/ */

/* /\*  Free the z_solver matrix coefficient tabular : coeftab and ucoeftab. *\/ */

/* /\*  WARNING: Call it with one unnique thread. *\/ */

/* /\*  Parameters: *\/ */
/* /\*  datacode   - solverMatrix *\/ */
/* /\*  factotype  - factorisation type (<API_FACT>) *\/ */

/* /\*  *\\/ *\/ */
/* /\* void z_CoefMatrix_Free(z_SopalinParam *sopar, *\/ */
/* /\*                      z_SolverMatrix *datacode, *\/ */
/* /\*                      pastix_int_t           factotype) *\/ */
/* /\* { *\/ */
/* /\*   pastix_int_t i; *\/ */

/* /\*   if ( (factotype == API_FACT_LU) *\/ */
/* /\*        || ( (factotype == API_FACT_LDLT) && sopar->iparm[IPARM_ESP]) ) *\/ */
/* /\*   { *\/ */
/* /\*     for (i=0 ; i < SYMB_CBLKNBR; i++) *\/ */
/* /\*       if (SOLV_UCOEFTAB(i) != NULL) *\/ */
/* /\*         memFree_null(SOLV_UCOEFTAB(i)); *\/ */
/* /\*   } *\/ */
/* /\*   for (i=0 ; i < SYMB_CBLKNBR; i++) *\/ */
/* /\*     if (SOLV_COEFTAB(i) != NULL) *\/ */
/* /\*       memFree_null(SOLV_COEFTAB(i)); *\/ */
/* /\* } *\/ */
