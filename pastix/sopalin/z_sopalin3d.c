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
 File: sopalin3d.c

 sopalin 3d main program.

 Authors:
 Mathieu Faverge - faverge@labri.fr
 Xavier  Lacoste - lacoste@labri.fr
 Pierre  Ramet   - ramet@labri.fr

 Date:
 Version 0.0 - february 2003
 */

/***********************************************/
/*                HEADERS                      */
/***********************************************/

/*
 * Redefinition du test pour madeleine, car les
 * requetes madeleine sont des structures
 */

#ifndef MAD_MPI
#define MPI_Request_is_equal(r1, r2) ((r1) == (r2))
#endif

/* Utilisation des blas IBM sur power */

#if (defined X_ARCHpower_ibm_aix)
#if (!(defined X_INCLUDE_ESSL))
#define X_INCLUDE_ESSL
#endif
#endif

/* Include pour compaq */
#if (defined X_ARCHalpha_compaq_osf1)
#include <sys/types.h>
#include <sys/resource.h>
#include <sys/processor.h>
#include <sys/sysinfo.h>
#include <machine/hal_sysinfo.h>
#define X_INCLUDE_CXML
#endif

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <pthread.h>

#ifdef X_INCLUDE_ESSL
#include <essl.h>
/*#include <pessl.h>*/
#endif /* X_INCLUDE_ESSL */

#include <signal.h>
#include "common.h"
#include "z_tools.h"
#ifdef PASTIX_EZTRACE
#include "pastix_eztrace.h"
#else
#include "trace.h"
#endif
#include "sopalin_define.h"
#include "symbol.h"
#include "z_ftgt.h"
#include "z_csc.h"
#include "z_updown.h"
#include "queue.h"
#include "bulles.h"
#include "z_solver.h"
#include "sopalin_thread.h"
#include "stack.h"
#include "z_sopalin3d.h"
#include "z_sopalin_init.h"
#include "perf.h"
#include "out.h"
#include "z_coefinit.h"
#include "z_ooc.h"
#include "order.h"
#include "z_debug_dump.h"
#ifdef PASTIX_WITH_STARPU
#include "starpu_zsubmit_tasks.h"
#endif

/**************************************/
/*          A Lire                    */
/**************************************/

/* ATTENTION cast alignement des buffers de reception (facto+updo) complexe */
/* SMP : mutex_queue necessaire ??? (si les files sont par thread) */
/* SMP : PB_SMP commenter SYMB_CBLKNUM(-t)=... 2D sans SOLVER_UPDOWN ??? */
/* SMP : PB_SMP free blocktarget ??? */
/* EXACT_TAG risque debordement SEPFB -> EXACT_THREAD pb dans updown ??? */


/**************************************/
/*           Options de compil        */
/**************************************/


/* Variables pour la version SMP */
#define SYNCHRO_THREAD  SYNCHRO_X_THREAD(SOLV_THRDNBR, sopalin_data->barrier)
#define MAXTHRDS        SOLV_THRDNBR

/* Définition d'un décalage pour les tags fanin et block */
#define SEPFB MAXTHRDS
#if (defined EXACT_TAG)
#undef  SEPFB
#define SEPFB 40000
#endif


/*static pastix_int_t iterator;*/
#define DO_ITER_MAX 10
/*#define DO_ITER(x) {for(iterator=0;iterator<DO_ITER_MAX;iterator++){(x);}}*/
#define DO_ITER(x) {if (SOLV_PROCNBR > 1) {(x);};}


#include "sopalin_time.h"

#ifdef COMPUTE

#include "z_sopalin_compute.h"



#else /* COMPUTE */

/* BLAS computations are disabled */
#define SOPALIN_GEAM
#define SOPALIN_GESM
#define SOPALIN_COPY
#define SOPALIN_SCAL
#define SOPALIN_PPF
#define SOPALIN_POF
#define SOPALIN_TRSM
#define SOPALIN_TRSV
#define SOPALIN_GEMM
#define SOPALIN_GEMV
#define SOPALIN_AXPY
#define SOPALIN_GER
#define SOPALIN_SYR

#endif /* COMPUTE */

/* Definition des debug */
#include "sopalin_debug.h"

/* Acces aux champs de datacode */
#include "sopalin_acces.h"

/***********************************/
/*      Affichage                  */
/***********************************/
/*
 Section: Macros

 Macros: Printing maccros.

 print_onempi(fmt, ...) - Print message by one MPI task.
 print_one(fmt, ...)    - Print message by one thread of one MPI task
 print_error(...)       - Print error messages (ignored).

 */
#define print_onempi(fmt, ...) if( SOLV_PROCNUM == 0 )           fprintf(stdout, fmt, __VA_ARGS__)
#define print_one(fmt, ...)    if( me == 0 && SOLV_PROCNUM == 0) fprintf(stdout, fmt, __VA_ARGS__)
#define print_error(...)


#ifdef VERIF_MPI
int err_mpi;
#endif

#define LOCAL_ALLOC_BTAG -100

/************************************************/
/*       Déclaration des fonctions              */
/************************************************/
/* Section: Prototypes */
#define z_dump_all                  API_CALL(z_dump_all)
#define z_init_struct_sopalin       API_CALL(z_init_struct_sopalin)
#define z_sopalin_launch            API_CALL(z_sopalin_launch)
#define z_sopalin_updo_comm         API_CALL(z_sopalin_updo_comm)
#define z_sopalin_thread            API_CALL(z_sopalin_thread)
#define z_sopalin_smp               API_CALL(z_sopalin_smp)
#define z_sopalin_updo_thread       API_CALL(z_sopalin_updo_thread)
#define z_sopalin_updo_smp          API_CALL(z_sopalin_updo_smp)
#define z_sopalin_updo_gmres_thread API_CALL(z_sopalin_updo_gmres_thread)
#define z_sopalin_updo_gmres_smp    API_CALL(z_sopalin_updo_gmres_smp)
#define z_sopalin_updo_grad_thread  API_CALL(z_sopalin_updo_grad_thread)
#define z_sopalin_updo_grad_smp     API_CALL(z_sopalin_updo_grad_smp)
#define z_sopalin_updo_pivot_thread API_CALL(z_sopalin_updo_pivot_thread)
#define z_sopalin_updo_pivot_smp    API_CALL(z_sopalin_updo_pivot_smp)
#define z_up_down                   API_CALL(z_up_down)
#define z_up_down_smp               API_CALL(z_up_down_smp)



void  z_dump_all                 (z_SolverMatrix *, z_CscMatrix * cscmtx, int);
void  z_init_struct_sopalin      (z_Sopalin_Data_t *sopalin_data, z_SolverMatrix *m,
                                z_SopalinParam *sopar);
void  z_sopalin_launch           (z_SolverMatrix *m, z_SopalinParam *sopaparam, pastix_int_t cas);
void* z_sopalin_updo_comm        (void *arg);
void  z_sopalin_thread           (z_SolverMatrix *m, z_SopalinParam *sopaparam);
void* z_sopalin_smp              (void *arg);
void  z_sopalin_updo_thread      (z_SolverMatrix *m, z_SopalinParam *sopaparam);
void* z_sopalin_updo_smp         (void *arg);
void  z_sopalin_updo_gmres_thread(z_SolverMatrix *m, z_SopalinParam *sopaparam);
void* z_sopalin_updo_gmres_smp   (void *arg);
void  z_sopalin_updo_grad_thread (z_SolverMatrix *m, z_SopalinParam *sopaparam);
void* z_sopalin_updo_grad_smp    (void *arg);
void  z_sopalin_updo_pivot_thread(z_SolverMatrix *m, z_SopalinParam *sopaparam);
void* z_sopalin_updo_pivot_smp   (void *arg);
void  z_up_down                  (void);
void* z_up_down_smp              (void * arg);

/************************************************/
/*              Dump des matrices               */
/************************************************/
/* Section: Debug functions */
/*
 Function: API_CALL(z_dump_all)

 Dumps the matrix and right-hand-side on disk.

 This function can dump the internal distributed CSC matrix,
 or the solvermatrix, or the Up-down vector.

 This function must be called by only one thread for
 each MPI process.

 The *x* value can be defined using *DUMP_CSC*, *DUMP_SOLV* and *DUMP_SMB*,
 combined like *DUMP_CSC | DUMP_SOLV | DUMP_SMB*

 Parameters:
 datacode - z_SolverMatrix
 x        - value indicating what to dump.

 Returns:
 void
 */
void z_dump_all(z_SolverMatrix *datacode,
              z_CscMatrix    *cscmtx,
              int           x)
{
    /*
     ** Fonction appellée par un seul thread par processus MPI
     **   instance : étape où sont dumpées les matrices
     **   x        : vecteurs a dumper.
     */
    static int instance = 0;
    FILE *file;
    char  filename[250];
#ifdef SOPALIN_LU
    FILE *fileu;
    char  filenameu[250];
#endif

    printf("Dump CSC SOLV SMB (%ld)...\n",(long) instance);

    /* CSC */
#ifdef USE_CSC
    if (x & DUMP_CSC)
    {
        sprintf(filename, "csc%ld.%ld",
                (long) instance,
                (long) SOLV_PROCNUM);
        file = fopen(filename, "w");
        z_dump2(datacode, cscmtx, NULL, file);
        fclose(file);
    }
#endif

    /* SOLV */
    if (x & DUMP_SOLV)
    {
        sprintf(filename, "solv%ld.%ld",(long) instance,(long) SOLV_PROCNUM);
        file = fopen(filename, "w");
#ifdef SOPALIN_LU
        sprintf(filenameu, "solvu%ld.%ld",(long) instance,(long) SOLV_PROCNUM);
        fileu = fopen(filenameu, "w");
        z_dump3_LU(datacode, file, fileu);
        fclose(fileu);
#else
        z_dump3(datacode, file);
#endif
        fclose(file);
    }

    /* SMB */
    if (x & DUMP_SMB)
    {
        sprintf(filename, "smb%ld.%ld",(long) instance,(long) SOLV_PROCNUM);
        file = fopen( filename, "w");
        z_dump5(datacode, file);
        fclose(file);
    }
    instance++;
}

/****************************************************************************/
/*                     COMMUNICATION ROUTINES                               */
/****************************************************************************/

#include "z_sopalin_sendrecv.c"

/****************************************************************************/
/*                       COMPUTE ROUTINES                                   */
/****************************************************************************/

/* Section : Computation routines prototypes */
void API_CALL(z_compute_diag)  (z_Sopalin_Data_t *sopalin_data, pastix_int_t me, pastix_int_t task);
void API_CALL(z_compute_1d)    (z_Sopalin_Data_t *sopalin_data, pastix_int_t me, pastix_int_t task);
void API_CALL(z_compute_1dgemm)(z_Sopalin_Data_t *sopalin_data, pastix_int_t me, pastix_int_t task, pastix_int_t i, pastix_int_t b2);
void API_CALL(z_compute_e1)    (z_Sopalin_Data_t *sopalin_data, pastix_int_t me, pastix_int_t task);
void API_CALL(z_compute_e2)    (z_Sopalin_Data_t *sopalin_data, pastix_int_t me, pastix_int_t task);

#include "z_sopalin_compute.c"

#if defined(PASTIX_DYNSCHED)
#include "z_dynsched.h"
#endif
#ifdef SOLVER_UPDOWN
/****************************************************************************/
/*                       UP AND DOWN STEP                                   */
/****************************************************************************/
/* Section : Up-down step routines prototypes */
/* Initialisation / Nettoyage */
void  API_CALL(z_updo_init)(z_Sopalin_Data_t *sopalin_data, z_SolverMatrix *datacode, z_SopalinParam *sopaparam);

/* Thread de communication */
void* API_CALL(z_updo_thread_comm)(void *);

/* Lancement de updo seul */
void  API_CALL(z_updo_thread)(z_SolverMatrix *datacode, z_SopalinParam *sopaparam);

#include "z_updo.c"
/* Section : Reffinement step routines prototypes */
/* Raffinement du second membre */
void* API_CALL(z_pivotstatique_smp)(void *arg);
void* API_CALL(z_gmres_smp)        (void *arg);
void* API_CALL(z_grad_smp)         (void *arg);
void* API_CALL(z_bicgstab_smp)      (void *arg);

/* Lancement d'une des fonctions seules */
void  API_CALL(z_pivot_thread)   (z_SolverMatrix *datacode, z_SopalinParam *sopaparam);
void  API_CALL(z_gmres_thread)   (z_SolverMatrix *datacode, z_SopalinParam *sopaparam);
void  API_CALL(z_grad_thread)    (z_SolverMatrix *datacode, z_SopalinParam *sopaparam);
void  API_CALL(z_bicgstab_thread)(z_SolverMatrix *datacode, z_SopalinParam *sopaparam);

#include "z_csc_intern_compute.h"

#define RAFF_CLOCK_INIT {clockInit(raff_clk);clockStart(raff_clk);}
#define RAFF_CLOCK_STOP {clockStop((raff_clk));}
#define RAFF_CLOCK_GET  clockGet()

/* #define DEBUG_RAFF */
#include "z_raff_functions.h"
#include "z_raff_grad.c"
#include "z_raff_gmres.c"
#include "z_raff_pivot.c"
#include "z_raff_bicgstab.c"

#endif

/****************************************************************************/
/*                  INITIALIZATION ROUTINES                                 */
/****************************************************************************/

/* Initialisation de la CSC */
/* Section: Sopalin functions */

/*
 * Function: z_init_struct_sopalin
 *
 * Initialization routine.
 *
 * Set the sopalin_data->critere, allocate FANIN_COEFTAB.
 *
 * Parameters:
 *       sopalin_data - Structure used during factorisation and resolution.
 *       datacode     - z_SolverMatrix structure.
 *	sopar        - Factorisation parameters.
 */
#define z_init_struct_sopalin API_CALL(z_init_struct_sopalin)
void z_init_struct_sopalin (z_Sopalin_Data_t *sopalin_data,
                          z_SolverMatrix   *datacode,
                          z_SopalinParam   *sopar)
{
    MPI_Comm pastix_comm = PASTIX_COMM;
    pastix_int_t      i;
#if (defined DEADCODE) && !(defined ALLOC_FTGT)
    pastix_int_t j;
#endif

    if (sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        print_onempi("%s", OUT2_SOP_TABG);

    /* A la louche ??? */
#ifdef SOPALIN_LU
    PACKAREA     *= 2;
#endif

#if (defined COMM_REORDER) && (defined PASTIX_NMAD_AGGREG)
    PACKMAX = 1;
#endif

    /*PACKMAX=(PACKMAX<32)?PACKMAX:32;*/

    print_debug(DBG_SOPALIN_MAIN,"--- PACKMAX=%ld PACKAREA=%ld\n", (long) PACKMAX, (long) PACKAREA);
    print_debug(DBG_SOPALIN_MAIN,"--- ** sopalin coefnbr=%ld ftgtnbr=%ld coefmax=%ld cblknbr=%ld bloknbr=%ld procnum=%ld procnbr=%ld\n",
                (long) SOLV_COEFNBR, (long) SOLV_FTGTNBR,
                (long) SOLV_COEFMAX, (long) SYMB_CBLKNBR,
                (long) SYMB_BLOKNBR, (long) SOLV_PROCNUM, (long) SOLV_PROCNBR);

    /* Redéfinition des erreurs pour Compaq */
#ifdef X_INCLUDE_ESSL
    {
        int ierno,inoal,inomes,itrace,iusadr,irange,dummy;

        /* Special pour le pere Goudin ... */

#define NUM_ERROR_PPF 2104
#define NUM_ERROR_POF 2115
#define NOT_ALTERED      0
#define LIMIT_ERROR    255
#define UNLIMITED_ERROR LIMIT_ERROR+1
#define IGNORE_ERROR    -1
#define IGNORE_MESSAGE  -1

        /*
         printf("les modifs du pere goudin\n");

         einfo(0,&dummy,&dummy);

         ierno=NUM_ERROR_PPF;
         inoal=UNLIMITED_ERROR;
         inomes=NOT_ALTERED;
         itrace=NOT_ALTERED;
         iusadr=NOT_ALTERED;
         irange=NUM_ERROR_PPF;

         errset(&ierno,&inoal,&inomes,&itrace,&iusadr,&irange);

         ierno=NUM_ERROR_POF;
         inoal=UNLIMITED_ERROR;
         inomes=NOT_ALTERED;
         itrace=NOT_ALTERED;
         iusadr=NOT_ALTERED;
         irange=NUM_ERROR_POF;

         errset(&ierno,&inoal,&inomes,&itrace,&iusadr,&irange);
         */
    }
#endif /* X_INCLUDE_ESSL */

    /* Statistiques d'allocation */
#ifdef ALLOC_FTGT
    {
        double factor = 0.0;
        pastix_int_t    alloc_init =
            SYMB_CBLKNBR*3*sizeof(pastix_int_t)+
            SYMB_BLOKNBR*3*sizeof(pastix_int_t)+
            SYMB_CBLKNBR*1*sizeof(pastix_int_t)+
            SYMB_BLOKNBR*1*sizeof(pastix_int_t)+
            SOLV_TASKNBR  *sizeof(z_Task)+
            SOLV_FTGTNBR  *sizeof(z_FanInTarget)+
            SOLV_COEFNBR  *sizeof(pastix_complex64_t)+
            SOLV_INDNBR   *sizeof(pastix_int_t);

        factor = 100.0 / (double)alloc_init;
        (void)factor;
        if (sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_CHATTERBOX) {
            printf("symbol.cblk %12ld %2.2lf %%\n",
                   (long)  (SYMB_CBLKNBR*3*sizeof(pastix_int_t)),
                   (double)(SYMB_CBLKNBR*3*sizeof(pastix_int_t))*factor);
            printf("symbol.blok %12ld %2.2lf %%\n",
                   (long)  (SYMB_BLOKNBR*3*sizeof(pastix_int_t)),
                   (double)(SYMB_BLOKNBR*3*sizeof(pastix_int_t))*factor);
            printf("z_solver.cblk %12ld %2.2lf %%\n",
                   (long)  (SYMB_CBLKNBR*1*sizeof(pastix_int_t)),
                   (double)(SYMB_CBLKNBR*1*sizeof(pastix_int_t))*factor);
            printf("z_solver.blok %12ld %2.2lf %%\n",
                   (long)  (SYMB_BLOKNBR*1*sizeof(pastix_int_t)),
                   (double)(SYMB_BLOKNBR*1*sizeof(pastix_int_t))*factor);
            printf("z_solver.task %12ld %2.2lf %%\n",
                   (long)  (SOLV_TASKNBR*1*sizeof(z_Task)),
                   (double)(SOLV_TASKNBR*1*sizeof(z_Task))*factor);
            printf("z_solver.ftgt %12ld %2.2lf %%\n",
                   (long)  (SOLV_FTGTNBR*1*sizeof(z_FanInTarget)),
                   (double)(SOLV_FTGTNBR*1*sizeof(z_FanInTarget))*factor);
            printf("z_solver.coef %12ld %2.2lf %%\n",
                   (long)  (SOLV_COEFNBR*1*sizeof(pastix_complex64_t)),
                   (double)(SOLV_COEFNBR*1*sizeof(pastix_complex64_t))*factor);
            printf("z_solver.ind  %12ld %2.2lf %%\n",
                   (long)  (SOLV_INDNBR *1*sizeof(pastix_int_t)),
                   (double)(SOLV_INDNBR *1*sizeof(pastix_int_t))*factor);
        }
    }
#endif

    /* Computing sopalin_data->critere */
    /* Computing sopalin_data->critere */
#ifdef COMPUTE
#ifdef USE_CSC

    sopalin_data->critere = sopar->espilondiag;
    if (sopalin_data->critere<0.0)
    {
        /* sopalin_data->critere absolu */
        sopalin_data->critere=-sopalin_data->critere;
        if (sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_YES) {
            fprintf(stdout, "Pivoting criterium (epsilon) = %g\n", sopalin_data->critere);
        }
    } else {
        if (sopar->usenocsc == 1) {
            sopalin_data->critere = sopar->espilondiag;
        } else {
            if (sopar->fakefact == 1) {
                sopalin_data->critere = (UPDOWN_GNODENBR*UPDOWN_GNODENBR+UPDOWN_GNODENBR)*sqrt(sopar->espilondiag);
            } else {
                sopalin_data->critere = z_CscNorm1(sopalin_data->sopar->cscmtx, pastix_comm)*sqrt(sopar->espilondiag);
            }
        }
        if (sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_YES) {
            fprintf(stdout,"Pivoting criterium (||A||*sqrt(epsilon)) = %g\n", sopalin_data->critere);
        }
    }

    /* Allocating FANIN_COEFTAB */
#endif /* USE_CSC */
#endif

    for (i=0;i<SOLV_FTGTNBR;i++)
    {
        pastix_int_t ftgtsize;
        (void)ftgtsize;

#ifdef DEBUG_SOPALIN_INIT
        printf("ftgt %ld : %ld %ld %ld %ld\n",i,FANIN_LROWNUM(i),FANIN_FROWNUM(i),FANIN_LCOLNUM(i),FANIN_FCOLNUM(i));
#endif

        ftgtsize = (FANIN_LROWNUM(i)-FANIN_FROWNUM(i)+1)
            *(FANIN_LCOLNUM(i)-FANIN_FCOLNUM(i)+1);

#ifdef SOPALIN_LU
        ftgtsize *= 2;
#endif
#if (defined ALLOC_FTGT )
#if !(defined OOC_FTGT)
        FANIN_COEFTAB(i) = NULL;
#endif
#else
        MALLOC_INTERN(FANIN_COEFTAB(i), ftgtsize, pastix_complex64_t);
        for (j=0;j<ftgtsize;j++)
            FANIN_COEFTAB(i)[j] = 0.0;
#endif
    }
#ifdef DEBUG_SOPALIN_INIT
    printf("end init ftgttab\n");
#endif

    /* ???
     SYMB_NODENBR = UPDOWN_GNODENBR;
     fprintf(stdout, "Node NBR : %ld\n",SYMB_NODENBR);
     */

}

#include "z_contrib.c"

/****************************************************************************/
/*                    FACTORIZATION ROUTINE                                 */
/****************************************************************************/
/*
 * Function: z_sopalin_smp
 *
 * Factorization routine.
 *
 * This function is meant to be used when launching a thread.
 *
 * Parameters:
 *       arg - Pointer to a data structure <sopthread_data_t> with a
 *                <z_Sopalin_Data_t> pointer as *data*.
 *
 */
#define z_sopalin_smp API_CALL(z_sopalin_smp)
void* z_sopalin_smp(void *arg)
{
    sopthread_data_t *argument     = (sopthread_data_t *)arg;
    z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
    z_SolverMatrix     *datacode     = sopalin_data->datacode;
    z_SopalinParam     *sopar        = sopalin_data->sopar;
    z_Thread_Data_t    *thread_data;
    pastix_int_t               me           = argument->me;
    pastix_int_t               i            = 0;
#ifndef PASTIX_DYNSCHED
    pastix_int_t               ii           = 0;
#endif
    pastix_int_t               nbpivotT     = 0;
    int               init;
    double            maxtime;
    /*   pastix_int_t               cptinv = 0; */
#if (!(defined FORCE_NOMPI))
    MPI_Comm          pastix_comm = PASTIX_COMM;
#ifdef TEST_IRECV
    pastix_int_t               size;
#endif
#endif

#if (defined PASTIX_DYNSCHED)
    pastix_int_t bloknum;
    pastix_int_t itasktab  = me;
    pastix_int_t itasktab2 = me;
    int stolen = 0;
#endif

    MONOTHREAD_BEGIN;
    trace_begin_task(sopalin_data->tracefile,
                     SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 0,
                     STATE_L0_FACTOINIT, 0);
    MONOTHREAD_END;

    /* Initialisation des données propres à chaque thread */
    print_debug(DBG_SOPALIN_MAIN, "----- %2d : START init sopalin smp\n",
                (int)me);
    init = INIT_COMPUTE;
    if (THREAD_FUNNELED_OFF)
    {
        init = init | INIT_SEND;
    }
    if (THREAD_COMM_OFF)
    {
        init = init | INIT_RECV;
    }

    z_sopalin_init_smp(sopalin_data, me, 1, init);
    thread_data = sopalin_data->thread_data[me];
    print_debug(DBG_SOPALIN_MAIN, "----- %2d : END init sopalin smp\n", (int)me);

#ifdef TEST_IRECV
    size = PACKMAX*(sizeof(pastix_int_t)*MAXINFO)+PACKAREA*sizeof(pastix_complex64_t);
    for (i=0;i<MAX_R_REQUESTS;i++)
    {
        CALL_MPI MPI_Irecv(thread_data->recv_fanin_buffer[i],size,MPI_BYTE,
                           MPI_ANY_SOURCE,me,pastix_comm,
                           &(thread_data->recv_fanin_request[i]));
        TEST_MPI("MPI_Irecv");
    }
#endif /* TEST_IRECV */

    /* Synchro de fin d'initialisation */
    SYNCHRO_THREAD;
    MONOTHREAD_BEGIN;

#ifdef PASTIX_DUMP_FACTO
    API_CALL(z_dump_all)(datacode, sopar->cscmtx,
                       ((datacode->updovct.sm2xtab!=NULL)?
                        (DUMP_CSC | DUMP_SOLV | DUMP_SMB):
                        (DUMP_CSC | DUMP_SOLV)));
#endif

    if (THREAD_FUNNELED_OFF)
    {
        CALL_MPI MPI_Barrier(pastix_comm);
        TEST_MPI("MPI_Barrier");
    }

    if (THREAD_COMM_ON)
    {
        MUTEX_LOCK(&(sopalin_data->mutex_comm));
        sopalin_data->step_comm = COMMSTEP_FACTO;
        print_debug(DBG_THCOMM, "%s:%d FACTO\n", __FILE__, __LINE__);
        MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
        pthread_cond_broadcast(&(sopalin_data->cond_comm));
    }

    z_ooc_set_step(sopalin_data, OOCSTEP_SOPALIN);

    trace_begin_task(sopalin_data->tracefile,
                     SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 0,
                     STATE_L0_FACTO, 0);

    MONOTHREAD_END;
    SYNCHRO_THREAD;
    SOPALIN_CLOCK_INIT; /* Debut du compteur pour le temps de facto */
    COMM_CLOCK_INIT;
    if (sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        print_one("%s", OUT2_SOP_BSOP);
    print_debug(DBG_SOPALIN_MAIN,"----- [%d]%2d: sopalin starting...\n",
                (int) SOLV_PROCNUM, (int)me);

#ifdef COMPUTE_ALLOC
    ALLOC_CLOCK_INIT;
#endif

    /*****************************************************/
    /*            Main Loop                              */
    /*****************************************************/

#if defined PASTIX_DYNSCHED
    while(1){
        trace_begin_task(thread_data->tracefile,
                         SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                         STATE_WAITTSK, i);

        i = API_CALL(z_sopalin_dynsched_getNexTask)( sopalin_data, datacode, thread_data,
                                                   &itasktab, &itasktab2, &bloknum, me );

        stolen = itasktab != itasktab2;
        if ( i == -1 )
            break;

        /* GEMM tasks from ESP option */
        else if (i < -1)
        {
            i = TASK_ESP2TASK( i );
/* #ifdef ESP_WRITE */
/*             trace_begin_task1(thread_data->tracefile, */
/*                               SOPALIN_CLOCK_TRACE, */
/*                               SOLV_PROCNUM, me, */
/*                               STATE_COMP1DGEMM, */
/*                               SOLV_PROCDIAG( SYMB_CBLKNUM( bloknum ) ), */
/*                               SOLV_TASKTAB[i], */
/*                               stolen ); */
/* #else */
/*             trace_begin_task1(thread_data->tracefile, */
/*                               SOPALIN_CLOCK_TRACE, */
/*                               SOLV_PROCNUM, me, */
/*                               STATE_COMP1DGEMM, */
/*                               TASK_PROC( i ), */
/*                               SOLV_TASKTAB[i], */
/*                               stolen ); */
/* #endif */

            API_CALL(z_compute_1dgemm)(sopalin_data, me, i, bloknum, -1);
            trace_end_task1();

            MUTEX_LOCK(&(sopalin_data->tasktab_mutex[itasktab2]));
            sopalin_data->tasktab_indice[itasktab2]++;
            MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[itasktab2]));
            continue;
        }

#elif (defined SMP_SOPALIN)
        for (ii=0;ii<SOLV_TTSKNBR;ii++){
            i = SOLV_TTSKTAB(ii);

#else /* DYNSCHED, SMP_SOPALIN */
            for (ii=0;ii<SOLV_TASKNBR;ii++){
                i = queueGet(&(sopalin_data->taskqueue));
#endif /* DYNSCHED, SMP_SOPALIN */

#ifdef COMPUTE_ALLOC
                ALLOC_CLOCK_STOP;
                printf("Step %lf memsize %lf\n",ALLOC_CLOCK_GET,
                       100.0*((double)(sopalin_data->current_alloc))/((double)(SOLV_COEFNBR)));
#endif /* COMPUTE_ALLOC */

                print_debug(DBG_SOPALIN_MAIN,
                            "[%ld]%ld: z_Task %ld\n"
                            "[%ld]%ld: taskid prionum cblknum bloknum ctrcnt btagptr"
                            " indnum tasknext\n"
                            "[%ld]%ld: %ld %ld %ld %ld %ld %ld (%ld %ld %ld %ld)\n",
                            (long)SOLV_PROCNUM,(long)me,
                            (long)i,(long)SOLV_PROCNUM,(long)me,
                            (long)SOLV_PROCNUM,(long)me,
                            (long)TASK_TASKID(i),(long)TASK_PRIONUM(i),
                            (long)TASK_CBLKNUM(i),(long)TASK_BLOKNUM(i),
                            (long)TASK_CTRBCNT(i),
                            (long)TASK_INDNUM(i),
                            (long)SYMB_FCOLNUM(TASK_CBLKNUM(i)),
                            (long)SYMB_LCOLNUM(TASK_CBLKNUM(i)),
                            (long)SYMB_FROWNUM(TASK_BLOKNUM(i)),
                            (long)SYMB_LROWNUM(TASK_BLOKNUM(i)));

                trace_begin_task(thread_data->tracefile,
                                 SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                                 STATE_WAITREM, i);

                /* Compute task */
                switch(TASK_TASKID(i))
                {
                case COMP_1D:
                    /* Wait for contributions */
                    API_CALL(z_wait_contrib_comp_1d)(sopalin_data, me, i);

                    z_ooc_wait_task(sopalin_data,i, me);
                    z_ooc_wait_for_cblk(sopalin_data, TASK_CBLKNUM(i),me);


                    /* trace_begin_task1(thread_data->tracefile, */
                    /*                   SOPALIN_CLOCK_TRACE, */
                    /*                   SOLV_PROCNUM, me, */
                    /*                   STATE_COMP1D, */
                    /*                   TASK_PROC( i ), */
                    /*                   SOLV_TASKTAB[i], */
                    /*                   stolen ); */

                    /* Compute */
                    API_CALL(z_compute_1d)(sopalin_data, me, i);

                    trace_end_task1();

                    z_ooc_save_coef(sopalin_data, i, TASK_CBLKNUM(i), me);

                    break;
                default:
                    errorPrint("Taskid unknown for task %ld\n", (long)i);
                    EXIT(MOD_SOPALIN,INTERNAL_ERR);
                }

#ifdef FORCE_CONSO
                API_CALL(z_rcsd_testall_fab)(sopalin_data, me);
#endif

#ifdef PASTIX_DYNSCHED
                MUTEX_LOCK(&(sopalin_data->tasktab_mutex[itasktab2]));
                sopalin_data->tasktab_indice[itasktab2]++;
                MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[itasktab2]));
#endif
                trace_begin_task(thread_data->tracefile,
                                 SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                                 STATE_IDLE, i);
            } /* FIN boucle Principale */

#ifdef _UNUSED_
        }}
#endif
    /* Sauvegarde du temps de facto */
    SOPALIN_CLOCK_STOP;

    trace_begin_task(thread_data->tracefile,
                     SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                     STATE_IDLE, i);

    print_debug(DBG_SOPALIN_MAIN,"----- %2d-%2d : sopalin time %lf\n",
                (int)SOLV_PROCNUM, (int)me, SOPALIN_CLOCK_GET);

    /*   fprintf(stdout, "%d - %d : Nombre d'inversion %d\n", (int)SOLV_PROCNUM, (int)me, (int)cptinv); */

    /* Suppression des Comms lancées inutilement */
#if (defined TEST_IRECV) && !(defined FORCE_NOMPI)
    for (i=0; i<MAX_R_REQUESTS; i++)
    {
        int flag;
        MPI_Status status;

        CALL_MPI MPI_Cancel(&thread_data->recv_fanin_request[i]);
        TEST_MPI("MPI_Cancel");
        CALL_MPI MPI_Test(&thread_data->recv_fanin_request[i], &flag, &status);
        TEST_MPI("MPI_Test");
    }
#endif

#if (defined TEST_ISEND) && !(defined FORCE_NOMPI)
    if (THREAD_FUNNELED_OFF)
    {
        /* Attente de la fin des communications en envoi */
        if (SOLV_PROCNBR > 1)
            API_CALL(z_send_waitall_fab)(sopalin_data, me);
    }
#endif /* Fin attente comm */

#ifdef TRYLOCK
    print_debug(DBG_SOPALIN_MAIN,"----- %2d : TRYLOCK free = %4ld / busy = %4ld / wait = %4ld\n",
                (int)me, (int)thread_data->ptfree, (int)thread_data->ptbusy, (int)thread_data->ptwait);
#endif

    trace_finish(thread_data->tracefile, SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me);
    /* Nettoyage des variables associées aux threads */
    z_sopalin_clean_smp(sopalin_data, me);

    /* Synchro de tous les threads avant calcul globaux par le thread principal */
    SYNCHRO_THREAD;

    print_debug(DBG_SOPALIN_MAIN,"%d-%d : Synchro Avant reduction resultats\n",
                (int)SOLV_PROCNUM, (int)me);

    /*******************************************/
    /*           Partie Monothread             */
    /*******************************************/
    MONOTHREAD_BEGIN;

    trace_begin_task(sopalin_data->tracefile,
                     SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 0,
                     STATE_L0_FACTOCLEAN, 0);

    /* Envoi du message de fin aux threads de comm */
#ifndef FORCE_NOMPI
    if (THREAD_COMM_ON && THREAD_FUNNELED_OFF)
    {
        for (i=0; i < SOLV_PROCNBR; i++)
        {
            int tag, iterator;
            int miun = -1;
            if (i == SOLV_PROCNUM) continue;
            for (iterator=0; iterator < sopar->nbthrdcomm; iterator++)
            {
                tag = (sopar->type_comm == 3) ? iterator : TAG_FANIN;
                CALL_MPI MPI_Send(&miun, 1, MPI_INT, i, tag, PASTIX_COMM);
                TEST_MPI("MPI_Send");
                /*fprintf(stderr," %d : Envoi %d à %d\n", SOLV_PROCNUM, iterator, i);*/
            }
        }
    }
#endif

    /* Calcul du nombre de pivotage réalisé */
    for(i= 0; i< SOLV_THRDNBR; i++)
        nbpivotT += sopalin_data->thread_data[i]->nbpivot;
    sopar->diagchange = nbpivotT;

    /* Calcul du temps de facto */
    maxtime = thread_data->sop_clk;
    for(i=1; i<SOLV_THRDNBR; i++)
    {
        maxtime = MAX(sopalin_data->thread_data[i]->sop_clk, maxtime);
    }
    sopar->dparm[DPARM_FACT_TIME] = maxtime;

    /* WARNING : Don't put one (All)Reduce before thread synchronization FACTOEND */
    if (THREAD_COMM_ON)
    {
        MUTEX_LOCK(&(sopalin_data->mutex_comm));
        while(sopalin_data->step_comm != COMMSTEP_FACTOEND)
            COND_WAIT(&(sopalin_data->cond_comm), &(sopalin_data->mutex_comm));
        MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
    }

    /* Calcul de l'inertie de la matrice (pour CROUT seulement, sans pivotage) */
    sopar->iparm[IPARM_INERTIA] = -1;
#if (!defined TYPE_COMPLEX) && (!defined CHOL_SOPALIN) && (!defined OOC)
    {
        pastix_complex64_t *ga;
        pastix_int_t c, k, stride, size, inertia;
        inertia=0;
        for (c=0;c<SYMB_CBLKNBR;c++)
        {
            ga     =&(SOLV_COEFTAB(c)[SOLV_COEFIND(SYMB_BLOKNUM(c))]);
            stride = SOLV_STRIDE(c);
            size   = SYMB_LCOLNUM(c)-SYMB_FCOLNUM(c)+1;
            for (k=0;k<size;k++)
                if (ga[k+k*stride]>fzero) inertia++;
        }
        MyMPI_Allreduce(&inertia,&(sopar->iparm[IPARM_INERTIA]),1,
                        PASTIX_MPI_INT,MPI_SUM,pastix_comm);
    }
#endif

#ifdef PASTIX_DYNSCHED
    if (sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    {
        for(i=1; i<SOLV_THRDNBR; i++)
        {
            thread_data->esp += sopalin_data->thread_data[i]->esp;
        }
        MyMPI_Reduce(&(thread_data->esp), &(sopar->iparm[IPARM_ESP_NBTASKS]), 1,
                     PASTIX_MPI_INT, MPI_SUM, 0, pastix_comm);
    }
#endif

    if (THREAD_FUNNELED_OFF)
    {
        CALL_MPI MPI_Barrier(pastix_comm);
        TEST_MPI("MPI_Barrier");
    }

#ifdef STATS_SOPALIN
    if (sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
    {
        double overhead     = (double)(sopalin_data->max_alloc + SOLV_COEFNBR)/(double)(SOLV_COEFNBR);

        fprintf(stdout, "%2d - Local number of terms allocated\t Cblk+Ftgt : %10ld,\t Cblk : %10ld,\t Overhead : %.2lf (%.2lf%%)\n",
                (int)SOLV_PROCNUM,
                (long)(sopalin_data->max_alloc + SOLV_COEFNBR),
                (long)(SOLV_COEFNBR),
                overhead,
                (overhead - 1.0) * 100.0);
        {
            pastix_int_t    tmp_max_alloc = sopalin_data->max_alloc;
            pastix_int_t    tmp_coefnbr   = SOLV_COEFNBR;
            pastix_int_t    max_max_alloc = 0;
            pastix_int_t    max_coefnbr   = 0;
            pastix_int_t    sum_max_alloc = 0;
            pastix_int_t    sum_coefnbr   = 0;
            double overhead2;

            MyMPI_Allreduce(&tmp_max_alloc,&max_max_alloc,1,PASTIX_MPI_INT,MPI_MAX,pastix_comm);
            MyMPI_Allreduce(&tmp_coefnbr,  &max_coefnbr,  1,PASTIX_MPI_INT,MPI_MAX,pastix_comm);
            MyMPI_Allreduce(&tmp_max_alloc,&sum_max_alloc,1,PASTIX_MPI_INT,MPI_SUM,pastix_comm);
            MyMPI_Allreduce(&tmp_coefnbr,  &sum_coefnbr,  1,PASTIX_MPI_INT,MPI_SUM,pastix_comm);

            overhead2 = (double)(max_max_alloc+max_coefnbr)/(double)(max_coefnbr);
            print_one("Maximum number of terms allocated\t Cblk+Ftgt : %10ld,\t Cblk : %10ld,\t Overhead : %.2lf (%.2lf%%)\n",
                      (long)(max_max_alloc+max_coefnbr),
                      (long) max_coefnbr,
                      overhead2,
                      (overhead2 - 1.0) * 100.0);

            overhead2 = (double)(sum_max_alloc+sum_coefnbr)/(double)(sum_coefnbr);
            print_one("Total number of terms allocated\t\t Cblk+Ftgt : %10ld,\t Cblk : %10ld,\t Overhead : %.2lf (%.2lf%%)\n",
                      (long)(sum_max_alloc+sum_coefnbr),
                      (long) sum_coefnbr,
                      overhead2,
                      (overhead2 - 1.0) * 100.0);
            sopar->iparm[IPARM_ALLOCATED_TERMS]=sum_max_alloc+sum_coefnbr;
        }
    }
#endif

    /* free file structures */
    z_sopalin_clean(sopalin_data, 1);

#ifdef PASTIX_DUMP_FACTO
    API_CALL(z_dump_all)(datacode, sopar->cscmtx, DUMP_SOLV);
#endif

    /* Fin des threads de comms et d'OOC */
    if (THREAD_COMM_ON)
    {
        if ((sopar->iparm[IPARM_END_TASK] == API_TASK_NUMFACT)
            || (sopar->iparm[IPARM_DISTRIBUTION_LEVEL] != 0))
        {
            MUTEX_LOCK(&(sopalin_data->mutex_comm));
            sopalin_data->step_comm = COMMSTEP_END;
            print_debug(DBG_THCOMM, "%s:%d END\n", __FILE__, __LINE__);
            MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
            pthread_cond_broadcast(&(sopalin_data->cond_comm));
        }
    }

#ifdef OOC
    if (sopalin_data->sopar->iparm[IPARM_END_TASK] < API_TASK_SOLVE)
        z_ooc_stop_thread(sopalin_data);
#endif

    trace_begin_task(sopalin_data->tracefile,
                     SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 0,
                     STATE_L0_IDLE, 0);

    MONOTHREAD_END;
    if (sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        print_one("%s", OUT2_SOP_ESOP);

    if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_CHATTERBOX)
        fprintf(stdout, OUT4_FACT_COMM_TIME,
                (int)SOLV_PROCNUM, (int)me, COMM_CLOCK_GET);
    return 0;
}

/******************************************************************************/
/*         Fonction pour les threads de communication                         */
/******************************************************************************/
/*
 Function: API_CALL(z_sopalin_updo_comm)

 Updown function used for the communication thread.

 Parameters:
 arg - Pointer to a data structure <sopthread_data_t> with a
 <z_Sopalin_Data_t> pointer as *data*.
 */
#define z_sopalin_updo_comm API_CALL(z_sopalin_updo_comm)
void *z_sopalin_updo_comm ( void *arg )
{
#ifndef FORCE_NOMPI
    sopthread_data_t *argument     = (sopthread_data_t *)arg;
    z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
    z_SolverMatrix     *datacode     = sopalin_data->datacode;
    pastix_int_t               me           = argument->me;
    if (THREAD_COMM_ON)
    {
        /* z_Thread_Data_t    *thread_data  = sopalin_data->thread_data[me]; */

        /* On se met en attente du debut de la descente
         * ou on quitte si on ne reitere pas */
        MUTEX_LOCK(&(sopalin_data->mutex_comm));
        while(1)
        {
            switch(sopalin_data->step_comm)
            {
                /* Il n'y a plus de comm */
            case COMMSTEP_END:
                sopalin_data->step_comm = COMMSTEP_INIT;
                print_debug(DBG_THCOMM, "%s:%d INIT\n", __FILE__, __LINE__);
                MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
                print_debug(DBG_SOPALIN_THREADCOMM,
                            "%d : je quitte\n", (int)SOLV_PROCNUM);
                return (void *)1;
                break;

                /* Factorisation */
            case COMMSTEP_FACTO:
                MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
                API_CALL(sendrecv_smp)(arg);
                /* On ne lance qu'une facto et
                 * le reste ne necessite qu'un thread pour l'instant */
                if (me > SOLV_THRDNBR) return (void *)1;
                MUTEX_LOCK(&(sopalin_data->mutex_comm));
                break;

                /* Udpo */
            case COMMSTEP_DOWN:
                MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
                API_CALL(z_updo_thread_comm)(arg);
                MUTEX_LOCK(&(sopalin_data->mutex_comm));
                break;

                /* Un AllReduce est a faire en funneled */
            case COMMSTEP_ALLREDUCE:
                if (THREAD_FUNNELED_ON)
                {
                    {
                        Pastix_Allreduce_t *allreduce = &(sopalin_data->allreduce);
                        MPI_Allreduce(allreduce->sendbuf,
                                      allreduce->recvbuf,
                                      allreduce->count,
                                      allreduce->datatype,
                                      allreduce->op,
                                      PASTIX_COMM);
                        sopalin_data->step_comm = COMMSTEP_INIT;
                        print_debug(DBG_THCOMM, "%s:%d INIT\n", __FILE__, __LINE__);
                        pthread_cond_broadcast(&(sopalin_data->cond_comm));
                    }
                    break;
                }
            case COMMSTEP_REDUCE:
                if (THREAD_FUNNELED_ON)
                {
                    {
                        Pastix_Allreduce_t *reduce = &(sopalin_data->allreduce);
                        MPI_Reduce(reduce->sendbuf,
                                   reduce->recvbuf,
                                   reduce->count,
                                   reduce->datatype,
                                   reduce->op,
                                   0,
                                   PASTIX_COMM);
                        sopalin_data->step_comm = COMMSTEP_INIT;
                        print_debug(DBG_THCOMM, "%s:%d INIT\n", __FILE__, __LINE__);
                        pthread_cond_broadcast(&(sopalin_data->cond_comm));
                    }
                    break;
                }
                /* Sortie spontannée du wait */
            default:
                COND_WAIT(&(sopalin_data->cond_comm),
                          &(sopalin_data->mutex_comm));
            }
        }
        MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
    }
#endif
    (void)arg;
    return (void*)NULL;
}

/******************************************************************************/
/*         Fonction pour les threads de calculs                               */
/******************************************************************************/

/*
 * Function: z_sopalin_thread
 *
 * Function launching computing, communictating and out of core threads on
 * the factorization step only.
 *
 * Initiate the <z_Sopalin_Data_t> structure, launch threads, clean and restore.
 *
 * Parameters:
 *       m         - The <z_SolverMatrix> structure.
 *	sopaparam - Sopalin parameters in the <z_SopalinParam> stucture.
 */
#define z_sopalin_thread API_CALL(z_sopalin_thread)
void z_sopalin_thread(z_SolverMatrix *m,
                    z_SopalinParam *sopaparam)
{
    z_Backup b;
    z_Sopalin_Data_t *sopalin_data = NULL;
    z_SolverMatrix  *datacode = NULL;

    MALLOC_INTERN(sopalin_data, 1, z_Sopalin_Data_t);

    z_sopalin_backup(m,&b);
    z_sopalin_init(sopalin_data, m, sopaparam, 1);
    API_CALL(z_init_struct_sopalin)(sopalin_data, m, sopaparam);

    datacode = sopalin_data->datacode;

#ifdef PASTIX_WITH_STARPU
    if (sopalin_data->sopar->iparm[IPARM_STARPU] == API_YES)
    {
        starpu_zsubmit_tasks(sopalin_data);
    }
    else
#endif
    {
        sopalin_launch_thread(sopalin_data,
                              SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree,
                              sopalin_data->sopar->iparm[IPARM_VERBOSE],
                              SOLV_THRDNBR,          API_CALL(z_sopalin_smp),       sopalin_data,
                              sopaparam->nbthrdcomm, API_CALL(z_sopalin_updo_comm), sopalin_data,
                              OOC_THREAD_NBR,        z_ooc_thread,                  sopalin_data);
    }
    z_sopalin_clean(sopalin_data, 2);
    z_sopalin_restore(m,&b);

    memFree_null(sopalin_data);
}

/*
 Function: z_sopalin_updo_smp

 Function used for computing thread creation to compute factorisation and
 resolution.

 Parameters:
 arg - Pointer to a data structure <sopthread_data_t> with a
 <z_Sopalin_Data_t> pointer as *data*.
 */
void* API_CALL(z_sopalin_updo_smp)(void *arg)
{
    sopthread_data_t *argument     = (sopthread_data_t *)arg;
    z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
    pastix_int_t               me           = argument->me;

    API_CALL(z_sopalin_smp)(argument);
    if (sopalin_data->sopar->iparm[IPARM_DISTRIBUTION_LEVEL] != 0)
    {
        if ((sopalin_data->datacode->clustnum == 0) && (me == 0))
            errorPrintW("Updown incompatible with 2D distribution");
        return 0;
    }

    MONOTHREAD_BEGIN;
    z_sopalin_init(sopalin_data, NULL, NULL, 0);
    MONOTHREAD_END;
    API_CALL(z_up_down_smp)(argument);

    return 0;
}
/*
 Function: API_CALL(z_sopalin_updo_thread)

 Function launching computing, communicating and out of core threads on
 the factorization and solve steps.

 Initiate the <z_Sopalin_Data_t> structure, launch threads, clean and restore.

 Parameters:
 m         - The <z_SolverMatrix> structure.
 sopaparam - Sopalin parameters in the <z_SopalinParam> stucture.
 */
void API_CALL(z_sopalin_updo_thread)(z_SolverMatrix *m,
                                   z_SopalinParam *sopaparam)
{
    z_Backup b;
    z_Sopalin_Data_t *sopalin_data = NULL;
    z_SolverMatrix   *datacode = NULL;

    MALLOC_INTERN(sopalin_data, 1, z_Sopalin_Data_t);

    z_sopalin_backup(m,&b);
    z_sopalin_init(sopalin_data, m, sopaparam, 1);
    API_CALL(z_init_struct_sopalin)(sopalin_data, m, sopaparam);
    datacode = sopalin_data->datacode;
#ifdef PASTIX_WITH_STARPU
    if (sopalin_data->sopar->iparm[IPARM_STARPU] == API_YES)
    {
        starpu_zsubmit_tasks(sopalin_data);

    }
    else
#endif
    {
        sopalin_launch_thread(sopalin_data,
                              SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree, sopalin_data->sopar->iparm[IPARM_VERBOSE],
                              SOLV_THRDNBR,          API_CALL(z_sopalin_updo_smp),  sopalin_data,
                              sopaparam->nbthrdcomm, API_CALL(z_sopalin_updo_comm), sopalin_data,
                              OOC_THREAD_NBR,        z_ooc_thread,                  sopalin_data);
    }

    z_sopalin_clean(sopalin_data, 2);
    z_sopalin_restore(m,&b);

    memFree_null(sopalin_data);
}

/*
 Function: API_CALL(z_sopalin_updo_gmres_smp)

 Function used for computing thread creation to compute factorisation,
 resolution and gmres.

 Parameters:
 arg - Pointer to a data structure <sopthread_data_t> with a
 <z_Sopalin_Data_t> pointer as *data*.
 */
void* API_CALL(z_sopalin_updo_gmres_smp)(void *arg)
{
    sopthread_data_t *argument     = (sopthread_data_t *)arg;
    z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
    pastix_int_t               me           = argument->me;

    API_CALL(z_sopalin_smp)(argument);
    if (sopalin_data->sopar->iparm[IPARM_DISTRIBUTION_LEVEL] != 0)
    {
        if ((sopalin_data->datacode->clustnum == 0) && (me == 0))
            errorPrintW("Updown incompatible with 2D distribution");
        return 0;
    }

    MONOTHREAD_BEGIN;
    z_sopalin_init(sopalin_data, NULL, NULL, 0);
    MONOTHREAD_END;
    API_CALL(z_gmres_smp)(argument);

    return 0;
}
/*
 Function: API_CALL(z_sopalin_updo_gmres_thread)

 Function launching computing, communicating and out of core threads on
 the factorization, solve and reffinement (using GMRES) steps.

 Initiate the <z_Sopalin_Data_t> structure, launch threads, clean and restore.

 Parameters:
 m         - The <z_SolverMatrix> structure.
 sopaparam - Sopalin parameters in the <z_SopalinParam> stucture.
 */
void API_CALL(z_sopalin_updo_gmres_thread)(z_SolverMatrix *m, z_SopalinParam *sopaparam)
{
    z_Backup b;
    z_Sopalin_Data_t *sopalin_data;
    z_SolverMatrix   *datacode = m;



    MALLOC_INTERN(sopalin_data, 1, z_Sopalin_Data_t);

    z_sopalin_backup(m,&b);
    z_sopalin_init(sopalin_data, m, sopaparam, 1);
    API_CALL(z_init_struct_sopalin)(sopalin_data, m, sopaparam);
    datacode = sopalin_data->datacode;
#ifdef PASTIX_WITH_STARPU
    if (sopalin_data->sopar->iparm[IPARM_STARPU] == API_YES)
    {

        starpu_zsubmit_tasks(sopalin_data);

    }
    else
#endif
    {
        sopalin_launch_thread(sopalin_data,
                              SOLV_PROCNUM,          SOLV_PROCNBR,                     datacode->btree,
                              sopalin_data->sopar->iparm[IPARM_VERBOSE],
                              SOLV_THRDNBR,          API_CALL(z_sopalin_updo_gmres_smp), sopalin_data,
                              sopaparam->nbthrdcomm, API_CALL(z_sopalin_updo_comm),      sopalin_data,
                              OOC_THREAD_NBR,        z_ooc_thread,                       sopalin_data);
    }
    z_sopalin_clean(sopalin_data, 2);
    z_sopalin_restore(m,&b);

    memFree_null(sopalin_data);
}

/*
 Function: API_CALL(z_sopalin_updo_grad_smp)

 Function used for computing thread creation to compute factorisation,
 resolution and conjugate gradient.

 Parameters:
 arg - Pointer to a data structure <sopthread_data_t> with a
 <z_Sopalin_Data_t> pointer as *data*.
 */
void* API_CALL(z_sopalin_updo_grad_smp)(void *arg)
{
    sopthread_data_t *argument     = (sopthread_data_t *)arg;
    z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
    pastix_int_t               me           = argument->me;

    API_CALL(z_sopalin_smp)(argument);
    if (sopalin_data->sopar->iparm[IPARM_DISTRIBUTION_LEVEL] != 0)
    {
        if ((sopalin_data->datacode->clustnum == 0) && (me == 0))
            errorPrintW("Updown incompatible with 2D distribution");
        return 0;
    }

    MONOTHREAD_BEGIN;
    z_sopalin_init(sopalin_data, NULL, NULL, 0);
    MONOTHREAD_END;
    API_CALL(z_up_down_smp)(argument);
    API_CALL(z_grad_smp)   (argument);

    return 0;
}

/*
 Function: API_CALL(z_sopalin_updo_grad_thread)

 Function launching computing, communicating and out of core threads on
 the factorization, solve and reffinement (using conjugate grandient) steps.

 Initiate the <z_Sopalin_Data_t> structure, launch threads, clean and restore.

 Parameters:
 m         - The <z_SolverMatrix> structure.
 sopaparam - Sopalin parameters in the <z_SopalinParam> stucture.
 */
void API_CALL(z_sopalin_updo_grad_thread)(z_SolverMatrix *m, z_SopalinParam *sopaparam)
{
    z_Backup b;
    z_Sopalin_Data_t *sopalin_data = NULL;
    z_SolverMatrix   *datacode = NULL;

    MALLOC_INTERN(sopalin_data, 1, z_Sopalin_Data_t);

    z_sopalin_backup(m,&b);
    z_sopalin_init(sopalin_data, m, sopaparam, 1);
    API_CALL(z_init_struct_sopalin)(sopalin_data, m, sopaparam);
    datacode = sopalin_data->datacode;
#ifdef PASTIX_WITH_STARPU
    if (sopalin_data->sopar->iparm[IPARM_STARPU] == API_YES)
    {

        starpu_zsubmit_tasks(sopalin_data);

    }
    else
#endif
    {
        sopalin_launch_thread(sopalin_data,
                              SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree, sopalin_data->sopar->iparm[IPARM_VERBOSE],
                              SOLV_THRDNBR,          API_CALL(z_sopalin_updo_grad_smp), sopalin_data,
                              sopaparam->nbthrdcomm, API_CALL(z_sopalin_updo_comm),     sopalin_data,
                              OOC_THREAD_NBR,        z_ooc_thread,                      sopalin_data);
    }
    z_sopalin_clean(sopalin_data, 2);
    z_sopalin_restore(m,&b);

    memFree_null(sopalin_data);
}

/*
 Function: API_CALL(z_sopalin_updo_pivot_smp)

 Function used for computing thread creation to compute factorisation,
 resolution and pivoting refinement.

 Parameters:
 arg - Pointer to a data structure <sopthread_data_t> with a
 <z_Sopalin_Data_t> pointer as *data*.
 */
void* API_CALL(z_sopalin_updo_pivot_smp)(void *arg)
{
    sopthread_data_t *argument     = (sopthread_data_t *)arg;
    z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
    pastix_int_t               me           = argument->me;

    API_CALL(z_sopalin_smp)(argument);
    if (sopalin_data->sopar->iparm[IPARM_DISTRIBUTION_LEVEL] != 0)
    {
        if ((sopalin_data->datacode->clustnum == 0) && (me == 0))
            errorPrintW("Updown incompatible with 2D distribution");
        return 0;
    }

    MONOTHREAD_BEGIN;
    z_sopalin_init(sopalin_data, NULL, NULL, 0);
    MONOTHREAD_END;
    API_CALL(z_up_down_smp)(argument);
    API_CALL(z_pivotstatique_smp)(argument);

    return 0;
}

/*
 Function: API_CALL(z_sopalin_updo_pivot_thread)

 Function launching computing, communicating and out of core threads on
 the factorization, solve and reffinement (using pivoting refinement) steps.

 Initiate the <z_Sopalin_Data_t> structure, launch threads, clean and restore.

 Parameters:
 m         - The <z_SolverMatrix> structure.
 sopaparam - Sopalin parameters in the <z_SopalinParam> stucture.
 */
void API_CALL(z_sopalin_updo_pivot_thread)(z_SolverMatrix *m, z_SopalinParam *sopaparam)
{
    z_Backup b;
    z_Sopalin_Data_t *sopalin_data = NULL;
    z_SolverMatrix   *datacode = NULL;

    MALLOC_INTERN(sopalin_data, 1, z_Sopalin_Data_t);

    z_sopalin_backup(m,&b);
    z_sopalin_init(sopalin_data, m, sopaparam, 1);
    datacode = sopalin_data->datacode;
#ifdef PASTIX_WITH_STARPU
    if (sopalin_data->sopar->iparm[IPARM_STARPU] == API_YES)
    {

        starpu_zsubmit_tasks(sopalin_data);

    }
    else
#endif
    {
        sopalin_launch_thread(sopalin_data,
                              SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree, sopalin_data->sopar->iparm[IPARM_VERBOSE],
                              SOLV_THRDNBR,          API_CALL(z_sopalin_updo_pivot_smp), sopalin_data,
                              sopaparam->nbthrdcomm, API_CALL(z_sopalin_updo_comm),      sopalin_data,
                              OOC_THREAD_NBR,        z_ooc_thread,                       sopalin_data);
    }
    z_sopalin_clean(sopalin_data, 2);
    z_sopalin_restore(m,&b);
}


/*
 Function: API_CALL(sopalin_updo_bicgstab_smp)

 Function used for computing thread creation to compute factorisation,
 resolution and bicgstab refinement.

 Parameters:
 arg - Pointer to a data structure <sopthread_data_t> with a
 <z_Sopalin_Data_t> pointer as *data*.
 */
void* API_CALL(sopalin_updo_bicgstab_smp)(void *arg)
{
    sopthread_data_t *argument     = (sopthread_data_t *)arg;
    z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
    pastix_int_t               me           = argument->me;

    API_CALL(z_sopalin_smp)(argument);
    if (sopalin_data->sopar->iparm[IPARM_DISTRIBUTION_LEVEL] != 0)
    {
        if ((sopalin_data->datacode->clustnum == 0) && (me == 0))
            errorPrintW("Updown incompatible with 2D distribution");
        return 0;
    }

    MONOTHREAD_BEGIN;
    z_sopalin_init(sopalin_data, NULL, NULL, 0);
    MONOTHREAD_END;
    API_CALL(z_up_down_smp)(argument);
    API_CALL(z_bicgstab_smp)(argument);

    return 0;
}
/*
 Function: API_CALL(sopalin_updo_bicgstab_thread)

 Function launching computing, communicating and out of core threads on
 the factorization, solve and reffinement (using bicgstab refinement) steps.

 Initiate the <z_Sopalin_Data_t> structure, launch threads, clean and restore.

 Parameters:
 m         - The <z_SolverMatrix> structure.
 sopaparam - Sopalin parameters in the <z_SopalinParam> stucture.
 */
void API_CALL(sopalin_updo_bicgstab_thread)(z_SolverMatrix *m, z_SopalinParam *sopaparam)
{
    z_Backup b;
    z_Sopalin_Data_t *sopalin_data = NULL;
    z_SolverMatrix   *datacode = NULL;

    MALLOC_INTERN(sopalin_data, 1, z_Sopalin_Data_t);

    z_sopalin_backup(m,&b);
    z_sopalin_init(sopalin_data, m, sopaparam, 1);
    datacode = sopalin_data->datacode;
#ifdef PASTIX_WITH_STARPU
    if (sopalin_data->sopar->iparm[IPARM_STARPU] == API_YES)
    {

        starpu_zsubmit_tasks(sopalin_data);

    }
    else
#endif
    {
        sopalin_launch_thread(sopalin_data,
                              SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree, sopalin_data->sopar->iparm[IPARM_VERBOSE],
                              SOLV_THRDNBR,          API_CALL(sopalin_updo_bicgstab_smp), sopalin_data,
                              sopaparam->nbthrdcomm, API_CALL(z_sopalin_updo_comm),      sopalin_data,
                              OOC_THREAD_NBR,        z_ooc_thread,                       sopalin_data);
    }
    z_sopalin_clean(sopalin_data, 2);
    z_sopalin_restore(m,&b);
}

/*
 Function: API_CALL(z_sopalin_launch)

 TODO: Comment (unused ?)
 */
void API_CALL(z_sopalin_launch)(z_SolverMatrix *m,
                              z_SopalinParam *sopaparam,
                              pastix_int_t cas)
{
    z_Backup b;
    z_Sopalin_Data_t *sopalin_data = NULL;
    z_SolverMatrix   *datacode     = NULL;

    MALLOC_INTERN(sopalin_data, 1, z_Sopalin_Data_t);

    if (cas < UPDO_ONLY)
    {
        z_sopalin_backup(m,&b);
        z_sopalin_init(sopalin_data, m, sopaparam, 1);
        API_CALL(z_init_struct_sopalin)(sopalin_data, m, sopaparam);
    }
    else
    {
        z_sopalin_init(sopalin_data, m, sopaparam, 0);
    }

    datacode = sopalin_data->datacode;
    switch(cas){
    case SOPALIN_ONLY:
        sopalin_launch_thread(sopalin_data,
                              SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree, sopalin_data->sopar->iparm[IPARM_VERBOSE],
                              SOLV_THRDNBR,          API_CALL(z_sopalin_smp),       sopalin_data,
                              sopaparam->nbthrdcomm, API_CALL(z_sopalin_updo_comm), sopalin_data,
                              OOC_THREAD_NBR,        z_ooc_thread,                  sopalin_data);
        break;
    case SOPALIN_UPDO:
        sopalin_launch_thread(sopalin_data,
                              SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree, sopalin_data->sopar->iparm[IPARM_VERBOSE],
                              SOLV_THRDNBR,          API_CALL(z_sopalin_updo_smp),  sopalin_data,
                              sopaparam->nbthrdcomm, API_CALL(z_sopalin_updo_comm), sopalin_data,
                              OOC_THREAD_NBR,        z_ooc_thread,                  sopalin_data);
        break;
    case SOPALIN_UPDO_GMRES:
        sopalin_launch_thread(sopalin_data,
                              SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree, sopalin_data->sopar->iparm[IPARM_VERBOSE],
                              SOLV_THRDNBR,          API_CALL(z_sopalin_updo_gmres_smp), sopalin_data,
                              sopaparam->nbthrdcomm, API_CALL(z_sopalin_updo_comm),      sopalin_data,
                              OOC_THREAD_NBR,        z_ooc_thread,                       sopalin_data);
        break;
    case SOPALIN_UPDO_GRAD:
        sopalin_launch_thread(sopalin_data,
                              SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree, sopalin_data->sopar->iparm[IPARM_VERBOSE],
                              SOLV_THRDNBR,          API_CALL(z_sopalin_updo_grad_smp), sopalin_data,
                              sopaparam->nbthrdcomm, API_CALL(z_sopalin_updo_comm),     sopalin_data,
                              OOC_THREAD_NBR,        z_ooc_thread,                      sopalin_data);
        break;
    case SOPALIN_UPDO_PIVOT:
        sopalin_launch_thread(sopalin_data,
                              SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree, sopalin_data->sopar->iparm[IPARM_VERBOSE],
                              SOLV_THRDNBR,          API_CALL(z_sopalin_updo_pivot_smp), sopalin_data,
                              sopaparam->nbthrdcomm, API_CALL(z_sopalin_updo_comm),      sopalin_data,
                              OOC_THREAD_NBR,        z_ooc_thread,                       sopalin_data);
        break;
    case UPDO_ONLY:
        sopalin_launch_thread(sopalin_data,
                              SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree, sopalin_data->sopar->iparm[IPARM_VERBOSE],
                              SOLV_THRDNBR,          API_CALL(z_up_down_smp),       sopalin_data,
                              sopaparam->nbthrdcomm, API_CALL(z_sopalin_updo_comm), sopalin_data,
                              OOC_THREAD_NBR,        z_ooc_thread,                  sopalin_data);
        break;
    case RAFF_GMRES:
        sopalin_launch_thread(sopalin_data,
                              SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree, sopalin_data->sopar->iparm[IPARM_VERBOSE],
                              SOLV_THRDNBR,          API_CALL(z_gmres_smp),         sopalin_data,
                              sopaparam->nbthrdcomm, API_CALL(z_sopalin_updo_comm), sopalin_data,
                              OOC_THREAD_NBR,        z_ooc_thread,                  sopalin_data);
        break;
    case RAFF_GRAD:
        sopalin_launch_thread(sopalin_data,
                              SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree, sopalin_data->sopar->iparm[IPARM_VERBOSE],
                              SOLV_THRDNBR,          API_CALL(z_grad_smp),          sopalin_data,
                              sopaparam->nbthrdcomm, API_CALL(z_sopalin_updo_comm), sopalin_data,
                              OOC_THREAD_NBR,        z_ooc_thread,                  sopalin_data);
        break;
    case RAFF_PIVOT:
        sopalin_launch_thread(sopalin_data,
                              SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree, sopalin_data->sopar->iparm[IPARM_VERBOSE],
                              SOLV_THRDNBR,          API_CALL(z_pivotstatique_smp), sopalin_data,
                              sopaparam->nbthrdcomm, API_CALL(z_sopalin_updo_comm), sopalin_data,
                              OOC_THREAD_NBR,        z_ooc_thread,                  sopalin_data);
        break;
    default:
        if( SOLV_PROCNUM == 0 )
        {
            errorPrint("undefined case.");
            EXIT(MOD_SOPALIN,BADPARAMETER_ERR);
        }
    }

    z_sopalin_clean(sopalin_data, 2);
    if (cas < UPDO_ONLY)
        z_sopalin_restore(m,&b);

    memFree_null(sopalin_data);
}
