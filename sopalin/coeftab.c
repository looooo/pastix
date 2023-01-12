/**
 *
 * @file coeftab.c
 *
 * PaStiX coefficient array initialization and free routines.
 *
 * @copyright 2015-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Esragul Korkmaz
 * @author Gregoire Pichon
 * @author Tony Delarue
 * @date 2022-11-17
 *
 **/
#include "common.h"
#include "bcsc/bcsc.h"
#include "isched.h"
#include "blend/solver.h"
#include "sopalin/coeftab.h"
#include "pastix_zcores.h"
#include "pastix_ccores.h"
#include "pastix_dcores.h"
#include "pastix_scores.h"
#include "cpucblk_zpack.h"
#include "cpucblk_cpack.h"
#include "cpucblk_dpack.h"
#include "cpucblk_spack.h"

#include "pastix_zccores.h"
#include "pastix_dscores.h"

#if defined(PASTIX_WITH_PARSEC)
#include "sopalin/parsec/pastix_parsec.h"
#endif

#if defined(PASTIX_WITH_STARPU)
#include "sopalin/starpu/pastix_starpu.h"
#endif

coeftab_fct_memory_t coeftabMemory[4] =
{
    coeftab_smemory, coeftab_dmemory, coeftab_cmemory, coeftab_zmemory
};

static void (*initfunc[2][4])( pastix_coefside_t, const SolverMatrix*,
                               const pastix_bcsc_t*, pastix_int_t, const char *) =
{
    /* Normal precision functions */
    {
        cpucblk_sinit,
        cpucblk_dinit,
        cpucblk_cinit,
        cpucblk_zinit
    },
    /* Mixed-precision functions */
    {
        cpucblk_sinit,
        cpucblk_dsinit,
        cpucblk_cinit,
        cpucblk_zcinit
    }
};

/**
 * @brief Internal structure specific to the parallel call of pcoeftabInit()
 */
struct coeftabinit_s {
    const SolverMatrix  *datacode; /**< The solver matrix                         */
    const pastix_bcsc_t *bcsc;     /**< The internal block CSC                    */
    const char          *dirname;  /**< The pointer to the output directory       */
    pastix_coefside_t    side;     /**< The side of the matrix beeing initialized */
    pastix_int_t         mixed;    /**< The mixed-precision parameter             */
};

/**
 *******************************************************************************
 *
 * @brief Allocates the entire coeftab matrix with a single allocation.
 *
 * The different cblk coeftabs are assigned inside the
 * allocation, each one with their right size.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *          On exit, the internal coeftab matrix is allocated.
 *
 *******************************************************************************/
void
coeftabAlloc( pastix_data_t *pastix_data )
{
    SolverMatrix     *solvmatr = pastix_data->solvmatr;
    SolverCblk       *cblk     = solvmatr->cblktab;
    pastix_int_t      i, step  = 0;
    pastix_coeftype_t flttype  = solvmatr->flttype;
    size_t            size     = solvmatr->coefnbr * pastix_size_of( flttype );
    char             *workL    = NULL;
    char             *workU    = NULL;

    workL = pastix_malloc_pinned( size );
    memset( workL, 0, size );

    /* Only allocates the U part if necessary */
    if ( pastix_data->iparm[IPARM_FACTORIZATION] == PastixFactLU ) {
        workU = pastix_malloc_pinned( size );
        memset( workU, 0, size );
    }

    /*
     * Assign the cblks to their corresponding index in work
     * lcoeftabs and ucoeftabs are both one allocation
     */
    for ( i=0; i<solvmatr->cblknbr; i++, cblk++ ) {

        /* FANIN and RECV cblks are allocated independently */
        if ( cblk->cblktype & (CBLK_RECV|CBLK_FANIN) ) {
            continue;
        }

        assert( cblk->lcoeftab == NULL );
        assert( (size_t) (step) < size );
        cblk->lcoeftab = workL + step;

        if ( pastix_data->iparm[IPARM_FACTORIZATION] == PastixFactLU ) {
            assert( cblk->ucoeftab == NULL );
            cblk->ucoeftab = workU + step;
        }

        step += cblk_colnbr( cblk ) * cblk->stride * pastix_size_of( flttype );
    }
}

/**
 *******************************************************************************
 *
 * @brief Internal routine called by each static thread to Initialize the solver
 * matrix structure.
 *
 * This routine is the routine called by each thread in the static scheduler and
 * launched by the coeftabinit().
 *
 *******************************************************************************
 *
 * @param[inout] ctx
 *          The internal scheduler context
 *
 * @param[in] args
 *          The data structure specific to the function cpucblk_zinit()
 *
 *******************************************************************************/
void
pcoeftabInit( isched_thread_t *ctx,
              void            *args )
{
    struct coeftabinit_s *ciargs   = (struct coeftabinit_s*)args;
    const SolverMatrix   *datacode = ciargs->datacode;
    const pastix_bcsc_t  *bcsc     = ciargs->bcsc;
    const char           *dirname  = ciargs->dirname;
    pastix_coefside_t     side     = ciargs->side;
    pastix_int_t          mixed    = ciargs->mixed;
    pastix_int_t i, itercblk;
    pastix_int_t task;
    int rank = ctx->rank;

    for (i=0; i < datacode->ttsknbr[rank]; i++)
    {
        task = datacode->ttsktab[rank][i];
        itercblk = datacode->tasktab[task].cblknum;

        /* Init as full rank */
        initfunc[mixed][bcsc->flttype - 2]( side, datacode, bcsc, itercblk, dirname );
    }
}

/**
 *******************************************************************************
 *
 * @brief Initialize the solver matrix structure
 *
 * This routine is a parallel routine to initialize the solver matrix structure
 * through the internal static scheduler
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that hold the solver matrix to initialize.
 *
 * @param[in] side
 *          Describe the side(s) of the matrix that must be initialized.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 *******************************************************************************/
void
coeftabInit( pastix_data_t    *pastix_data,
             pastix_coefside_t side )
{
    struct coeftabinit_s args;

    pastix_data->solvmatr->globalalloc = pastix_data->iparm[IPARM_GLOBAL_ALLOCATION];
    args.datacode = pastix_data->solvmatr;
    args.bcsc     = pastix_data->bcsc;
    args.side     = side;
    args.mixed    = pastix_data->iparm[IPARM_MIXED];

    /* Allocates the coeftab matrix before multi-threading if global allocation is enabled */
    if ( args.datacode->globalalloc )
    {
        if ( pastix_data->iparm[IPARM_COMPRESS_WHEN] != PastixCompressNever ) {
            pastix_print_warning( "Global allocation is not allowed with compression. It is disabled\n" );
            pastix_data->solvmatr->globalalloc = 0;
        }
        else {
            /* Make sure there is no compression */
            assert( pastix_data->iparm[IPARM_COMPRESS_WHEN] == PastixCompressNever );
            coeftabAlloc( pastix_data );
        }
    }

#if defined(PASTIX_DEBUG_DUMP_COEFTAB)
    /* Make sure dir_local is initialized before calling it with multiple threads */
    pastix_gendirectories( pastix_data );
#endif
    args.dirname = pastix_data->dir_local;

    /* Initialize ILU level block variables */
    {
        SolverCblk *cblk = pastix_data->solvmatr->cblktab;
        SolverBlok *blok = pastix_data->solvmatr->bloktab;
        pastix_int_t i;

        for (i=0; i<pastix_data->solvmatr->cblknbr; i++, cblk++) {
            cblk->ctrbcnt = cblk[1].brownum - cblk[0].brownum;
        }

        for (i=0; i<pastix_data->solvmatr->bloknbr; i++, blok++) {
            blok->iluklvl = INT_MAX;
        }
    }
    isched_parallel_call( pastix_data->isched, pcoeftabInit, &args );
}

/**
 *******************************************************************************
 *
 * @brief Free the solver matrix structure
 *
 * This routine free all data structure refereing to the solver matrix L, even
 * the runtime descriptors if present.
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The solver matrix structure of the problem.
 *
 *******************************************************************************/
void
coeftabExit( SolverMatrix *solvmtx )
{
    pastix_int_t i;

#if defined(PASTIX_WITH_PARSEC)
    {
        if ( solvmtx->parsec_desc != NULL ) {
            parsec_sparse_matrix_destroy( solvmtx->parsec_desc );
            free( solvmtx->parsec_desc );
        }
        solvmtx->parsec_desc = NULL;
    }
#endif
#if defined(PASTIX_WITH_STARPU)
    {
        if ( solvmtx->starpu_desc != NULL ) {
            starpu_sparse_matrix_destroy( solvmtx->starpu_desc );
            free( solvmtx->starpu_desc );
        }
        solvmtx->starpu_desc = NULL;
    }
#endif

    /* If the coeftab is a single allocation, only free the first block */
    if ( solvmtx->globalalloc ) {
        pastix_free_pinned( solvmtx->cblktab->lcoeftab );
        if ( solvmtx->cblktab->ucoeftab ) {
            pastix_free_pinned( solvmtx->cblktab->ucoeftab );
        }
    }

    /* Free arrays of solvmtx */
    if ( solvmtx->cblktab )
    {
        SolverCblk *cblk = solvmtx->cblktab;

        for ( i = 0; i < solvmtx->cblknbr; i++, cblk++ ) {
            if ( cblk->cblktype & (CBLK_FANIN|CBLK_RECV) ) {
                continue;
            }
            /* If the blocks are already freed set them to NULL, else free them */
            if ( solvmtx->globalalloc ) {
                cblk->lcoeftab = NULL;
                if ( cblk->ucoeftab ) {
                    cblk->ucoeftab = NULL;
                }
            }
            /* Free is precision independent, so we can use any version */
            else {
                cpucblk_zfree( PastixLUCoef, cblk );
            }
        }
    }
}

/**
 * @brief Internal structure specific to the parallel call of pcoeftabComp()
 */
struct coeftabcomp_s {
    SolverMatrix        *solvmtx; /**< The solver matrix               */
    pastix_coeftype_t    flttype; /**< The arithmetic type             */
    pastix_atomic_lock_t lock;    /**< Lock to protect the gain update */
    pastix_int_t         gain;    /**< The memory gain on output       */
};

/**
 *******************************************************************************
 *
 * @brief Internal routine called by each static thread to Initialize the solver
 * matrix structure.
 *
 * This routine is the routine called by each thread in the static scheduler and
 * launched by the coeftabCompress().
 *
 *******************************************************************************
 *
 * @param[inout] ctx
 *          The internal scheduler context
 *
 * @param[in] args
 *          The data structure specific to the function cpucblk_zcompress()
 *
 *******************************************************************************/
static void
pcoeftabComp( isched_thread_t *ctx,
              void            *args )
{
    struct coeftabcomp_s *ccargs   = (struct coeftabcomp_s*)args;
    SolverMatrix         *solvmtx  = ccargs->solvmtx;
    pastix_coeftype_t     flttype  = ccargs->flttype;
    pastix_atomic_lock_t *lock     = &(ccargs->lock);
    pastix_int_t         *fullgain = &(ccargs->gain);
    SolverCblk           *cblk;
    pastix_coefside_t side = (solvmtx->factotype == PastixFactLU) ? PastixLUCoef : PastixLCoef;
    pastix_int_t i, itercblk;
    pastix_int_t task, gain = 0;
    int rank = ctx->rank;

    pastix_int_t (*compfunc)( const SolverMatrix*, pastix_coefside_t, int, SolverCblk* ) = NULL;

    switch( flttype ) {
    case PastixComplex32:
        compfunc = cpucblk_ccompress;
        break;
    case PastixComplex64:
        compfunc = cpucblk_zcompress;
        break;
    case PastixFloat:
        compfunc = cpucblk_scompress;
        break;
    case PastixDouble:
    case PastixPattern:
    default:
        compfunc = cpucblk_dcompress;
    }

    for (i=0; i < solvmtx->ttsknbr[rank]; i++)
    {
        task     = solvmtx->ttsktab[rank][i];
        itercblk = solvmtx->tasktab[task].cblknum;
        cblk     = solvmtx->cblktab + itercblk;

        if ( cblk->cblktype & CBLK_COMPRESSED ) {
            gain += compfunc( solvmtx, side, -1, cblk );
        }
    }

    pastix_atomic_lock( lock );
    *fullgain += gain;
    pastix_atomic_unlock( lock );
}

/**
 *******************************************************************************
 *
 * @brief Compress the factorized matrix structure if not already done.
 *
 * This routine compress all column blocks that are marked for compression, and
 * return the amount of memory saved by the compression.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that holds the problem
 *
 *******************************************************************************
 *
 * @return The memory gain resulting from the compression to low-rank format in
 *         number of elements.
 *
 *******************************************************************************/
pastix_int_t
coeftabCompress( pastix_data_t *pastix_data )
{
    struct coeftabcomp_s args;
    pastix_lr_t  *lr;

    args.solvmtx = pastix_data->solvmatr;
    args.flttype = pastix_data->bcsc->flttype;
    args.lock    = PASTIX_ATOMIC_UNLOCKED;
    args.gain    = 0;

    /* Set the lowrank properties */
    lr = &(pastix_data->solvmatr->lowrank);
    lr->compress_method     = pastix_data->iparm[IPARM_COMPRESS_METHOD];
    lr->compress_min_width  = pastix_data->iparm[IPARM_COMPRESS_MIN_WIDTH];
    lr->compress_min_height = pastix_data->iparm[IPARM_COMPRESS_MIN_HEIGHT];
    lr->tolerance           = pastix_data->dparm[DPARM_COMPRESS_TOLERANCE];

    isched_parallel_call( pastix_data->isched, pcoeftabComp, (void*)(&args) );

    return args.gain;
}

/**
 *******************************************************************************
 *
 * @brief Compute the ILU levels of a cblk
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix data structure.
 *
 * @param[in] cblk
 *          The column block to compute the ILU levels of.
 *
 *******************************************************************************/
void
coeftabComputeCblkILULevels( const SolverMatrix *solvmtx, SolverCblk *cblk )
{
    /* If there are off diagonal supernodes in the column */
    SolverBlok *blokB = cblk->fblokptr + 1; /* this diagonal block */
    SolverBlok *lblkB = cblk[1].fblokptr;   /* the next diagonal block */

    for (; blokB<lblkB; blokB++) {
        SolverCblk *fcblk = solvmtx->cblktab + blokB->fcblknm;
        SolverBlok *blokC = fcblk->fblokptr;
        SolverBlok *blokA;

        /* If there are off-diagonal supernodes in the column*/
        for (blokA=blokB; blokA<lblkB; blokA++) {
            int lvl_AB;

            /* Find the facing block */
            while ( !is_block_inside_fblock(blokA, blokC) ) {
                blokC++;
                assert( blokC < fcblk[1].fblokptr );
            }

            /* Compute the level k of the block */
            if ( (blokA->iluklvl == INT_MAX) ||
                 (blokB->iluklvl == INT_MAX) )
            {
                lvl_AB = INT_MAX;
            }
            else {
                lvl_AB = blokA->iluklvl + blokB->iluklvl + 1;
            }

            pastix_cblk_lock( fcblk );
            blokC->iluklvl = pastix_imin( blokC->iluklvl,
                                          lvl_AB );
            assert( blokC->iluklvl >= 0 );
            pastix_cblk_unlock( fcblk );
        }

        pastix_atomic_dec_32b( &(fcblk->ctrbcnt) );
    }
}

#if defined(PASTIX_WITH_MPI)
/**
 *******************************************************************************
 *
 * @brief Send all the cblk to their original nodes from the root node.
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The solver matrix to scatter.
 *
 * @param[in] comm
 *          MPI communicator
 *
 * @param[in] root
 *          Root node from which the solvmtx is scattered.
 *
 * @param[in] typesze
 *          Size of the data elements in the SolverMatrix.
 *
 ******************************************************************************/
void
coeftab_scatter( SolverMatrix     *solvmtx,
                 PASTIX_Comm       comm,
                 pastix_int_t      root,
                 pastix_coeftype_t flttype )
{
    static void (*free_fct_array[4])( pastix_coefside_t, SolverCblk * ) = {
        cpucblk_sfree, cpucblk_dfree, cpucblk_cfree, cpucblk_zfree
    };
    void (*free_fct)( pastix_coefside_t, SolverCblk * ) = free_fct_array[ flttype - 2 ];
    SolverCblk  *cblk;
    pastix_int_t i;
    MPI_Status   status;
    size_t eltsize = pastix_size_of( flttype );
    pastix_coefside_t side = PastixLCoef;

    if ( solvmtx->factotype == PastixFactLU ) {
        side = PastixLUCoef;
        eltsize *= 2;
    }

    cblk = solvmtx->cblktab;
    for(i=0; i<solvmtx->cblknbr; i++, cblk++)
    {
        size_t cblksize = eltsize * cblk->stride  * cblk_colnbr( cblk );

        /* Data which does not belong to the root node must be sent */
        if ( (solvmtx->clustnum == root) &&
             (solvmtx->clustnum != cblk->ownerid) )
        {
            MPI_Send( cblk->lcoeftab, cblksize, MPI_CHAR,
                      cblk->ownerid, i, comm );
            free_fct( side, cblk );
        }

        /* Data which belongs locally must be received */
        if ( (solvmtx->clustnum != root) &&
             (solvmtx->clustnum == cblk->ownerid) )
        {
            MPI_Recv( cblk->lcoeftab, cblksize, MPI_CHAR,
                      root, i, comm, &status );
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Gather all the column blocks on the root node.
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The solver matrix to gather.
 *
 * @param[in] comm
 *          MPI communicator
 *
 * @param[in] root
 *          Root node on which the solvmtx is gathered.
 *
 * @param[in] typesze
 *          Size of the data elements in the SolverMatrix.
 *
 ******************************************************************************/
void
coeftab_gather( SolverMatrix     *solvmtx,
                PASTIX_Comm       comm,
                pastix_int_t      root,
                pastix_coeftype_t flttype )
{
    static void (*alloc_fct_array[4])( pastix_coefside_t, SolverCblk * ) = {
        cpucblk_salloc, cpucblk_dalloc, cpucblk_calloc, cpucblk_zalloc
    };
    void (*alloc_fct)( pastix_coefside_t, SolverCblk * ) = alloc_fct_array[ flttype - 2 ];
    SolverCblk       *cblk;
    pastix_int_t      i;
    MPI_Status        status;
    size_t            eltsize = pastix_size_of( flttype );
    pastix_coefside_t side    = PastixLCoef;
    size_t            bufsize = 0;
    void             *recvbuf = NULL;

    if ( solvmtx->factotype == PastixFactLU ) {
        side = PastixLUCoef;
        eltsize *= 2;
    }

    cblk = solvmtx->cblktab;
    for( i=0; i<solvmtx->cblknbr; i++, cblk++ )
    {
        size_t cblksize = eltsize * cblk->stride  * cblk_colnbr( cblk );

        if ( (solvmtx->clustnum == root) &&
             (solvmtx->clustnum != cblk->ownerid) )
        {
            if ( cblk->cblktype & CBLK_COMPRESSED ) {

                cblksize += (cblk[1].fblokptr - cblk[0].fblokptr) * sizeof(int);
                if ( side != PastixLCoef ) {
                    cblksize *= 2;
                }
                MALLOC_INTERN( recvbuf, cblksize, char );
                MPI_Recv( recvbuf, cblksize, MPI_CHAR,
                          cblk->ownerid, i, comm, &status );

                switch ( flttype ) {
                    case PastixComplex64:
                        cpucblk_zunpack_lr( side, cblk, recvbuf );
                        break;
                    case PastixComplex32:
                        cpucblk_cunpack_lr( side, cblk, recvbuf );
                        break;
                    case PastixDouble:
                        cpucblk_dunpack_lr( side, cblk, recvbuf );
                        break;
                    case PastixFloat:
                        cpucblk_sunpack_lr( side, cblk, recvbuf );
                        break;
                    default:
                        assert( 0 );
                }
                free( recvbuf );
            }
            else {
                alloc_fct( side, cblk );
                MPI_Recv( cblk->lcoeftab, cblksize, MPI_CHAR,
                          cblk->ownerid, i, comm, &status );
            }
        }

        if ( (solvmtx->clustnum != root) &&
             (solvmtx->clustnum == cblk->ownerid) )
        {
            if ( cblk->cblktype & CBLK_COMPRESSED ) {
                void *buffer = NULL;

                switch ( flttype ) {
                case PastixComplex64:
                    bufsize = cpucblk_zcompute_size( side, cblk );
                    buffer = cpucblk_zpack_lr( side, cblk, bufsize );
                    break;

                case PastixComplex32:
                    bufsize = cpucblk_ccompute_size( side, cblk );
                    buffer = cpucblk_cpack_lr( side, cblk, bufsize );
                    break;

                case PastixDouble:
                    bufsize = cpucblk_dcompute_size( side, cblk );
                    buffer = cpucblk_dpack_lr( side, cblk, bufsize );
                    break;

                case PastixFloat:
                    bufsize = cpucblk_scompute_size( side, cblk );
                    buffer = cpucblk_spack_lr( side, cblk, bufsize );
                    break;

                default:
                    assert( 0 );
                }

                MPI_Send( buffer, bufsize, MPI_CHAR,
                          root, i, comm );
                free( buffer );
            } else {
                MPI_Send( cblk->lcoeftab, cblksize, MPI_CHAR,
                          root, i, comm );
            }
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Nullify all the distant cblks.
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The solver matrix.
 *
 ******************************************************************************/
void
coeftab_nullify( SolverMatrix *solvmtx )
{
    SolverCblk *cblk;
    pastix_int_t i;
    pastix_coefside_t side = ( solvmtx->factotype == PastixFactLU ) ? PastixLUCoef : PastixLCoef;

    cblk = solvmtx->cblktab;
    for( i=0; i < solvmtx->cblknbr; i++, cblk++ )
    {
        if ( solvmtx->clustnum != cblk->ownerid )
        {
            cpucblk_zfree( side, cblk );
        }
    }
}
#endif /* PASTIX_WITH_MPI */
