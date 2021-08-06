/**
 *
 * @file codelet_zgemmsp.c
 *
 * StarPU codelets for blas-like functions
 *
 * @copyright 2016-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Ian Masliah
 * @date 2021-06-21
 *
 * @precisions normal z -> z c d s
 *
 * @addtogroup pastix_starpu
 * @{
 *
 **/
#include "common.h"
#include "blend/solver.h"
#include "sopalin/sopalin_data.h"
#include "pastix_zcores.h"
#if defined(PASTIX_WITH_CUDA)
#include "pastix_zcuda.h"
#endif
#include "pastix_starpu.h"
#include "pastix_zstarpu.h"
#include "codelets.h"
#include "pastix_starpu_model.h"

/**
 * Cblk version
 */

#if defined( PASTIX_STARPU_PROFILING )
measure_t cblk_zgemmsp_perf[STARPU_NMAXWORKERS];
#endif

struct cl_cblk_zgemmsp_args_s {
    sopalin_data_t   *sopalin_data;
#if defined( PASTIX_STARPU_PROFILING )
    measure_t        *perf;
#endif
    double            flops;
    pastix_coefside_t sideA;
    pastix_coefside_t sideB;
    pastix_trans_t    trans;
    const SolverCblk *cblk;
    const SolverBlok *blok;
    SolverCblk       *fcblk;
};

static struct starpu_perfmodel starpu_cblk_zgemmsp_model =
{
    .type = STARPU_PER_ARCH,
    .symbol = "cblk_zgemmsp",
    .arch_cost_function = cblk_gemmsp_cost,
};


#if !defined(PASTIX_STARPU_SIMULATION)
static void fct_cblk_zgemmsp_cpu(void *descr[], void *cl_arg)
{
    
    const pastix_complex64_t *     A;
    const pastix_complex64_t *     B;
    pastix_complex64_t *           C;
    struct cl_cblk_zgemmsp_args_s *args = (struct cl_cblk_zgemmsp_args_s *) cl_arg;

    A = (const pastix_complex64_t *)STARPU_VECTOR_GET_PTR(descr[0]);
    B = (const pastix_complex64_t *)STARPU_VECTOR_GET_PTR(descr[1]);
    C = (pastix_complex64_t *)STARPU_VECTOR_GET_PTR(descr[2]);

    /* Check layout due to NULL workspace for now */
    assert(  args->cblk->cblktype & CBLK_LAYOUT_2D );
    assert( args->fcblk->cblktype & CBLK_LAYOUT_2D );

    args->flops = cpucblk_zgemmsp( args->sideA, args->sideB, args->trans,
                                   args->cblk, args->blok, args->fcblk,
                                   A, B, C, NULL, -1,
                                   &( args->sopalin_data->solvmtx->lowrank ) );
}

#if defined(PASTIX_WITH_CUDA)
static void fct_cblk_zgemmsp_gpu(void *descr[], void *cl_arg)
{
    const cuDoubleComplex         *A;
    const cuDoubleComplex         *B;
    cuDoubleComplex               *C;
    struct cl_cblk_zgemmsp_args_s *args  = (struct cl_cblk_zgemmsp_args_s *) cl_arg;

    A = (const cuDoubleComplex *)STARPU_VECTOR_GET_PTR(descr[0]);
    B = (const cuDoubleComplex *)STARPU_VECTOR_GET_PTR(descr[1]);
    C = (cuDoubleComplex *)STARPU_VECTOR_GET_PTR(descr[2]);

    /* Check layout due to NULL workspace for now */
    assert(  args->cblk->cblktype & CBLK_LAYOUT_2D );
    assert( args->fcblk->cblktype & CBLK_LAYOUT_2D );

    args->flops = gpucblk_zgemmsp( args->sideA, args->sideB, args->trans,
                                   args->cblk, args->blok, args->fcblk,
                                   A, B, C,
                                   &( args->sopalin_data->solvmtx->lowrank ),
                                   starpu_cuda_get_local_stream() );
}
#endif /* defined(PASTIX_WITH_CUDA) */
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_ANY( cblk_zgemmsp, 3, STARPU_CUDA_ASYNC );

void
starpu_task_cblk_zgemmsp( sopalin_data_t   *sopalin_data,
                          pastix_coefside_t sideA,
                          pastix_coefside_t sideB,
                          pastix_trans_t    trans,
                          const SolverCblk *cblk,
                          const SolverBlok *blok,
                          SolverCblk       *fcblk,
                          int               prio )
{
    struct cl_cblk_zgemmsp_args_s *cl_arg;
    struct starpu_task *task;
    int ret;

    /*
     * Create the arguments array
     */
    cl_arg = malloc( sizeof(struct cl_cblk_zgemmsp_args_s) );
    cl_arg->sopalin_data = sopalin_data;
#if defined(PASTIX_STARPU_PROFILING)
    cl_arg->perf         = cblk_zgemmsp_perf;
    cl_arg->flops        = NAN;
#endif
    cl_arg->sideA        = sideA;
    cl_arg->sideB        = sideB;
    cl_arg->trans        = trans;
    cl_arg->cblk         = cblk;
    cl_arg->blok         = blok;
    cl_arg->fcblk        = fcblk;

    /*
     * Create the task structure
     */
    task = starpu_task_create();
    //task->name;
    task->cl = &cl_cblk_zgemmsp_any;

#if defined(PASTIX_WITH_CUDA)
    if ( !(cblk->cblktype  & CBLK_COMPRESSED) &&
         !(fcblk->cblktype & CBLK_COMPRESSED) )
    {
        /* Disable CUDA */
        task->where = cl_cblk_zgemmsp_any.where & (~STARPU_CUDA);
    }
#endif

    /* Register the handles */
    task->handles[0] = cblk->handler[sideA];
    task->handles[1] = cblk->handler[sideB];
    task->handles[2] = fcblk->handler[sideA];

    /* Register appropriate modes of handling data */
    task->cl->modes[0] = STARPU_R;
    task->cl->modes[1] = STARPU_R;
    task->cl->modes[2] = STARPU_RW;

    /* Register the arguments */
    task->cl_arg      = cl_arg;
    task->cl_arg_size = sizeof( struct cl_cblk_zgemmsp_args_s );
    task->cl_arg_free = 1;

#if defined(PASTIX_STARPU_PROFILING)
    /* Register the callback function to profile the performance */
    task->callback_func     = cl_profiling_callback;
    task->callback_arg      = cl_arg;
    task->callback_arg_free = 0; /* Already freed by cl_arg_free */
#endif

#if defined(PASTIX_STARPU_SYNCHRONOUS)
    /* Synchronous call for debug */
    task->synchronous = 1;
#endif

    /* Select the priority */
#if defined(PASTIX_STARPU_HETEROPRIO)
    task->priority = BucketGEMM1D;
    (void) prio;
#else
    task->priority = prio;
#endif

    //task->flops; /* Check if that can be used in the kernels */

    /*
     * Submit the task for execution
     */
    ret = starpu_task_submit(task);
    if (ret == -ENODEV)
    {
        fprintf( stderr, "starpu_task_cblk_zgemmsp: could not submit the task\n" );
        task->destroy = 0;
        starpu_task_destroy(task);
    }

}

/**
 * Blok version
 */
#if defined( PASTIX_STARPU_PROFILING )
measure_t blok_zgemmsp_perf[STARPU_NMAXWORKERS];
#endif


struct cl_blok_zgemmsp_args_s {
    sopalin_data_t   *sopalin_data;
#if defined( PASTIX_STARPU_PROFILING )
    measure_t        *perf;
#endif
    double            flops;
    pastix_coefside_t sideA;
    pastix_coefside_t sideB;
    pastix_trans_t    trans;
    const SolverCblk *cblk;
    SolverCblk       *fcblk;
    pastix_int_t      blok_mk;
    pastix_int_t      blok_nk;
    pastix_int_t      blok_mn;
};

static struct starpu_perfmodel starpu_blok_zgemmsp_model =
{
    .type = STARPU_PER_ARCH,
    .symbol = "blok_zgemmsp",
    .arch_cost_function = blok_gemmsp_cost,
};

#if !defined(PASTIX_STARPU_SIMULATION)
static void fct_blok_zgemmsp_cpu( void *descr[], void *cl_arg )
{
    const pastix_complex64_t      *A;
    const pastix_complex64_t      *B;
    pastix_complex64_t            *C;
    struct cl_blok_zgemmsp_args_s *args = (struct cl_blok_zgemmsp_args_s *) cl_arg;

    A = (const pastix_complex64_t *)STARPU_VECTOR_GET_PTR(descr[0]);
    B = (const pastix_complex64_t *)STARPU_VECTOR_GET_PTR(descr[1]);
    C = (pastix_complex64_t *)STARPU_VECTOR_GET_PTR(descr[2]);

    assert( args->cblk->cblktype  & CBLK_TASKS_2D );
    assert( args->fcblk->cblktype & CBLK_TASKS_2D );

    args->flops = cpublok_zgemmsp( args->sideA, args->sideB, args->trans,
                                   args->cblk, args->fcblk,
                                   args->blok_mk, args->blok_nk, args->blok_mn,
                                   A, B, C,
                                   &(args->sopalin_data->solvmtx->lowrank) );
}

#if defined(PASTIX_WITH_CUDA)
static void fct_blok_zgemmsp_gpu( void *descr[], void *cl_arg )
{
    const cuDoubleComplex         *A;
    const cuDoubleComplex         *B;
    cuDoubleComplex               *C;
    struct cl_blok_zgemmsp_args_s *args = (struct cl_blok_zgemmsp_args_s *) cl_arg;

    A = (const cuDoubleComplex *)STARPU_VECTOR_GET_PTR(descr[0]);
    B = (const cuDoubleComplex *)STARPU_VECTOR_GET_PTR(descr[1]);
    C = (cuDoubleComplex *)STARPU_VECTOR_GET_PTR(descr[2]);

    assert( args->cblk->cblktype  & CBLK_TASKS_2D );
    assert( args->fcblk->cblktype & CBLK_TASKS_2D );

    args->flops = gpublok_zgemmsp( args->sideA, args->sideB, args->trans,
                                   args->cblk, args->fcblk,
                                   args->blok_mk, args->blok_nk, args->blok_mn,
                                   A, B, C,
                                   &(args->sopalin_data->solvmtx->lowrank),
                                   starpu_cuda_get_local_stream() );
}
#endif /* defined(PASTIX_WITH_CUDA) */
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_ANY( blok_zgemmsp, 3, STARPU_CUDA_ASYNC );

void
starpu_task_blok_zgemmsp( sopalin_data_t   *sopalin_data,
                          pastix_coefside_t sideA,
                          pastix_coefside_t sideB,
                          pastix_trans_t    trans,
                          const SolverCblk *cblk,
                          SolverCblk       *fcblk,
                          const SolverBlok *blokA,
                          const SolverBlok *blokB,
                          int               prio )
{
    SolverBlok  *blokC = fcblk->fblokptr;
    pastix_int_t frownum;
    pastix_int_t lrownum;
    pastix_int_t blok_mn = 0, j = 0;
    pastix_int_t blok_mk = blokA - cblk->fblokptr;
    pastix_int_t blok_nk = blokB - cblk->fblokptr;

    struct cl_blok_zgemmsp_args_s *cl_arg;
    struct starpu_task            *task;
    int                            ret;

    assert( blok_nk <= blok_mk );

    do {
        frownum = blokC->frownum;
        lrownum = blokC->lrownum;
        blok_mn += j;
        j = 1;

        /* Increase lrownum as long as blocks are facing the same cblk */
        while( (blokC < fcblk[1].fblokptr-1) &&
               (blokC[0].fcblknm == blokC[1].fcblknm) &&
               (blokC[0].lcblknm == blokC[1].lcblknm) )
        {
            blokC++; j++;
            lrownum = blokC->lrownum;
        }
        blokC++;
    }
    while( !((blokA->frownum >= frownum) &&
             (blokA->lrownum <= lrownum)) );

    blokC = fcblk->fblokptr + blok_mn;

    assert( blokA->lcblknm == blokB->lcblknm );
    assert( blokB->fcblknm == blokC->lcblknm );
    assert( blokC->frownum <= blokA->frownum );
    assert( blokA[-1].fcblknm != blokA[0].fcblknm );
    assert( blokB[-1].fcblknm != blokB[0].fcblknm );
    assert( (blok_mn == 0) || (blokC[-1].fcblknm != blokC[0].fcblknm) );

    /*
     * Create the arguments array
     */
    cl_arg = malloc( sizeof(struct cl_blok_zgemmsp_args_s) );
    cl_arg->sopalin_data = sopalin_data;
#if defined(PASTIX_STARPU_PROFILING)
    cl_arg->perf         = blok_zgemmsp_perf;
#endif
    cl_arg->flops        = NAN;
    cl_arg->sideA        = sideA;
    cl_arg->sideB        = sideB;
    cl_arg->trans        = trans;
    cl_arg->cblk         = cblk;
    cl_arg->fcblk        = fcblk;
    cl_arg->blok_mk      = blok_mk;
    cl_arg->blok_nk      = blok_nk;
    cl_arg->blok_mn      = blok_mn;

    /*
     * Create the task structure
     */
    task = starpu_task_create();
    //task->name;
    task->cl = &cl_blok_zgemmsp_any;

#if defined(PASTIX_WITH_CUDA)
    if ( !(cblk->cblktype  & CBLK_COMPRESSED) &&
         !(fcblk->cblktype & CBLK_COMPRESSED) )
    {
        /* Disable CUDA */
        task->where = cl_blok_zgemmsp_any.where & (~STARPU_CUDA);
    }
#endif

    /* Register the handles */
    task->handles[0] = blokA->handler[sideA];
    task->handles[1] = blokB->handler[sideB];
    task->handles[2] = blokC->handler[sideA];

    /* Register appropriate modes of handling data */
    task->cl->modes[0] = STARPU_R;
    task->cl->modes[1] = STARPU_R;
    task->cl->modes[2] = STARPU_RW;

    /* Register the arguments */
    task->cl_arg = cl_arg;
    task->cl_arg_size = sizeof(cl_arg);
    task->cl_arg_free = 1;

#if defined( PASTIX_STARPU_PROFILING )
    /* Register the callback function to profile the performance */
    task->callback_func     = cl_profiling_callback;
    task->callback_arg      = cl_arg;
    task->callback_arg_free = 0; /* Already freed by cl_arg_free */
#endif

#if defined(PASTIX_STARPU_SYNCHRONOUS)
    /* Synchronous call for debug */
    task->synchronous = 1;
#endif

    /* Select the priority */
#if defined(PASTIX_STARPU_HETEROPRIO)
    task->priority = BucketGEMM2D;
    (void) prio;
#else
    task->priority = prio;
#endif
    /*
     * Submit the task for execution
     */
    ret = starpu_task_submit(task);
    if (ret == -ENODEV)
    {
        fprintf( stderr, "starpu_task_blok_zgemmsp: could not submit the task\n" );
        task->destroy = 0;
        starpu_task_destroy(task);
    }
}

/**
 * @}
 */
