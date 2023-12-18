/**
 *
 * @file pastix_starpu.h
 *
 * StarPU support for the numerical factorization and solve of PaStiX.
 *
 * @copyright 2016-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Alycia Lisito
 * @author Florent Pruvost
 * @author Nolan Bredel
 * @author Tom Moenne-Loccoz
 * @date 2023-12-18
 *
 * @addtogroup pastix_starpu
 * @{
 *   This module describes the functionnality provided by the runtime system
 *   StarPU for the numerical factorization and solve.
 *
 **/
#ifndef _pastix_starpu_h_
#define _pastix_starpu_h_

#include "common.h"
#include "blend/solver.h"

#if defined(PASTIX_WITH_MPI)
#include <starpu_mpi.h>
#else
#include <starpu.h>
#endif

#include <starpu_profiling.h>

#if defined(PASTIX_WITH_CUDA) && !defined(PASTIX_STARPU_SIMULATION)
#include <starpu_scheduler.h>
#include <starpu_cuda.h>

#include <cublas.h>
#include <starpu_cublas.h>
#if defined(PASTIX_WITH_CUBLAS_V2)
#include <cublas_v2.h>
#include <starpu_cublas_v2.h>
#endif
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
typedef struct starpu_conf starpu_conf_t;

#if defined(PASTIX_WITH_MPI)

#if defined(PASTIX_STARPU_SYNC)
#define pastix_starpu_insert_task( _codelet_, ... )                         \
    starpu_mpi_insert_task( sopalin_data->solvmtx->solv_comm, _codelet_, STARPU_TASK_SYNCHRONOUS, 1, ##__VA_ARGS__ )
#else
#define pastix_starpu_insert_task( _codelet_, ... )                         \
    starpu_mpi_insert_task( sopalin_data->solvmtx->solv_comm, _codelet_, ##__VA_ARGS__ )
#endif

#else

#if defined(PASTIX_STARPU_SYNC)
#define pastix_starpu_insert_task( _codelet_, ... )                         \
    starpu_insert_task( _codelet_, STARPU_TASK_SYNCHRONOUS, 1, ##__VA_ARGS__ )
#else
#define pastix_starpu_insert_task( _codelet_, ... )                         \
    starpu_insert_task( _codelet_, ##__VA_ARGS__ )
#endif

#endif

#if defined( PASTIX_STARPU_HETEROPRIO )
typedef enum heteroprio_bucket_order_e {
    BucketSolveDiag = 0,
    BucketSolveGEMM = 0,
    BucketSolveTRSM = 0,
    BucketFacto1D   = 0,
    BucketFacto2D   = 0,
    BucketScalo     = 0,
    BucketTRSM1D    = 2,
    BucketTRSM2D    = 2,
    BucketGEMM1D    = 1,
    BucketGEMM2D    = 3,
    BucketNumber
} heteroprio_bucket_order_t;
#endif
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 * @brief Additional StarPU handlers for a column-block when using 2D kernels.
 *
 * Handle requirements for contiguous allocation of the block handlers when
 * using StarPU data partitioning.
 */
typedef struct starpu_cblk_s {
    pastix_int_t          handlenbr; /**< Number of 2D block handlers in the column-block */
    starpu_data_handle_t *handletab; /**< Array of 2D block handlers for the column-block */
} starpu_cblk_t;

/**
 * @brief StarPU descriptor stucture for the sparse matrix.
 */
typedef struct starpu_sparse_matrix_desc_s {
    int64_t          mpitag;         /**< MPI id of StarPU */
    int              typesze;        /**< Arithmetic size                                                                              */
    pastix_mtxtype_t mtxtype;        /**< Matrix structure: PastixGeneral, PastixSymmetric or PastixHermitian.                         */
    SolverMatrix    *solvmtx;        /**< Solver matrix structure that describes the problem and stores the original data              */
    starpu_cblk_t   *cblktab_handle; /**< Array of 2D column-block handlers (NULL when using 1D kernels only)                          */
    void           **gpu_blocktab;   /**< Pointer to GPU arrays that contains frownum,lrownum of each block for Fermi (NULL otherwise) */
} starpu_sparse_matrix_desc_t;

/**
 * @brief StarPU descriptor for the vectors linked to a given sparse matrix.
 */
typedef struct starpu_rhs_desc_s {
    int64_t               mpitag;    /**< MPI id of StarPU */
    int                   ncol;      /**< Number of columns of the matrix                                                 */
    int                   typesze;   /**< Arithmetic size                                                                 */
    SolverMatrix         *solvmtx;   /**< Solver matrix structure that describes the problem and stores the original data */
    starpu_data_handle_t *handletab; /**< Array of handlers for the blocks */
    void                 *dataptr;   /**< Store the main data pointer to check that the descriptor matches the reference  */
} starpu_rhs_desc_t;

void starpu_sparse_matrix_init( SolverMatrix      *solvmtx,
                                pastix_mtxtype_t   mtxtype,
                                int                nodes,
                                int                myrank,
                                pastix_coeftype_t  flttype );
void starpu_sparse_matrix_destroy( starpu_sparse_matrix_desc_t *desc );
void starpu_sparse_matrix_getoncpu( starpu_sparse_matrix_desc_t *desc );

void starpu_rhs_init( SolverMatrix *solvmtx,
                      pastix_rhs_t  rhsb,
                      int           typesze,
                      int           nodes,
                      int           myrank );
void starpu_rhs_destroy( starpu_rhs_desc_t *desc );
void starpu_rhs_getoncpu( starpu_rhs_desc_t *desc );

void starpu_sparse_cblk_wont_use( pastix_coefside_t coef, SolverCblk *cblk );
void pastix_starpu_init( pastix_data_t *pastix,
                         int           *argc,
                         char         **argv[],
                         const int     *bindtab );
void pastix_starpu_finalize( pastix_data_t *pastix );

/**
 * MPI tag management
 */
int     pastix_starpu_tag_init( pastix_data_t *pastix );
int64_t pastix_starpu_tag_book( int64_t nbtags );
void    pastix_starpu_tag_release( int64_t min );

#ifndef DOXYGEN_SHOULD_SKIP_THIS
struct measure_s;
typedef struct measure_s measure_t;

struct measure_s {
    double sum;
    double sum2;
    long   n;
};

/**
 * @brief Helper function and variable for the testings
 */
struct starpu_profile_s;
typedef struct starpu_profile_s starpu_profile_t;

/**
 * @brief Profiling data structure to register a codelet to profile
 */
struct starpu_profile_s {
    starpu_profile_t *next;                         /**< Link to the next implementation  */
    const char       *name;                         /**< Short name of the function       */
    measure_t         measures[STARPU_NMAXWORKERS]; /**< Pointer to the array of measures */
};
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 * @brief Base structure to all codelet arguments that include the profiling data
 */
typedef struct profile_data_s {
#if defined( PASTIX_STARPU_PROFILING )
    measure_t *measures;
#endif
    double     flops;
} profile_data_t;

#if defined( PASTIX_STARPU_PROFILING )
void cl_profiling_callback( void *callback_arg );
void profiling_register_cl( starpu_profile_t *codelet );
void profiling_display_allinfo();
#else
/**
 *******************************************************************************
 *
 * @brief Displays all profiling data collected into all measurements tables of
 * the profile_list.
 *
 ******************************************************************************/
static inline void profiling_display_allinfo() {}
#endif

#if defined( PASTIX_STARPU_PROFILING_LOG )
void profiling_log_init( const char* dirname );
void cl_profiling_log_register( const char *task_name, const char* cl_name,
                                int m, int n, int k, double flops, double speed );

void profiling_log_fini();
#else
#ifndef DOXYGEN_SHOULD_SKIP_THIS
static inline void profiling_log_init( const char* dirname ) {
    (void) dirname;
}
static inline void profiling_log_fini() {}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#ifdef PASTIX_STARPU_STATS
static inline void
print_stats( double sub, double com, __attribute__((unused)) SolverMatrix *solvmtx )
{
    int src = 0;
    MPI_Comm_rank( solvmtx->solv_comm, &src );
    fprintf( stderr, "    Time to submit tasks on node %d        %e s\n", src, clockVal( sub ) );
    fprintf( stderr, "    Time to execute tasks on node %d       %e s\n", src, clockVal( com ) );
    fprintf( stderr, "    Total time on node %d                  %e s\n", src, clockVal( sub ) + clockVal( com ) );
}
#endif
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 * @brief StarPU Interface to handle cblks and bloks
 */
extern struct starpu_data_interface_ops pastix_starpu_interface_ops;

/**
 * @brief Alias to get the Interface id
 */
#define PASTIX_STARPU_INTERFACE_ID pastix_starpu_interface_ops.interfaceid

/**
 * @brief Interface data structure to register the pieces of data in StarPU
 */
typedef struct pastix_starpu_interface_s {
    enum starpu_data_interface_id id;        /**< Identifier of the interface               */
    pastix_coeftype_t             flttype;   /**< Floating type of the elements             */
    int                           offset;    /**< -1 for cblk, blok offset for the subdatas */
    int                           nbblok;    /**< Number of blocks                          */
    size_t                        allocsize; /**< size currently allocated                  */
    const SolverCblk             *cblk;      /**< Internal structure used to store the cblk */
    void                         *dataptr;   /**< Pointer on data                           */
} pastix_starpu_interface_t;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
static inline void *
pastix_starpu_cblk_get_ptr( void *interf ) {
    return ((pastix_starpu_interface_t *)interf)->dataptr;
}

static inline void *
pastix_starpu_blok_get_ptr( void *interf ) {
    return ((pastix_starpu_interface_t *)interf)->dataptr;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 * @brief Register a cblk at the StarPU level
 *
 * @param[out] handleptr
 *      The StarPU data handle to the registered data. Space must be allocated on call.
 *
 * @param[in] home_node
 *      The StarPU memory node enum to specify where the initial data is located
 *      -1 if not local, STARPU_MAIN_RAM if local.
 *
 * @param[in] cblk
 *      The cblk to register
 *
 * @param[in] side
 *      Specify which part of the cblk (Upper or Lower) to register
 *
 * @param[in] flttype
 *      Specify the arithmetic floating type of the coefficients
 */
void pastix_starpu_register( starpu_data_handle_t *handleptr,
                             const SolverCblk     *cblk,
                             pastix_coefside_t     side,
                             pastix_coeftype_t     flttype );

void pastix_starpu_register_ws( starpu_data_handle_t *handleptr,
                                const SolverCblk     *cblk,
                                pastix_coeftype_t     flttype );

void pastix_starpu_register_blok( starpu_data_handle_t *handleptr,
                                  const SolverCblk     *cblk,
                                  const SolverBlok     *blok,
                                  pastix_coeftype_t     flttype );

/**
 * @brief Initialize the interface ID
 */
void pastix_starpu_interface_init();

/**
 * @brief Finalize the interface and reset the ID
 */
void pastix_starpu_interface_fini();

/**
 * @brief Main structure for all tasks of fanin_init type
 */
extern struct starpu_codelet cl_fanin_init_cpu;

/**
 * @brief Main structure for all tasks of rhs_init type
 */
extern struct starpu_codelet cl_rhs_init_cpu;

#endif /* _pastix_starpu_h_ */

/**
 * @}
 */
