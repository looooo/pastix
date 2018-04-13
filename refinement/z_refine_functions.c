/**
 *
 * @file z_refine_functions.c
 *
 * PaStiX refinement functions implementations.
 *
 * @copyright 2015-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Theophile Terraz
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "cblas.h"
#include "bcsc.h"
#include "z_bcsc.h"
#include "sopalin_data.h"
#include "z_refine_functions.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief Allocate a vector
 *
 *******************************************************************************
 *
 * @param[in] size
 *          The size of the vector
 *
 *******************************************************************************
 *
 * @return The allocated vector
 *
 *******************************************************************************/
void *z_Pastix_Malloc( size_t size )
{
    void *x = NULL;
    MALLOC_INTERN(x, size, char);
    memset(x, 0, size);
    return x;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief Free a vector
 *
 *******************************************************************************
 *
 * @param[inout] x
 *          The vector to be free
 *
 *******************************************************************************/
void z_Pastix_Free( void *x )
{
    memFree_null(x);
}


/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief Print statistics about one iteration
 *
 *******************************************************************************
 *
 * @param[in] t0
 *          The clock value at the beginning of the iteration
 *
 * @param[in] tf
 *          The clock value at the end of the iteration
 *
 * @param[in] err
 *          The backward error after the iteration
 *
 * @param[in] nb_iters
 *          Current number of refinement iterations
 *
 *******************************************************************************/
void z_Pastix_Verbose( double t0, double tf, double err, pastix_int_t nb_iters )
{
    double stt;

    stt = tf - t0;
    fprintf(stdout, OUT_ITERREFINE_ITER, (int)nb_iters);
    fprintf(stdout, OUT_ITERREFINE_TTT, stt);
    fprintf(stdout, OUT_ITERREFINE_ERR, err);
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief Final output
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The PaStiX data structure that describes the solver instance.
 *
 * @param[in] err
 *          The final backward error
 *
 * @param[in] nb_iters
 *          The final number of iterations
 *
 * @param[in] tf
 *          The final clock value
 *
 * @param[inout] x
 *          The vector that is to be overwritten by gmresx
 *
 * @param[in] gmresx
 *          The final solution
 *
 *******************************************************************************/
void z_Pastix_End( pastix_data_t *pastix_data, pastix_complex64_t err,
                   pastix_int_t nb_iters, double tf,
                   void *x, pastix_complex64_t *gmresx )
{
    pastix_complex64_t *xptr = (pastix_complex64_t *)x;
    pastix_int_t        n    = pastix_data->bcsc->gN;
    pastix_int_t i;
    (void)err;
    (void)nb_iters;
    (void)tf;

    for (i=0; i<n; i++) {
        xptr[i] = gmresx[i];
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief Initiate first solution depending on the use of a preconditioner
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The PaStiX data structure that describes the solver instance.
 *
 * @param[in] x
 *          The original solution provided by the user
 *
 * @param[inout] gmresx
 *          The starting point of iterative methods
 *
 *******************************************************************************/
void z_Pastix_X( pastix_data_t *pastix_data, void *x, pastix_complex64_t *gmresx )
{
    pastix_int_t        i;
    pastix_int_t        n = pastix_data->bcsc->gN;
    pastix_complex64_t *xptr = (pastix_complex64_t *)x;

    if (1 /*pastix_data->iparm[IPARM_ONLY_REFINE] == 0*/)
    {
        for (i=0; i<n; i++, xptr++) {
            gmresx[i]= *xptr;
        }
    }
    else
    {
        for (i=0; i<n; i++, xptr++) {
            gmresx[i] = 0.0;
        }
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief The number of unknowns
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The PaStiX data structure that describes the solver instance.
 *
 *******************************************************************************
 *
 * @return The number of unknowns
 *
 *******************************************************************************/
pastix_int_t z_Pastix_n( pastix_data_t *pastix_data )
{
    return pastix_data->bcsc->gN;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief Initiate the B vector used in iterative methods
 *
 *******************************************************************************
 *
 * @param[in] b
 *          The vector given by the user
 *
 * @param[out] refineb
 *          The vector used in iterative methods
 *
 * @param[in] n
 *          The number of elements of both b and refineb
 *
 *******************************************************************************/
void z_Pastix_B( const pastix_complex64_t *b, pastix_complex64_t *refineb, pastix_int_t n )
{
    memcpy( refineb, b, n * sizeof(pastix_complex64_t) );
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief The precision required
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The PaStiX data structure that describes the solver instance.
 *
 *******************************************************************************
 *
 * @return The precision required by the user
 *
 *******************************************************************************/
pastix_complex64_t z_Pastix_Eps( pastix_data_t *pastix_data )
{
    return pastix_data->dparm[DPARM_EPSILON_REFINEMENT];
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief The maximum number of iterations
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The PaStiX data structure that describes the solver instance.
 *
 *******************************************************************************
 *
 * @return The maximal number of iterations authorized
 *
 *******************************************************************************/
pastix_int_t z_Pastix_Itermax( pastix_data_t *pastix_data )
{
    return pastix_data->iparm[IPARM_ITERMAX];
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief The size of the Krylov space
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The PaStiX data structure that describes the solver instance.
 *
 *******************************************************************************
 *
 * @return The maximum size authorized for the Krylov space
 *
 *******************************************************************************/
pastix_int_t z_Pastix_Krylov_Space( pastix_data_t *pastix_data )
{
    return pastix_data->iparm[IPARM_GMRES_IM];
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief Perform the norm2 of a vector
 *
 *******************************************************************************
 *
 * @param[in] x
 *          The vector which norm2 is to be computed
 *
 * @param[in] n
 *          The number of elements of x
 *
 *******************************************************************************
 *
 * @return The frobenius norm of the vector
 *
 *******************************************************************************/
double z_Pastix_Norm2( pastix_int_t n, const pastix_complex64_t *x )
{
    double normx;
    normx = z_vectFrobeniusNorm( n, x );
    return normx;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief Solve A x = b with A the sparse matrix
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The PaStiX data structure that describes the solver instance.
 *
 * @param[inout] d
 *          On entry, the right hand side
 *          On exit, the solution of tha problem A x = b
 *
 *******************************************************************************/
void z_Pastix_spsv( pastix_data_t *pastix_data, pastix_complex64_t *b )
{
    pastix_int_t n = pastix_data->bcsc->gN;
    pastix_data->iparm[IPARM_VERBOSE]--;
    pastix_subtask_solve( pastix_data, 1, b, n );
    pastix_data->iparm[IPARM_VERBOSE]++;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief Compute y = \alpha A x + \beta y
 *
 *******************************************************************************
 *
 * @param[in] m
 *          The number of rows of the matrix A, and the size of y.
 *
 * @param[in] n
 *          The number of columns of the matrix A, and the size of x.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          The dense matrix A of size lda-by-n.
 *
 * @param[in] lda
 *          The leading dimension of the matrix A. lda >= max(1,m)
 *
 * @param[in] x
 *          The vector x of size n.
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[inout] y
 *          On entry, the initial vector y of size m.
 *          On exit, the updated vector.
 *
 *******************************************************************************/
void z_Pastix_gemv( pastix_int_t m,
                    pastix_int_t n,
                    pastix_complex64_t alpha,
                    const pastix_complex64_t *A,
                    pastix_int_t lda,
                    const pastix_complex64_t *x,
                    pastix_complex64_t  beta,
                    pastix_complex64_t *y )
{
    cblas_zgemv( CblasColMajor, CblasNoTrans, m, n,
                 CBLAS_SADDR(alpha), A, lda, x, 1,
                 CBLAS_SADDR(beta), y, 1 );
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief Scale a vector
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of elements of the vector
 *
 * @param[in] alpha
 *          The scaling parameter
 *
 * @param[inout] x
 *          The vector to be scaled
 *
 *******************************************************************************/
void z_Pastix_scal( pastix_int_t n, pastix_complex64_t alpha, pastix_complex64_t *x )
{
    cblas_zscal( n, CBLAS_SADDR(alpha), x, 1 );
}

#if defined(PRECISION_z) || defined(PRECISION_c)
/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief Compute a scalar product between complex vectors: x.conj(y)
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of elements of vectors x, y and r
 *
 * @param[in] y
 *          The first vector of the scalar product
 *
 * @param[in] n
 *          The second vector of the scalar product
 *
 * @param[out] r
 *          The result of the scalar product
 *
 *******************************************************************************/
pastix_complex64_t
z_Pastix_Dotu( pastix_int_t n,
               const pastix_complex64_t *x,
               const pastix_complex64_t *y )
{
    return z_bcscDotu( n, x, y );
}
#endif

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief Compute a regular scalar product x.y
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of elements of vectors x, y and r
 *
 * @param[in] x
 *          The first vector of the scalar product
 *
 * @param[in] y
 *          The second vector of the scalar product
 *
 * @param[out] r
 *          The result of the scalar product
 *
 *******************************************************************************
 *
 * @return The allocated vector
 *
 *******************************************************************************/
pastix_complex64_t
z_Pastix_Dotc( pastix_int_t n,
               const pastix_complex64_t *x,
               const pastix_complex64_t *y )
{
    return z_bcscDotc( n, x, y );
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief Perform r = Ax
 *
 *******************************************************************************
 *
 * @param[in] bcsc
 *          The Pastix bcsc
 *
 * @param[in] x
 *          The vector that multiplies the matrix A
 *
 * @param[in] r
 *          The result of the matrix-vector product
 *
 *******************************************************************************/
void z_Pastix_Ax( pastix_bcsc_t *bcsc, pastix_complex64_t *x, pastix_complex64_t *r )
{
    pastix_int_t alpha = 1.0;
    pastix_int_t beta = 0.0;
    void* xptr = (void*)x;
    void* yptr = (void*)r;

    z_bcscGemv(PastixNoTrans, alpha, bcsc, xptr, beta, yptr );
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief Perform r = b - Ax
 *
 *******************************************************************************
 *
 * @param[in] bcsc
 *          The Pastix bcsc
 *
 * @param[in] b
 *          The vector to be copied in r
 *
 * @param[in] x
 *          The vector that multiplies A
 *
 * @param[in] r
 *          The result b-Ax
 *
 *******************************************************************************/
void z_Pastix_bMAx( pastix_bcsc_t            *bcsc,
                    const pastix_complex64_t *b,
                    const pastix_complex64_t *x,
                    pastix_complex64_t       *r )
{
    pastix_complex64_t alpha = -1.0;
    pastix_complex64_t beta = 1.0;

    memcpy( r, b, bcsc->gN * sizeof( pastix_complex64_t ) );
    z_bcscGemv( PastixNoTrans, alpha, bcsc, x, beta, r );
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief Perform y = alpha A x + y
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure that holds the A matrix.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] x
 *          The vector x
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[inout] y
 *          On entry, the vector y
 *          On exit, alpha A x + y
 *
 *******************************************************************************/
void z_Pastix_spmv( pastix_data_t            *pastix_data,
                    pastix_complex64_t        alpha,
                    const pastix_complex64_t *x,
                    pastix_complex64_t        beta,
                    pastix_complex64_t       *y )
{
    pastix_bcsc_t *bcsc = pastix_data->bcsc;
    z_bcscGemv( PastixNoTrans, alpha, bcsc, x, beta, y );
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief Perform x = beta x + y
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of elements of x and y
 *
 * @param[in] beta
 *          The scaling parameter of x
 *
 * @param[in] y
 *          The vector added to beta x
 *
 * @param[inout] x
 *          The resulting vector
 *
 *******************************************************************************/
void z_Pastix_BYPX( pastix_int_t n, pastix_complex64_t *beta,
                    pastix_complex64_t *y, pastix_complex64_t *x )
{
    z_bcscScal( x, *beta, n, 1);
    z_bcscAxpy( n, 1, 1., y, x );
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief Coy a vector y = x
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of elements of vectors x and y
 *
 * @param[in] x
 *          The vector to be scaled
 *
 * @param[inout] y
 *          The resulting solution
 *
 *******************************************************************************/
static inline void
z_Pastix_copy( pastix_int_t              n,
               const pastix_complex64_t *x,
               pastix_complex64_t       *y )
{
    memcpy( y, x, n * sizeof(pastix_complex64_t) );
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief Perform y = alpha * x + y
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of elements of vectors x and y
 *
 * @param[in] alpha
 *          The scalar to scale x
 *
 * @param[in] x
 *          The vector to be scaled
 *
 * @param[inout] y
 *          The resulting solution
 *
 *******************************************************************************/
static inline void
z_Pastix_axpy( pastix_int_t              n,
               pastix_complex64_t        alpha,
               const pastix_complex64_t *x,
               pastix_complex64_t       *y )
{
    z_bcscAxpy( n, 1, alpha, x, y );
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief Rank of the thread
 *
 *******************************************************************************
 *
 * @param[in] arg
 *          The structure containing threads informations
 *
 *******************************************************************************
 *
 * @return The rank of the current thread
 *
 *******************************************************************************/
pastix_int_t z_Pastix_me( void *arg )
{
    (void)arg;
    //sopthread_data_t *argument = (sopthread_data_t *)arg;
    pastix_int_t        me       = 0; //argument->me;
    return me;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief Initiate functions pointers to define basic operations
 *
 *******************************************************************************
 *
 * @param[out] solveur
 *          The structure to be filled
 *
 *******************************************************************************/
void z_Pastix_Solveur( struct z_solver *solveur )
{
    /* Allocations */
    solveur->Malloc  = &z_Pastix_Malloc;
    solveur->Free    = &z_Pastix_Free;

    /* Interface functions */
    solveur->Verbose = &z_Pastix_Verbose;
    solveur->End     = &z_Pastix_End;
    solveur->X       = &z_Pastix_X;
    solveur->N       = &z_Pastix_n;
    solveur->B       = &z_Pastix_B;
    solveur->Eps     = &z_Pastix_Eps;
    solveur->Itermax = &z_Pastix_Itermax;
    solveur->me      = &z_Pastix_me;
    solveur->Krylov_Space = &z_Pastix_Krylov_Space;

    /* Basic operations */
    solveur->Norm    = &z_Pastix_Norm2;
    solveur->Precond = &z_Pastix_spsv;
    solveur->Scal    = &z_Pastix_scal;
    solveur->Ax      = &z_Pastix_Ax;
    solveur->bMAx    = &z_Pastix_bMAx;
    solveur->BYPX    = &z_Pastix_BYPX;
    solveur->AXPY    = &z_Pastix_axpy;

    solveur->output_oneiter = &z_Pastix_Verbose;
    solveur->dot     = &z_Pastix_Dotc;
    solveur->scal    = &z_Pastix_scal;
    solveur->copy    = &z_Pastix_copy;
    solveur->axpy    = &z_Pastix_axpy;
    solveur->spmv    = &z_Pastix_spmv;
    solveur->spsv    = &z_Pastix_spsv;
    solveur->norm    = &z_Pastix_Norm2;
    solveur->gemv    = &z_Pastix_gemv;
}
