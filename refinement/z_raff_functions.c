/**
 *
 * @file z_raff_functions.c
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
#include "bcsc.h"
#include "z_bcsc.h"
#include "sopalin_thread.h"
#include "sopalin_data.h"
#include "z_raff_functions.h"

/**
 *******************************************************************************
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
    fprintf(stdout, OUT_ITERRAFF_ITER, (int)nb_iters);
    fprintf(stdout, OUT_ITERRAFF_TTT, stt);
    fprintf(stdout, OUT_ITERRAFF_ERR, err);
}

/**
 *******************************************************************************
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
 * @param[in] gmres
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

    for (i=0; i<n; i++)
        xptr[i] = gmresx[i];
}

/**
 *******************************************************************************
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

    if (pastix_data->iparm[IPARM_ONLY_RAFF] == API_NO)
    {
        for (i=0; i<n; i++, xptr++)
            gmresx[i]= *xptr;
    }
    else
    {
        for (i=0; i<n; i++, xptr++)
            gmresx[i]=0.0;
    }
}

/**
 *******************************************************************************
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
 * @brief Initiate the B vector used in iterative methods
 *
 *******************************************************************************
 *
 * @param[in] b
 *          The vector given by the user
 *
 * @param[out] raffb
 *          The vector used in iterative methods
 *
 * @param[in] n
 *          The number of elements of both b and raffb
 *
 *******************************************************************************/
void z_Pastix_B( void *b, pastix_complex64_t *raffb, pastix_int_t n )
{
    pastix_complex64_t *bptr = (pastix_complex64_t *)b;
    pastix_int_t i;

    for (i=0; i<n; i++, bptr++)
    {
        raffb[i]= *bptr;
    }
}

/**
 *******************************************************************************
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
pastix_complex64_t z_Pastix_Norm2( pastix_complex64_t *x, pastix_int_t n )
{
    double normx;
    void *xptr = (void*)x;
    normx = z_vectFrobeniusNorm(xptr, n);
    return normx;
}

/**
 *******************************************************************************
 *
 * @brief Apply a preconditionner
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The PaStiX data structure that describes the solver instance.
 *
 * @param[in] s
 *          The vector on which preconditionner is to be applied
 *
 * @param[out] d
 *          On exit, the d vector contains s preconditionned
 *
 *******************************************************************************/
void z_Pastix_Precond( pastix_data_t *pastix_data, pastix_complex64_t *s, pastix_complex64_t *d )
{
    pastix_int_t n = pastix_data->bcsc->gN;
    pastix_int_t nrhs = 1;
    void* bptr = (void*)d;

    memcpy(d, s, n * sizeof( pastix_complex64_t ));
    if (pastix_data->iparm[IPARM_ONLY_RAFF] == API_NO)
    {
        sopalin_data_t sopalin_data;
        sopalin_data.solvmtx = pastix_data->solvmatr;

        switch ( pastix_data->iparm[IPARM_FACTORIZATION] ){
        case PastixFactLLT:
            sopalin_ztrsm( pastix_data, PastixLeft, PastixLower,
                           PastixNoTrans,   PastixNonUnit, &sopalin_data, nrhs, bptr, n );
            sopalin_ztrsm( pastix_data, PastixLeft, PastixLower,
                           PastixConjTrans, PastixNonUnit, &sopalin_data, nrhs, bptr, n );
            break;

        case PastixFactLDLT:
            sopalin_ztrsm( pastix_data, PastixLeft, PastixLower,
                           PastixNoTrans, PastixUnit, &sopalin_data, nrhs, bptr, n );
            sopalin_zdiag( pastix_data, &sopalin_data, nrhs, bptr, n );
            sopalin_ztrsm( pastix_data, PastixLeft, PastixLower,
                           PastixTrans,   PastixUnit, &sopalin_data, nrhs, bptr, n );
            break;

        case PastixFactLDLH:
            sopalin_ztrsm( pastix_data, PastixLeft, PastixLower,
                           PastixNoTrans,   PastixUnit, &sopalin_data, nrhs, bptr, n );
            sopalin_zdiag( pastix_data, &sopalin_data, nrhs, bptr, n );
            sopalin_ztrsm( pastix_data, PastixLeft, PastixLower,
                           PastixConjTrans, PastixUnit, &sopalin_data, nrhs, bptr, n );
            break;

        case PastixFactLU:
        default:
            sopalin_ztrsm( pastix_data, PastixLeft, PastixLower,
                           PastixNoTrans, PastixUnit,    &sopalin_data, nrhs, bptr, n );
            sopalin_ztrsm( pastix_data, PastixLeft, PastixUpper,
                           PastixNoTrans, PastixNonUnit, &sopalin_data, nrhs, bptr, n );
            break;
        }
    }
}

/**
 *******************************************************************************
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
void z_Pastix_Scal( pastix_int_t n, pastix_complex64_t alpha, pastix_complex64_t *x )
{
    z_bcscScal( x, alpha, n, 1);
}

#if defined(PRECISION_z) || defined(PRECISION_c)
/**
 *******************************************************************************
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
void z_Pastix_Dotc( pastix_int_t n, pastix_complex64_t *x,
                    pastix_complex64_t *y, pastix_complex64_t *r )
{
    *r = z_bcscDotc(n, x, y);
}
#endif

/**
 *******************************************************************************
 *
 * @brief Compute a regular scalar product x.y
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
 *******************************************************************************
 *
 * @return The allocated vector
 *
 *******************************************************************************/
void z_Pastix_Dotu( pastix_int_t n, pastix_complex64_t *x,
                    pastix_complex64_t *y, pastix_complex64_t *r )
{
    *r = z_bcscDotu(n, x, y);
}

/**
 *******************************************************************************
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
void z_Pastix_bMAx( pastix_bcsc_t *bcsc, pastix_complex64_t *b,
                    pastix_complex64_t *x, pastix_complex64_t *r )
{
    pastix_int_t alpha = -1.0;
    pastix_int_t beta = 1.0;
    void* xptr = (void*)x;
    void* yptr = (void*)r;

    memcpy(r, b, bcsc->gN * sizeof( pastix_complex64_t ));
    z_bcscGemv(PastixNoTrans, alpha, bcsc, xptr, beta, yptr );
}

/**
 *******************************************************************************
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
    void *yptr = (void*)y;
    void *xptr = (void*)x;

    z_bcscScal( xptr, *beta, n, 1);
    z_bcscAxpy( n, 1., 1., yptr, xptr );
}


/**
 *******************************************************************************
 *
 * @brief Perform y = alpha * coeff x + y
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of elements of vectors x and y
 *
 * @param[in] coeff
 *          The first scaling parameter
 *
 * @param[in] alpha
 *          The second scaling parameter
 *
 * @param[in] x
 *          The vector to be scaled
 *
 * @param[in] y
 *          The resulting solution
 *
 *******************************************************************************/
void z_Pastix_AXPY( pastix_int_t n, double coeff,
                    pastix_complex64_t *alpha, pastix_complex64_t *x,
                    pastix_complex64_t *y )
{
    void *yptr = (void*)y;
    void *xptr = (void*)x;
    z_bcscAxpy( n, 1, coeff*(*alpha), yptr, xptr );
}


/**
 *******************************************************************************
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
    sopthread_data_t *argument = (sopthread_data_t *)arg;
    pastix_int_t        me       = argument->me;
    return me;
}

/**
 *******************************************************************************
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
    solveur->Malloc      = &z_Pastix_Malloc;
    solveur->Free        = &z_Pastix_Free;

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
    solveur->Precond = &z_Pastix_Precond;
    solveur->Scal    = &z_Pastix_Scal;
    solveur->Dotc    = &z_Pastix_Dotc;
    solveur->Ax      = &z_Pastix_Ax;
    solveur->AXPY    = &z_Pastix_AXPY;
    solveur->bMAx    = &z_Pastix_bMAx;
    solveur->BYPX    = &z_Pastix_BYPX;
}
