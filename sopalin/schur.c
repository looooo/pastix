/**
 *
 * @file sopalin/schur.c
 *
 * PaStiX schur interface functions
 *
 * @copyright 2017-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2021-01-03
 *
 * @addtogroup pastix_schur
 * @{
 *
 **/
#include "common.h"
#include <spm.h>
#include <lapacke.h>
#include "blend/solver.h"
#include "sopalin/coeftab_z.h"
#include "sopalin/coeftab_c.h"
#include "sopalin/coeftab_d.h"
#include "sopalin/coeftab_s.h"

/**
 *******************************************************************************
 *
 * @brief Set the list of unknowns that belongs to the schur complement.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix data structure of the solver to store the list of Schur
 *          unknowns.
 *
 * @param[in] n
 *          The number of unknowns in the Schur complement.
 *
 * @param[in] list
 *          Array of integer of size n.
 *          The list of unknowns belonging to the Schur complement with the same
 *          baseval as the associated spm.
 *
 *******************************************************************************/
void
pastixSetSchurUnknownList( pastix_data_t      *pastix_data,
                           pastix_int_t        n,
                           const pastix_int_t *list)
{
    if ( n > 0 ) {
        pastix_data->schur_n    = n;
        pastix_data->schur_list = (pastix_int_t*)malloc(n * sizeof(pastix_int_t));
        memcpy( pastix_data->schur_list, list, n * sizeof(pastix_int_t) );
    }
}

/**
 *******************************************************************************
 *
 * @brief Return the Schur complement.
 *
 * The Schur complement is returned in the column major layout used by the
 * classic linear algebra libraries such as Blas or Lapack.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix data structure of the problem solved.
 *
 * @param[inout] S
 *          Array of size spm->n -by- lds of arithmetic spm->flttype, where spm
 *          is the spm of the original problem.
 *          On exit, the array contains the Schur complement of the factorized
 *          matrix.
 *
 * @param[in] lds
 *          The leading dimension of the S array.
 *
 ********************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
pastixGetSchur( const pastix_data_t *pastix_data,
                void                *S,
                pastix_int_t         lds )
{
    pastix_int_t *iparm;

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        pastix_print_error( "pastix_getSchur: wrong pastix_data parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }
    if (S == NULL) {
        pastix_print_error( "pastix_getSchur: S parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }
    if (lds <= 0) {
        pastix_print_error( "pastix_getSchur: lds parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_NUMFACT) ) {
        pastix_print_error( "pastix_getSchur: All steps from pastix_task_init() to pastix_task_numfact() have to be called before calling this function" );
        return PASTIX_ERR_BADPARAMETER;
    }
#if defined(PASTIX_WITH_MPI)
    if (pastix_data->inter_node_procnbr > 1) {
        if ( pastix_data->inter_node_procnum == 0 ) {
            pastix_print_error( "pastix_getSchur: Schur complement is not available yet with multiple MPI processes\n" );
        }
        return -1;
    }
#endif

    iparm = pastix_data->iparm;
    switch(iparm[IPARM_FLOAT])
    {
    case PastixPattern:
        break;
    case PastixFloat:
        coeftab_sgetschur( pastix_data->solvmatr, S, lds );
        break;
    case PastixComplex32:
        coeftab_cgetschur( pastix_data->solvmatr, S, lds );
        break;
    case PastixComplex64:
        coeftab_zgetschur( pastix_data->solvmatr, S, lds );
        break;
    case PastixDouble:
    default:
        coeftab_dgetschur( pastix_data->solvmatr, S, lds );
    }
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_solve
 *
 * @brief Get the vector in an RHS data structure.
 *
 *******************************************************************************
 *
 * @param[in] rhsB
 *          The pastix_rhs_t data structure used to solve the system.
 *
 * @param[in] m
 *          The number of rows of the vector b, must be equal to the number of
 *          unknowns in the Schur complement.
 *
 * @param[in] n
 *          The number of columns of the vector b.
 *
 * @param[inout] b
 *          On entry, a vector of size ldb-by-n.
 *          On exit, the m-by-n leading part contains the right hand side
 *          related to the Schur part.
 *
 * @param[in] ldb
 *          The leading dimension of the vector b.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
pastixRhsSchurGet( const pastix_data_t *pastix_data,
                   pastix_int_t         m,
                   pastix_int_t         n,
                   pastix_rhs_t         rhsB,
                   void                *B,
                   pastix_int_t         ldb )
{
    const SolverMatrix *solvmtx;
    const SolverCblk   *cblk;
    pastix_int_t        mschur;
    char               *bptr;

    if ( pastix_data == NULL ) {
        pastix_print_error( "pastixRhsSchurGet: wrong pastix_data parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( rhsB == NULL ) {
        pastix_print_error( "pastixRhsSchurGet: wrong rhsB parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( B == NULL ) {
        pastix_print_error( "pastixRhsSchurGet: wrong b parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }

    solvmtx = pastix_data->solvmatr;
    cblk    = solvmtx->cblktab + solvmtx->cblkschur;
    mschur  = solvmtx->nodenbr - cblk->fcolnum;

    if ( m != mschur ) {
        pastix_print_error( "pastixRhsSchurGet: wrong m parameter expecting %ld but was %ld\n",
                            (long)mschur, (long)m );
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( n != rhsB->n ) {
        pastix_print_error( "pastixRhsSchurGet: wrong n parameter expecting %ld but was %ld\n",
                            (long)rhsB->n, (long)n );
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( ldb < m ) {
        pastix_print_error( "pastixRhsSchurGet: wrong ldb parameter\n" );
        return PASTIX_ERR_BADPARAMETER;
    }

    bptr = rhsB->b;
    bptr += cblk->lcolidx * pastix_size_of( rhsB->flttype );

    switch( rhsB->flttype ) {
    case SpmComplex64:
        LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', mschur, n, (pastix_complex64_t *)bptr, rhsB->ld, B, ldb );
        break;
    case SpmComplex32:
        LAPACKE_clacpy_work( LAPACK_COL_MAJOR, 'A', mschur, n, (pastix_complex32_t *)bptr, rhsB->ld, B, ldb );
        break;
    case SpmDouble:
        LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', mschur, n, (double *)bptr, rhsB->ld, B, ldb );
        break;
    case SpmFloat:
        LAPACKE_slacpy_work( LAPACK_COL_MAJOR, 'A', mschur, n, (float *)bptr, rhsB->ld, B, ldb );
        break;
    default:
        pastix_print_error( "pastixRhsSchurGet: unknown flttype\n" );
        return PASTIX_ERR_BADPARAMETER;
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_solve
 *
 * @brief Set the vector in an RHS data structure.
 *
 *******************************************************************************
 *
 * @param[in] m
 *          The number of rows of the vector b.
 *
 * @param[in] n
 *          The number of columns of the vector b.
 *
 * @param[in] b
 *          The vector b.
 *
 * @param[in] ldb
 *          The leading dimension of the vector b.
 *
 * @param[out] rhsB
 *          The pastix_rhs_t data structure which contains the vector b.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
pastixRhsSchurSet( const pastix_data_t *pastix_data,
                   pastix_int_t         m,
                   pastix_int_t         n,
                   void                *B,
                   pastix_int_t         ldb,
                   pastix_rhs_t         rhsB )
{
    const SolverMatrix *solvmtx;
    const SolverCblk   *cblk;
    pastix_int_t        mschur;
    char               *bptr;

    if ( pastix_data == NULL ) {
        pastix_print_error( "pastixRhsSchurSet: wrong pastix_data parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( rhsB == NULL ) {
        pastix_print_error( "pastixRhsSchurSet: wrong rhsB parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( B == NULL ) {
        pastix_print_error( "pastixRhsSchurSet: wrong b parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }

    solvmtx = pastix_data->solvmatr;
    cblk    = solvmtx->cblktab + solvmtx->cblkschur;
    mschur  = solvmtx->nodenbr - cblk->fcolnum;

    if ( m != mschur ) {
        pastix_print_error( "pastixRhsSchurSet: wrong m parameter expecting %ld but was %ld\n",
                            (long)mschur, (long)m );
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( n != rhsB->n ) {
        pastix_print_error( "pastixRhsSchurSet: wrong n parameter expecting %ld but was %ld\n",
                            (long)rhsB->n, (long)n );
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( ldb < m ) {
        pastix_print_error( "pastixRhsSchurSet: wrong ldb parameter\n" );
        return PASTIX_ERR_BADPARAMETER;
    }

    bptr = rhsB->b;
    bptr += cblk->lcolidx * pastix_size_of( rhsB->flttype );

    switch( rhsB->flttype ) {
    case SpmComplex64:
        LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', mschur, n, B, ldb, (pastix_complex64_t *)bptr, rhsB->ld );
        break;
    case SpmComplex32:
        LAPACKE_clacpy_work( LAPACK_COL_MAJOR, 'A', mschur, n, B, ldb, (pastix_complex32_t *)bptr, rhsB->ld );
        break;
    case SpmDouble:
        LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', mschur, n, B, ldb, (double *)bptr, rhsB->ld );
        break;
    case SpmFloat:
        LAPACKE_slacpy_work( LAPACK_COL_MAJOR, 'A', mschur, n, B, ldb, (float *)bptr, rhsB->ld );
        break;
    default:
        pastix_print_error( "pastixRhsSchurSet: unknown flttype\n" );
        return PASTIX_ERR_BADPARAMETER;
    }

    return PASTIX_SUCCESS;
}

/**
 * @}
 */
