/**
 *
 * @file workspace_z.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Azzam Haidar
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "workspace.h"

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_zgeev - Allocates workspace for PLASMA_zgeev or PLASMA_zgeev_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The order of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile Hessenberg.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_zgeev(int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(N, N, PLASMA_FUNC_ZGEEV, PlasmaComplexDouble, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_zgehrd - Allocates workspace for PLASMA_zgehrd or PLASMA_zgehrd_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile Hessenberg.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_zgehrd(int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(N, N, PLASMA_FUNC_ZGEHRD, PlasmaComplexDouble, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_zgebrd - Allocates workspace for PLASMA_zgebrd or PLASMA_zgebrd_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile BRD.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_zgebrd(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_ZGEBRD, PlasmaComplexDouble, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_zgels - Allocates workspace for PLASMA_zgels or PLASMA_zgels_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile QR
 *          or the tile LQ factorization.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_zgels(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_ZGELS, PlasmaComplexDouble, descT); }

/***************************************************************************//**
@@ -192,33 +120,6 @@
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_zgels_Tile - Allocates tile workspace for PLASMA_zgels_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, PLASMA_desc on workspace handle for storage of the extra T factors required by the tile QR
 *          or the tile LQ factorization.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_zgels_Tile(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_ZGELS, PlasmaComplexDouble, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_zgeqrf - Allocates workspace for PLASMA_zgeqrf or PLASMA_zgeqrf_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile QR
 *          factorization.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_zgeqrf(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_ZGELS, PlasmaComplexDouble, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_zgeqrf_Tile - Allocates tile workspace for PLASMA_zgels_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, PLASMA_desc on workspace handle for storage of the extra T factors required by the tile QR
 *          or the tile LQ factorization.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_zgeqrf_Tile(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_ZGELS, PlasmaComplexDouble, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_zgelqf - Allocates workspace for PLASMA_zgelqf or PLASMA_zgelqf_Tile routines.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile LQ
 *          factorization.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_zgelqf(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_ZGELS, PlasmaComplexDouble, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_zgelqf_Tile - Allocates tile workspace for PLASMA_zgels_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, PLASMA_desc on workspace handle for storage of the extra T factors required by the tile QR
 *          or the tile LQ factorization.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_zgelqf_Tile(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_ZGELS, PlasmaComplexDouble, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_zgesv - Allocates workspace for PLASMA_zgesv or PLASMA_zgesv_Tile routines.
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The number of linear equations, i.e., the order of the matrix A. N >= 0.
 *
 * @param[out] descL
 *          On exit, workspace handle for storage of the extra L factors required by the tile LU
 *          factorization.
 *
 * @param[out] IPIV
 *          On exit, workspace handle for storage of pivot indexes required by the tile LU
 *          factorization (not equivalent to LAPACK).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_zgesv_incpiv(int N, PLASMA_desc **descL, int **IPIV) {
    int status = plasma_alloc_ibnb_tile(N, N, PLASMA_FUNC_ZGESV, PlasmaComplexDouble, descL);
    if (status != PLASMA_SUCCESS)
        return status;
    return plasma_alloc_ipiv(N, N, PLASMA_FUNC_ZGESV, (void**)IPIV); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_zgesv_Tile - Allocates workspace for PLASMA_zgesv_Tile routines.
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The number of linear equations, i.e., the order of the matrix A. N >= 0.
 *
 * @param[out] descL
 *          On exit, PLASMA descriptor on workspace handle for storage of the extra
 *          L factors required by the tile LU factorization.
 *
 * @param[out] IPIV
 *          On exit, workspace handle for storage of pivot indexes required by the tile LU
 *          factorization (not equivalent to LAPACK).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_zgesv_incpiv_Tile(int N, PLASMA_desc **descL, int **IPIV)
{
    int status = plasma_alloc_ibnb_tile(N, N, PLASMA_FUNC_ZGESV, PlasmaComplexDouble, descL);
    if (status != PLASMA_SUCCESS)
        return status;
    return plasma_alloc_ipiv(N, N, PLASMA_FUNC_ZGESV, (void **)IPIV);
}
/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_zgesvd - Allocates workspace for PLASMA_zgesvd or PLASMA_zgesvd_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile BRD.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_zgesvd(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_ZGESVD, PlasmaComplexDouble, descT); }

/***************************************************************************//**
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_zgetrf_incpiv - Allocates workspace for
 *  PLASMA_zgetrf_incpiv or PLASMA_zgetrf_incpiv_Tile or
 *  PLASMA_zgetrf_incpiv_Tile_Async routines.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descL
 *          On exit, workspace handle for storage of the extra L factors required by the tile LU
 *          factorization.
 *
 * @param[out] IPIV
 *          On exit, workspace handle for storage of pivot indexes required by the tile LU
 *          factorization (not equivalent to LAPACK).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************
 *
 * @sa PLASMA_zgetrf_incpiv
 * @sa PLASMA_zgetrf_incpiv_Tile
 * @sa PLASMA_zgetrf_incpiv_Tile_Async
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_zgetrf_incpiv(int M, int N, PLASMA_desc **descL, int **IPIV) {
    int status = plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_ZGESV, PlasmaComplexDouble, descL);
    if (status != PLASMA_SUCCESS)
        return status;
    return plasma_alloc_ipiv(M, N, PLASMA_FUNC_ZGESV, (void**)IPIV); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_zgetrf_incpiv_Tile - Allocates workspace for
 *  PLASMA_zgesv_incpiv_Tile or PLASMA_zgesv_incpiv_Tile_Async routines.
 *
 *******************************************************************************
 *
 * @param[in] N
 *          The number of linear equations, i.e., the order of the matrix A. N >= 0.
 *
 * @param[out] descL
 *          On exit, PLASMA descriptor on workspace handle for storage of the extra
 *          L factors required by the tile LU factorization.
 *
 * @param[out] IPIV
 *          On exit, workspace handle for storage of pivot indexes required by the tile LU
 *          factorization (not equivalent to LAPACK).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_zgetrf_incpiv_Tile(int N, PLASMA_desc **descL, int **IPIV)
{
    int status = plasma_alloc_ibnb_tile(N, N, PLASMA_FUNC_ZGESV, PlasmaComplexDouble, descL);
    if (status != PLASMA_SUCCESS)
        return status;
    return plasma_alloc_ipiv(N, N, PLASMA_FUNC_ZGESV, (void **)IPIV);
}
/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_zheev - Allocates workspace for PLASMA_zheev or PLASMA_zheev_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile TRD.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_zheev(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_ZHEEV, PlasmaComplexDouble, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_zheevd - Allocates workspace for PLASMA_zheevd or PLASMA_zheevd_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile TRD.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_zheevd(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_ZHEEVD, PlasmaComplexDouble, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_zhegv - Allocates workspace for PLASMA_zhegv or PLASMA_zhegv_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile TRD.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_zhegv(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_ZHEGV, PlasmaComplexDouble, descT); }

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_zhegvd - Allocates workspace for PLASMA_zhegvd or PLASMA_zhegvd_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile TRD.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_zhegvd(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_ZHEGVD, PlasmaComplexDouble, descT); }

 /***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Alloc_Workspace_zhetrd - Allocates workspace for PLASMA_zhetrd or PLASMA_zhetrd_Tile routine.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A.  N >= 0.
 *
 * @param[out] descT
 *          On exit, workspace handle for storage of the extra T factors required by the tile TRD.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Alloc_Workspace_zhetrd(int M, int N, PLASMA_desc **descT) {
    return plasma_alloc_ibnb_tile(M, N, PLASMA_FUNC_ZHETRD, PlasmaComplexDouble, descT); }
