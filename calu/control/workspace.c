/**
 *
 * @file workspace.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Jakub Kurzak
 * @date 2010-11-15
 *
 **/
#include "common.h"
#include "allocate.h"
#include "auxiliary.h"
#include "workspace.h"

#include <stdlib.h>

/***************************************************************************//**
 *
 **/
int plasma_alloc_ibnb(int M, int N, PLASMA_enum func, int type, void **memptr)
{
    size_t size;
    int status;
    int IB, NB, MT, NT;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("plasma_alloc_ibnb", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Tune NB & IB depending on M & N; Set IBNBSIZE */
    status = plasma_tune(func, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("plasma_alloc_ibnb", "plasma_tune() failed");
        return PLASMA_ERR_UNEXPECTED;
    }
    /* Set MT & NT & allocate */
    NB = PLASMA_NB;
    IB = PLASMA_IB;
    MT = (M%NB==0) ? (M/NB) : (M/NB+1);
    NT = (N%NB==0) ? (N/NB) : (N/NB+1);

    /* Size is doubled for RH QR to store the reduction T */
    if ((plasma->householder != PLASMA_FLAT_HOUSEHOLDER) &&
        (func == PLASMA_FUNC_SGELS  ||
         func == PLASMA_FUNC_DGELS  ||
         func == PLASMA_FUNC_CGELS  ||
         func == PLASMA_FUNC_ZGELS  ||
         func == PLASMA_FUNC_SGESVD ||
         func == PLASMA_FUNC_DGESVD ||
         func == PLASMA_FUNC_CGESVD ||
         func == PLASMA_FUNC_ZGESVD ))
        NT *= 2;

    size = (size_t)MT*NT*IB*NB * plasma_element_size(type);
    if (size <= 0) {
        *memptr = NULL;
        return PLASMA_SUCCESS;
    }
//  status = posix_memalign(memptr, STANDARD_PAGE_SIZE, size);
    *memptr = malloc(size);
//  if (status != 0) {
    if (*memptr == NULL) {
        plasma_error("plasma_alloc_ibnb_tile", "malloc() failed");
        return PLASMA_ERR_OUT_OF_RESOURCES;
    }
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 **/
int plasma_alloc_ibnb_tile(int M, int N, PLASMA_enum func, int type, PLASMA_desc **desc)
{
    int status;
    int IB, NB, MT, NT;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("plasma_alloc_ibnb_tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Tune NB & IB depending on M & N; Set IBNBSIZE */
    status = plasma_tune(func, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("plasma_alloc_ibnb_tile", "plasma_tune() failed");
        return PLASMA_ERR_UNEXPECTED;
    }

    /* Set MT & NT & allocate */
    NB = PLASMA_NB;
    IB = PLASMA_IB;
    MT = (M%NB==0) ? (M/NB) : (M/NB+1);
    NT = (N%NB==0) ? (N/NB) : (N/NB+1);

    /* Size is doubled for RH QR to store the reduction T */
    if ((plasma->householder != PLASMA_FLAT_HOUSEHOLDER) &&
        ((func == PLASMA_FUNC_SGELS)  ||
         (func == PLASMA_FUNC_DGELS)  ||
         (func == PLASMA_FUNC_CGELS)  ||
         (func == PLASMA_FUNC_ZGELS)  ||
         (func == PLASMA_FUNC_SGESVD) ||
         (func == PLASMA_FUNC_DGESVD) ||
         (func == PLASMA_FUNC_CGESVD) ||
         (func == PLASMA_FUNC_ZGESVD)))
        NT *= 2;

    /* Allocate and initialize descriptor */
    *desc = (PLASMA_desc*)malloc(sizeof(PLASMA_desc));
    if (*desc == NULL) {
        plasma_error("plasma_alloc_ibnb_tile", "malloc() failed");
        return PLASMA_ERR_OUT_OF_RESOURCES;
    }
    **desc = plasma_desc_init(type, IB, NB, IB*NB, MT*IB, NT*NB, 0, 0, MT*IB, NT*NB);

    /* Allocate matrix */
    if (plasma_desc_mat_alloc(*desc)) {
        plasma_error("plasma_alloc_ibnb_tile", "malloc() failed");
        return PLASMA_ERR_OUT_OF_RESOURCES;
    }

    /* Check that everything is ok */
    status = plasma_desc_check(*desc);
    if (status != PLASMA_SUCCESS) {
        plasma_error("plasma_alloc_ibnb_tile", "invalid descriptor");
        return status;
    }
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 **/
int plasma_alloc_ipiv(int M, int N, PLASMA_enum func, void **memptr)
{
    size_t size;
    int status;
    int NB, MT, NT;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("plasma_alloc_ipiv", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Tune NB & IB depending on M & N; Set IBNBSIZE */
    status = plasma_tune(func, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("plasma_alloc_ipiv", "plasma_tune() failed");
        return PLASMA_ERR_UNEXPECTED;
    }
    /* Set MT & NT & allocate */
    NB = PLASMA_NB;
    NT = (N%NB==0) ? (N/NB) : ((N/NB)+1);
    MT = (M%NB==0) ? (M/NB) : ((M/NB)+1);
    size = (size_t)MT*NT * NB * sizeof(int);
    if (size <= 0) {
        *memptr = NULL;
        return PLASMA_SUCCESS;
    }
//  status = posix_memalign(memptr, CACHE_LINE_SIZE, size);
    *memptr = malloc(size);
//  if (status != 0) {
    if (*memptr == NULL) {
        plasma_error("plasma_alloc_ipiv", "malloc() failed");
        return PLASMA_ERR_OUT_OF_RESOURCES;
    }
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Dealloc_Handle - Deallocate workspace handle allocated by any workspace allocation routine.
 *
 *******************************************************************************
 *
 * @param[in] handle
 *          Workspace handle
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Dealloc_Handle(void **handle)
{
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_Dealloc_Handle", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (*handle == NULL) {
        plasma_error("PLASMA_Dealloc_Handle", "attempting to deallocate a NULL handle");
        return PLASMA_ERR_UNALLOCATED;
    }
    free(*handle);
    *handle = NULL;
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Dealloc_Handle_Tile - Deallocate Tile workspace handle allocated by any tile workspace allocation routine.
 *
 *******************************************************************************
 *
 * @param[in] desc
 *          Descriptot handle
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Dealloc_Handle_Tile(PLASMA_desc **desc)
{
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_Dealloc_Handle_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (*desc == NULL) {
        plasma_error("PLASMA_Dealloc_Handle_Tile", "attempting to deallocate a NULL descriptor");
        return PLASMA_ERR_UNALLOCATED;
    }
    if ((*desc)->mat == NULL) {
        plasma_error("PLASMA_Dealloc_Handle_Tile", "attempting to deallocate a NULL pointer");
        return PLASMA_ERR_UNALLOCATED;
    }
    free((*desc)->mat);
    free(*desc);
    *desc = NULL;
    return PLASMA_SUCCESS;
}
