/**
 *
 * @file control/auxiliary.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Jakub Kurzak
 * @author Piotr Luszczek
 * @author Emmanuel Agullo
 * @date 2010-11-15
 *
 **/
#include "common.h"
#include "auxiliary.h"

#include <stdio.h>
#include <stdlib.h>

/***************************************************************************//**
 *
 *  Indicates a recoverable problem.
 *  User's erroneous action without severe consequences.
 *  Problems occuring while PLASMA is being used correctly.
 *  Context aware.
 *
 * @param[in] func_name
 *          Function location where warning occurred
 *
 * @param[in] msg_text
 *          Warning message to display.
 *
 ******************************************************************************/
void plasma_warning(const char *func_name, char* msg_text)
{
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL)
        plasma_fatal_error("plasma_warning", "PLASMA not initialized");
    if (plasma->warnings_enabled)
        fprintf(stderr, "PLASMA WARNING: %s(): %s\n", func_name, msg_text);
}

/***************************************************************************//**
 *
 *  Indicates a recoverable problem.
 *  User's erroneous action with potentially severe consequences.
 *  Problems occuring due to incorrect use of PLASMA.
 *  Context aware.
 *
 * @param[in] func_name
 *          Function location where warning occurred
 *
 * @param[in] msg_text
 *          Warning message to display.
 *
 ******************************************************************************/
void plasma_error(const char *func_name, char* msg_text)
{
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL)
        plasma_fatal_error("plasma_error", "PLASMA not initialized");
    if (plasma->errors_enabled)
        fprintf(stderr, "PLASMA ERROR: %s(): %s\n", func_name, msg_text);

}

/***************************************************************************//**
 *
 *  Unexpected behavior within the library.
 *  Unrecoverable user errors.
 *  Context oblivious.
 *
 * @param[in] func_name
 *          Function location where warning occurred
 *
 * @param[in] msg_text
 *          Warning message to display.
 *
 ******************************************************************************/
void plasma_fatal_error(const char *func_name, char* msg_text)
{
    fprintf(stderr, "PLASMA FATAL ERROR: %s(): %s\n", func_name, msg_text);
    exit(0);
}

/***************************************************************************//**
 *
 **/
int plasma_element_size(int type)
{
    switch(type) {
        case PlasmaByte:          return          1;
        case PlasmaInteger:       return   sizeof(int);
        case PlasmaRealFloat:     return   sizeof(float);
        case PlasmaRealDouble:    return   sizeof(double);
        case PlasmaComplexFloat:  return 2*sizeof(float);
        case PlasmaComplexDouble: return 2*sizeof(double);
        default: plasma_fatal_error("plasma_element_size", "undefined type");
                 return PLASMA_ERR_ILLEGAL_VALUE;

    }
}

/***************************************************************************//**
 *
 **/
void plasma_memcpy(void *dst, void *src, PLASMA_size size, int type)
{
    memcpy(dst, src, size*plasma_element_size(type));
}

/***************************************************************************//**
 *
 **/
void plasma_memzero(void *memptr, PLASMA_size size, int type)
{
    memset(memptr, 0, size*plasma_element_size(type));
}

/***************************************************************************//**
 *
 **/
void plasma_memset_int(int *mem, int size, int value)
{
    int i;

    for (i = 0; i < size; i++)
        mem[i] = value;
}

/***************************************************************************//**
 *  Returns core id
 **/
int plasma_rank(plasma_context_t *plasma)
{
    int rank;
    pthread_t thread_id;

    thread_id = pthread_self();
    for (rank = 0; rank < plasma->world_size; rank++)
        if (pthread_equal(plasma->thread_id[rank], thread_id))
            return rank;
    return PLASMA_ERR_NOT_FOUND;
}

/***************************************************************************//**
 *  Tune block size nb and internal block size ib
 **/
int plasma_tune(PLASMA_enum func, int M, int N, int NRHS)
{
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("plasma_tune", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }

    if (!plasma->autotuning_enabled)
        return PLASMA_SUCCESS;

    switch (func) {
        case PLASMA_FUNC_SGEMM:
        case PLASMA_FUNC_DGEMM:
        case PLASMA_FUNC_CGEMM:
        case PLASMA_FUNC_ZGEMM:
        case PLASMA_FUNC_CHERK:
        case PLASMA_FUNC_ZHERK:
        case PLASMA_FUNC_SSYRK:
        case PLASMA_FUNC_DSYRK:
        case PLASMA_FUNC_CSYRK:
        case PLASMA_FUNC_ZSYRK:
        case PLASMA_FUNC_CHEMM:
        case PLASMA_FUNC_ZHEMM:
        case PLASMA_FUNC_SSYMM:
        case PLASMA_FUNC_DSYMM:
        case PLASMA_FUNC_CSYMM:
        case PLASMA_FUNC_ZSYMM:
            plasma->nb = 120;
            plasma->ib = 120;  // not used in GEMM
            break;
        case PLASMA_FUNC_SPOSV:
        case PLASMA_FUNC_DPOSV:
        case PLASMA_FUNC_CPOSV:
        case PLASMA_FUNC_ZPOSV:
        case PLASMA_FUNC_ZCPOSV:
        case PLASMA_FUNC_DSPOSV:
            plasma->nb = 120;
            plasma->ib = 120;  // not used in Cholesky
            break;
        case PLASMA_FUNC_SGELS:
        case PLASMA_FUNC_DGELS:
        case PLASMA_FUNC_CGELS:
        case PLASMA_FUNC_ZGELS:
        case PLASMA_FUNC_ZCGELS:
        case PLASMA_FUNC_DSGELS:
            plasma->nb = 144;
            plasma->ib =  48;
            break;
        case PLASMA_FUNC_SGESV:
        case PLASMA_FUNC_DGESV:
        case PLASMA_FUNC_CGESV:
        case PLASMA_FUNC_ZGESV:
        case PLASMA_FUNC_ZCGESV:
        case PLASMA_FUNC_DSGESV:
            plasma->nb = 200;
            plasma->ib =  40;
            break;
        case PLASMA_FUNC_SGEEV:
        case PLASMA_FUNC_DGEEV:
        case PLASMA_FUNC_CGEEV:
        case PLASMA_FUNC_ZGEEV:
        case PLASMA_FUNC_SGEHRD:
        case PLASMA_FUNC_DGEHRD:
        case PLASMA_FUNC_CGEHRD:
        case PLASMA_FUNC_ZGEHRD:
            plasma->nb = 120;
            plasma->ib =  20;
            break;
        case PLASMA_FUNC_ZHEEV:
        case PLASMA_FUNC_CHEEV:
        case PLASMA_FUNC_DSYEV:
        case PLASMA_FUNC_SSYEV:
        case PLASMA_FUNC_ZHEEVD:
        case PLASMA_FUNC_CHEEVD:
        case PLASMA_FUNC_DSYEVD:
        case PLASMA_FUNC_SSYEVD:
        case PLASMA_FUNC_ZHETRD:
        case PLASMA_FUNC_CHETRD:
        case PLASMA_FUNC_DSYTRD:
        case PLASMA_FUNC_SSYTRD:
        case PLASMA_FUNC_ZHEGV:
        case PLASMA_FUNC_CHEGV:
        case PLASMA_FUNC_DSYGV:
        case PLASMA_FUNC_SSYGV:
        case PLASMA_FUNC_ZHEGVD:
        case PLASMA_FUNC_CHEGVD:
        case PLASMA_FUNC_DSYGVD:
        case PLASMA_FUNC_SSYGVD:
        case PLASMA_FUNC_ZHEGST:
        case PLASMA_FUNC_CHEGST:
        case PLASMA_FUNC_DSYGST:
        case PLASMA_FUNC_SSYGST:
            plasma->nb = 120;
            plasma->ib =  20;
            break;
        case PLASMA_FUNC_ZGESVD:
        case PLASMA_FUNC_CGESVD:
        case PLASMA_FUNC_DGESVD:
        case PLASMA_FUNC_SGESVD:
        case PLASMA_FUNC_ZGEBRD:
        case PLASMA_FUNC_CGEBRD:
        case PLASMA_FUNC_DGEBRD:
        case PLASMA_FUNC_SGEBRD:
            plasma->nb = 120;
            plasma->ib =  20;
            break;
        default:
            plasma_fatal_error("plasma_tune", "illegal parameter value");
            return PLASMA_ERR_ILLEGAL_VALUE;
    }
    /* Calculate A, B tile size and round up to cache line size */
    /* round up for the smallest type (float) - will hold for all */
    plasma->nbnbsize = plasma->nb * plasma->nb; // * sizeof(float);
//  plasma->nbnbsize = roundup(plasma->nbnbsize, CACHE_LINE_SIZE);
//  plasma->nbnbsize /= sizeof(float);
    /* Calculate T, L tile size and round up to cache line size */
    /* round up for the smallest type (float) - will hold for all */
    plasma->ibnbsize = plasma->ib * plasma->nb; // * sizeof(float);
//  plasma->ibnbsize = roundup(plasma->ibnbsize, CACHE_LINE_SIZE);
//  plasma->ibnbsize /= sizeof(float);
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 * @ingroup Auxiliary
 *
 *  PLASMA_Version - Reports PLASMA version number.
 *
 *******************************************************************************
 *
 * @param[out] ver_major
 *          PLASMA major version number.
 *
 * @param[out] ver_minor
 *          PLASMA minor version number.
 *
 * @param[out] ver_micro
 *          PLASMA micro version number.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 ******************************************************************************/
int PLASMA_Version(int *ver_major, int *ver_minor, int *ver_micro)
{
    if (! ver_major && ! ver_minor && ! ver_micro)
        return  PLASMA_ERR_ILLEGAL_VALUE;

    if (ver_major)
        *ver_major = PLASMA_VERSION_MAJOR;

    if (ver_minor)
        *ver_minor = PLASMA_VERSION_MINOR;

    if (ver_micro)
        *ver_micro = PLASMA_VERSION_MICRO;

    return PLASMA_SUCCESS;
}
