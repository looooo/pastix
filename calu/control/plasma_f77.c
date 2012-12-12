/**
 *
 * @file plasma_f77.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Bilel Hadri
 * @date 2010-11-15
 *
 **/
#include <stdlib.h>
#include "common.h"

#ifdef ADD_
    #define PLASMA_INIT             plasma_init_
    #define PLASMA_FINALIZE         plasma_finalize_
    #define PLASMA_ENABLE           plasma_enable_
    #define PLASMA_DISABLE          plasma_disable_
    #define PLASMA_SET              plasma_set_
    #define PLASMA_GET              plasma_get_
    #define PLASMA_DEALLOC_HANDLE   plasma_dealloc_handle_
    #define PLASMA_VERSION          plasma_version_
    #define PLASMA_DESC_CREATE      plasma_desc_create_
    #define PLASMA_DESC_DESTROY     plasma_desc_destroy_
    #define PLASMA_LAPACK_TO_TILE   plasma_lapack_to_tile_
    #define PLASMA_TILE_TO_LAPACK   plasma_tile_to_lapack_
#elif defined (NOCHANGE)
    #define PLASMA_INIT             plasma_init
    #define PLASMA_FINALIZE         plasma_finalize
    #define PLASMA_ENABLE           plasma_enable
    #define PLASMA_DISABLE          plasma_disable
    #define PLASMA_SET              plasma_set
    #define PLASMA_GET              plasma_get
    #define PLASMA_DEALLOC_HANDLE   plasma_dealloc_handle
    #define PLASMA_VERSION          plasma_version
    #define PLASMA_DESC_CREATE      plasma_desc_create
    #define PLASMA_DESC_DESTROY     plasma_desc_destroy
 #define PLASMA_LAPACK_TO_TILE   plasma_lapack_to_tile
 #define PLASMA_TILE_TO_LAPACK   plasma_tile_to_lapack
#endif

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  FORTRAN API - auxiliary function prototypes
 **/
void PLASMA_INIT(int *CORES, int *INFO)
{   *INFO = PLASMA_Init(*CORES); }

void PLASMA_FINALIZE(int *INFO)
{   *INFO = PLASMA_Finalize(); }

void PLASMA_ENABLE(PLASMA_enum *lever, int *INFO)
{   *INFO = PLASMA_Enable(*lever); }

void PLASMA_DISABLE(PLASMA_enum *lever, int *INFO)
{   *INFO = PLASMA_Disable(*lever); }

void PLASMA_SET(PLASMA_enum *param, int *value, int *INFO)
{   *INFO = PLASMA_Set(*param, *value); }

void PLASMA_GET(PLASMA_enum *param, int *value, int *INFO)
{   *INFO = PLASMA_Get(*param, value); }

void PLASMA_DEALLOC_HANDLE(size_t *sp, int *INFO)
{   free((void *)(*sp));
    *INFO = PLASMA_SUCCESS; }

void PLASMA_VERSION(int *VER_MAJOR, int *VER_MINOR, int *VER_MICRO, int *INFO)
{
    *VER_MAJOR = PLASMA_VERSION_MAJOR;
    *VER_MINOR = PLASMA_VERSION_MINOR;
    *VER_MICRO = PLASMA_VERSION_MICRO;
    *INFO = PLASMA_SUCCESS;
}

/***************************************************************************//**
 *  FORTRAN API - descriptor allocation and deallocation
 **/
void PLASMA_DESC_CREATE(PLASMA_desc **desc, void *mat, PLASMA_enum *dtyp, int *mb, int *nb, int *bsiz, int *lm, int *ln, int *i, int *j, int *m, int *n, int *INFO)
{   *INFO = PLASMA_Desc_Create(desc, mat, *dtyp, *mb, *nb, *bsiz, *lm, *ln, *i, *j, *m, *n); }

void PLASMA_DESC_DESTROY(PLASMA_desc **desc, int *INFO)
{   *INFO = PLASMA_Desc_Destroy(desc); }

/***************************************************************************//**
 *  FORTRAN API - conversion from LAPACK F77 matrix layout to tile layout
 **/
void PLASMA_LAPACK_TO_TILE(intptr_t *Af77, int *LDA, intptr_t *A, int *INFO)
{   *INFO = PLASMA_Lapack_to_Tile( (void *)Af77, *LDA, (PLASMA_desc *)(*A)); }

void PLASMA_TILE_TO_LAPACK(intptr_t *A, intptr_t *Af77, int *LDA, int *INFO)
{   *INFO = PLASMA_Tile_to_Lapack((PLASMA_desc *)(*A), (void *)Af77, *LDA); }

#ifdef __cplusplus
}
#endif
