/**
 *
 * @file pastix_starpu_interface.c
 *
 * Interface used by StarPU to handle factorization.
 *
 * @copyright 2021-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.2
 * @author Nolan Bredel
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @author Florent Pruvost
 * @author Tom Moenne-Loccoz
 * @date 2023-12-01
 *
 **/
#include "common/common.h"
#include "include/pastix/datatypes.h"
#include "blend/solver.h"
#include "kernels/pastix_zlrcores.h"
#include "kernels/pastix_clrcores.h"
#include "kernels/pastix_dlrcores.h"
#include "kernels/pastix_slrcores.h"
#include "pastix_starpu.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#if defined(PASTIX_STARPU_INTERFACE_DEBUG)
#define pastix_starpu_logger fprintf( stderr, "pastix_starpu: %s\n", __func__ )
#else
#define pastix_starpu_logger do {} while(0)
#endif

static inline void
psi_register_data_handle( starpu_data_handle_t handle, int home_node, void *data_interface )
{
    pastix_starpu_interface_t *interf = (pastix_starpu_interface_t *)data_interface;
    int                        node;

    pastix_starpu_logger;

    for ( node = 0; node < STARPU_MAXNODES; node++ ) {
        pastix_starpu_interface_t *local_interface =
            (pastix_starpu_interface_t *)starpu_data_get_interface_on_node( handle, node );

        memcpy( local_interface, interf, sizeof( pastix_starpu_interface_t ) );

        if ( node != home_node ) {
            local_interface->dataptr = NULL;
        }
    }
}

static inline starpu_ssize_t
psi_allocate_data_on_node( void *data_interface, unsigned node )
{
    pastix_starpu_interface_t *interf = (pastix_starpu_interface_t *)data_interface;
    starpu_ssize_t             allocated_memory;
    uintptr_t                  addr = 0;
    uintptr_t                  handle;

    pastix_starpu_logger;

    allocated_memory = interf->allocsize;
    if ( allocated_memory <= 0 ) {
        return 0;
    }

    handle = starpu_malloc_on_node( node, allocated_memory );
    if ( !handle ) {
        return -ENOMEM;
    }

#if defined(PASTIX_DEBUG_STARPU)
    {
        const SolverCblk *cblk = interf->cblk;
        if ( cblk->cblktype & CBLK_FANIN ) {
            fprintf( stderr, "allocate fanin %d [%p](%ld)\n",
                     cblk->gfaninnum, (void*)handle, allocated_memory );
        }
        if ( cblk->cblktype & CBLK_RECV ) {
            fprintf( stderr, "allocate recv %d [%p](%ld)\n",
                     cblk->gfaninnum, (void*)handle, allocated_memory );
        }
    }
#endif

    if ( starpu_node_get_kind( node ) != STARPU_OPENCL_RAM ) {
        addr = handle;
    }

    /* update the data properly */
    interf->dataptr = (void *)addr;

    /* /\* Allocate the workspace for the low-rank blocks *\/ */
    /* if ( interf->cblk->cblktype & CBLK_COMPRESSED ) */
    /* { */
    /*     SolverBlok       *blok    = interf->cblk->fblokptr; */
    /*     pastix_lrblock_t *LRblock = interf->dataptr; */
    /*     int               offset  = pastix_imax( 0, interf->offset ); */
    /*     int               i, ncols, M; */

    /*     assert( node == STARPU_MAIN_RAM ); */

    /*     ncols = cblk_colnbr( interf->cblk ); */
    /*     blok += offset; */
    /*     for ( i = 0; i < interf->nbblok; i++, blok++, LRblock++ ) { */
    /*         M = blok_rownbr( blok ); */

    /*         /\* Allocate the LR block to its max space *\/ */
    /*         switch ( interf->flttype ) { */
    /*             case PastixComplex64: */
    /*                 core_zlralloc( M, ncols, -1, LRblock ); */
    /*                 break; */
    /*             case PastixComplex32: */
    /*                 core_clralloc( M, ncols, -1, LRblock ); */
    /*                 break; */
    /*             case PastixDouble: */
    /*                 core_dlralloc( M, ncols, -1, LRblock ); */
    /*                 break; */
    /*             case PastixFloat: */
    /*                 core_slralloc( M, ncols, -1, LRblock ); */
    /*                 break; */
    /*             default: */
    /*                 assert( 0 ); */
    /*         } */
    /*     } */
    /* } */

    return allocated_memory;
}

static inline void
psi_free_data_on_node( void *data_interface, unsigned node )
{
    pastix_starpu_interface_t *interf = (pastix_starpu_interface_t *)data_interface;

    /* SolverCblk *cblk = interf->cblk; */

    pastix_starpu_logger;

    starpu_free_on_node( node, (uintptr_t)interf->dataptr, interf->allocsize );

    interf->dataptr = NULL;
}

static inline void
psi_init( void *data_interface )
{
    pastix_starpu_interface_t *interf = data_interface;
    interf->id                        = PASTIX_STARPU_INTERFACE_ID;
    interf->allocsize                 = -1;

    pastix_starpu_logger;
}

static inline void *
psi_to_pointer( void *data_interface, unsigned node )
{
    (void)node;
    pastix_starpu_interface_t *interf = (pastix_starpu_interface_t *)data_interface;

    pastix_starpu_logger;

    return interf->dataptr;
}

static inline size_t
psi_get_size( starpu_data_handle_t handle )
{
    pastix_starpu_interface_t *interf =
        starpu_data_get_interface_on_node( handle, STARPU_MAIN_RAM );
    const SolverCblk *cblk  = interf->cblk;
    pastix_int_t      ncols = cblk_colnbr( cblk );
    size_t            nrows;
    size_t            size;

    if ( interf->offset == -1 ) {
        nrows = cblk->stride;
    }
    else {
        SolverBlok *fblok = interf->cblk->fblokptr + interf->offset;
        SolverBlok *lblok = fblok + interf->nbblok;

        nrows = 0;
        for( ; fblok < lblok; fblok++ ) {
            nrows += blok_rownbr( fblok );
        }
    }

    size = ncols * nrows;

#ifdef STARPU_DEBUG
    STARPU_ASSERT_MSG( interf->id == PASTIX_STARPU_INTERFACE_ID,
                       "psi_get_size: The given data is not a pastix interface for starpu." );
#endif

    return size;
}

static inline size_t
psi_get_alloc_size( starpu_data_handle_t handle )
{
    pastix_starpu_interface_t *interf =
        starpu_data_get_interface_on_node( handle, STARPU_MAIN_RAM );

    pastix_starpu_logger;

#ifdef STARPU_DEBUG
    STARPU_ASSERT_MSG( interf->id == PASTIX_STARPU_INTERFACE_ID,
                       "psi_get_alloc_size: The given data is not a pastix interface for starpu." );
#endif

    STARPU_ASSERT_MSG( interf->allocsize != (size_t)-1,
                       "psi_get_alloc_size: The allocation size needs to be defined" );

    return interf->allocsize;
}

static inline uint32_t
psi_footprint( starpu_data_handle_t handle )
{
    pastix_starpu_interface_t *interf =
        starpu_data_get_interface_on_node( handle, STARPU_MAIN_RAM );
    const SolverCblk *cblk = interf->cblk;

    pastix_starpu_logger;

    return starpu_hash_crc32c_be( cblk->gcblknum, interf->offset + 1 );
}

static inline uint32_t
psi_alloc_footprint( starpu_data_handle_t handle )
{
    pastix_starpu_interface_t *interf =
        starpu_data_get_interface_on_node( handle, STARPU_MAIN_RAM );

    pastix_starpu_logger;

    return starpu_hash_crc32c_be( interf->allocsize, 0 );
}

static inline int
psi_compare( void *data_interface_a, void *data_interface_b )
{
    pastix_starpu_interface_t *pastix_interface_a = (pastix_starpu_interface_t *)data_interface_a;
    pastix_starpu_interface_t *pastix_interface_b = (pastix_starpu_interface_t *)data_interface_b;

    const SolverCblk *solva = pastix_interface_a->cblk;
    const SolverCblk *solvb = pastix_interface_b->cblk;

    pastix_starpu_logger;

    /* Two interfaces are considered compatible if they point to the same solver */
    return ( solva == solvb );
}

static inline int
psi_alloc_compare( void *data_interface_a, void *data_interface_b )
{
    pastix_starpu_interface_t *pastix_interface_a = (pastix_starpu_interface_t *)data_interface_a;
    pastix_starpu_interface_t *pastix_interface_b = (pastix_starpu_interface_t *)data_interface_b;

    pastix_starpu_logger;

    /* Two matrices are considered compatible if they have the same allocated size */
    return ( pastix_interface_a->allocsize == pastix_interface_b->allocsize );
}

static inline void
psi_display( starpu_data_handle_t handle, FILE *f )
{
    pastix_starpu_interface_t *interf =
        (pastix_starpu_interface_t *)starpu_data_get_interface_on_node( handle, STARPU_MAIN_RAM );

    const SolverCblk *cblk = interf->cblk;

    pastix_starpu_logger;

    if ( interf->offset == -1 ) {
        fprintf( f, "Cblk%ld", (long)( cblk->gcblknum ) );
    }
    else {
        fprintf( f, "Cblk%ldBlok%ld", (long)( cblk->gcblknum ), (long)( interf->offset ) );
    }
}

static inline size_t
psi_compute_size_lr( pastix_starpu_interface_t *interf )
{
    assert( interf->cblk->cblktype & CBLK_COMPRESSED );

    size_t            elemsize = pastix_size_of( interf->flttype );
    pastix_lrblock_t *LRblock  = interf->dataptr;
    pastix_int_t      N        = cblk_colnbr( interf->cblk );
    pastix_int_t      suv, M;

    SolverBlok *blok = interf->cblk->fblokptr + pastix_imax( 0, interf->offset );

    pastix_starpu_logger;

    suv = 0;
    int    i;
    for ( i = 0; i < interf->nbblok; i++, blok++, LRblock++ ) {
        M = blok_rownbr( blok );
        suv += sizeof( int ) + core_zlrgetsize( M, N, LRblock ) * elemsize;
    }
    return suv;
}

static inline size_t
psi_compute_size_fr( pastix_starpu_interface_t *interf )
{
    assert( !( interf->cblk->cblktype & CBLK_COMPRESSED ) );

    pastix_starpu_logger;

    return interf->allocsize;
}

static inline void
psi_pack_lr( pastix_starpu_interface_t *interf, void **ptr )
{
    assert( interf->cblk->cblktype & CBLK_COMPRESSED );

    char             *tmp      = *ptr;
    pastix_lrblock_t *LRblock  = interf->dataptr;
    int               N        = cblk_colnbr( interf->cblk );
    int               j        = 0;
    int               M;

    SolverBlok *blok = interf->cblk->fblokptr + pastix_imax( 0, interf->offset );

    pastix_starpu_logger;

    for ( ; j < interf->nbblok; j++, blok++, LRblock++ ) {
        M = blok_rownbr( blok );

        switch ( interf->flttype ) {
        case PastixComplex64:
            tmp = core_zlrpack( M, N, LRblock, tmp );
            break;
        case PastixComplex32:
            tmp = core_clrpack( M, N, LRblock, tmp );
            break;
        case PastixDouble:
            tmp = core_dlrpack( M, N, LRblock, tmp );
            break;
        case PastixFloat:
            tmp = core_slrpack( M, N, LRblock, tmp );
            break;
        default:
            assert( 0 );
        }
    }
}

static inline void
psi_pack_fr( pastix_starpu_interface_t *interf, void **ptr, starpu_ssize_t *count )
{
    assert( !( interf->cblk->cblktype & CBLK_COMPRESSED ) );

    pastix_starpu_logger;

    memcpy( *ptr, interf->dataptr, *count );
}

static inline int
psi_pack_data( starpu_data_handle_t handle, unsigned node, void **ptr, starpu_ssize_t *count )
{
    STARPU_ASSERT( starpu_data_test_if_allocated_on_node( handle, node ) );

    pastix_starpu_interface_t *interf =
        (pastix_starpu_interface_t *)starpu_data_get_interface_on_node( handle, node );

    const SolverCblk *cblk = interf->cblk;

    pastix_starpu_logger;

    if ( cblk->cblktype & CBLK_COMPRESSED ) {
        *count = psi_compute_size_lr( interf );
    }
    else {
        *count = psi_compute_size_fr( interf );
    }

    if ( ptr != NULL ) {
        *ptr = (void *)starpu_malloc_on_node_flags( node, *count, 0 );

        if ( cblk->cblktype & CBLK_COMPRESSED ) {
            psi_pack_lr( interf, ptr );
        }
        else {
            psi_pack_fr( interf, ptr, count );
        }
    }

    return 0;
}

static inline void
psi_unpack_lr( pastix_starpu_interface_t *interf, unsigned node, const void *ptr, size_t count )
{
    SolverBlok       *blok    = interf->cblk->fblokptr + pastix_imax( 0, interf->offset );
    pastix_lrblock_t *LRblock;
    const char       *input   = ptr;
    char             *output;
    size_t            lrsize  = interf->nbblok * sizeof( pastix_lrblock_t );
    int               N       = cblk_colnbr( interf->cblk );
    int               i;

    pastix_starpu_logger;

    assert( interf->cblk->cblktype & CBLK_COMPRESSED );
    assert( interf->allocsize == 0 );

    /* Remove the size of all the rk */
    count -= interf->nbblok * sizeof( int );

    interf->allocsize = count + lrsize;
    psi_allocate_data_on_node( interf, node );

    LRblock = interf->dataptr;
    output  = interf->dataptr;
    output += lrsize;

    for ( i=0; i < interf->nbblok; i++, blok++, LRblock++ ) {
        int M = blok_rownbr( blok );

        /* Allocate the LR block to its tight space */
        switch ( interf->flttype ) {
            case PastixComplex64:
                input = core_zlrunpack2( M, N, LRblock, input, &output );
                break;
            case PastixComplex32:
                input = core_clrunpack2( M, N, LRblock, input, &output );
                break;
            case PastixDouble:
                input = core_dlrunpack2( M, N, LRblock, input, &output );
                break;
            case PastixFloat:
                input = core_slrunpack2( M, N, LRblock, input, &output );
                break;
            default:
                assert( 0 );
        }
    }
}

static inline void
psi_unpack_fr( pastix_starpu_interface_t *interf, void *ptr, size_t count )
{
    pastix_starpu_logger;

    assert( count == interf->allocsize );
    memcpy( interf->dataptr, ptr, count );
}

static inline int
psi_peek_data( starpu_data_handle_t handle, unsigned node, void *ptr, size_t count )
{
    pastix_starpu_interface_t *interf =
        (pastix_starpu_interface_t *)starpu_data_get_interface_on_node( handle, node );

    STARPU_ASSERT( starpu_data_test_if_allocated_on_node( handle, node ) );

    pastix_starpu_logger;

    if ( interf->cblk->cblktype & CBLK_COMPRESSED ) {
        psi_unpack_lr( interf, node, ptr, count );
    }
    else {
        psi_unpack_fr( interf, ptr, count );
    }

    return 0;
}

static inline int
psi_unpack_data( starpu_data_handle_t handle, unsigned node, void *ptr, size_t count )
{
    pastix_starpu_logger;

    psi_peek_data( handle, node, ptr, count );

    /* Free the received information */
    starpu_free_on_node_flags( node, (uintptr_t)ptr, count, 0 );

    return 0;
}

static inline starpu_ssize_t
psi_describe( void *data_interface, char *buf, size_t size )
{
    pastix_starpu_interface_t *interf = (pastix_starpu_interface_t *)data_interface;
    const SolverCblk          *cblk   = interf->cblk;

    pastix_starpu_logger;

    return snprintf( buf, size, "Cblk%ld", (long)( cblk->gcblknum ) );
}

static inline int
psi_copy_any_to_any( void    *src_interface,
                     unsigned src_node,
                     void    *dst_interface,
                     unsigned dst_node,
                     void    *async_data )
{
    pastix_starpu_interface_t *pastix_src = (pastix_starpu_interface_t *)src_interface;
    pastix_starpu_interface_t *pastix_dst = (pastix_starpu_interface_t *)dst_interface;

    int ret = 0;

    pastix_starpu_logger;

    assert( !( pastix_src->cblk->cblktype & CBLK_COMPRESSED ) );

    if ( starpu_interface_copy( (uintptr_t)pastix_src->dataptr, 0, src_node,
                                (uintptr_t)pastix_dst->dataptr, 0, dst_node,
                                pastix_src->allocsize,
                                async_data ) )
    {
        ret = -EAGAIN;
    }
    starpu_interface_data_copy( src_node, dst_node, pastix_src->allocsize );

    return ret;
}

static const struct starpu_data_copy_methods psi_copy_methods = {
    .any_to_any = psi_copy_any_to_any,
};
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 * @ingroup pastix_starpu
 * @{
 *
 * @brief TODO
 */
struct starpu_data_interface_ops pastix_starpu_interface_ops = {
    .register_data_handle  = psi_register_data_handle,
    .allocate_data_on_node = psi_allocate_data_on_node,
    .free_data_on_node     = psi_free_data_on_node,
    .init                  = psi_init,
    .copy_methods          = &psi_copy_methods,
    .to_pointer            = psi_to_pointer,
    .get_size              = psi_get_size,
    .get_alloc_size        = psi_get_alloc_size,
    .footprint             = psi_footprint,
    .alloc_footprint       = psi_alloc_footprint,
    .compare               = psi_compare,
    .alloc_compare         = psi_alloc_compare,
    .display               = psi_display,
    .describe              = psi_describe,
    .pack_data             = psi_pack_data,
#if defined( HAVE_STARPU_DATA_PEEK )
    .peek_data             = psi_peek_data,
#endif
    .unpack_data           = psi_unpack_data,
    .interfaceid           = STARPU_UNKNOWN_INTERFACE_ID,
    .interface_size        = sizeof( pastix_starpu_interface_t ),
    .dontcache             = 0, /* Should we set it to 1 */
    .name                  = "PASTIX_STARPU_INTERFACE"
};

/**
 *******************************************************************************
 *
 * @brief Initialize the interface ID
 *
 ********************************************************************************/
void
pastix_starpu_interface_init()
{
    if ( pastix_starpu_interface_ops.interfaceid == STARPU_UNKNOWN_INTERFACE_ID ) {
        pastix_starpu_interface_ops.interfaceid = starpu_data_interface_get_next_id();

        /*
         * TODO : Add mpi datatype
         */

        /* #if defined(PASTIX_STARPU_USE_MPI_DATATYPES)                                          */
        /* #if defined(HAVE_STARPU_MPI_INTERFACE_DATATYPE_NODE_REGISTER)                         */
        /* starpu_mpi_interface_datatype_node_register( pastix_starpu_interface_ops.interfaceid, */
        /*                                              psi_allocate_datatype_node,              */
        /*                                              psi_free_datatype );                     */
        /* #else                                                                                 */
        /* starpu_mpi_interface_datatype_register( pastix_starpu_interface_ops.interfaceid,      */
        /*                                         psi_allocate_datatype,                        */
        /*                                         psi_free_datatype );                          */
        /* #endif                                                                                */
        /* #endif                                                                                */
    }
}

/**
 *******************************************************************************
 *
 * @brief Finalize the interface and reset the ID
 *
 ********************************************************************************/
void
pastix_starpu_interface_fini()
{
    if ( pastix_starpu_interface_ops.interfaceid != STARPU_UNKNOWN_INTERFACE_ID ) {
        /* #if defined(PASTIX_STARPU_USE_MPI_DATATYPES)                                         */
        /* starpu_mpi_interface_datatype_unregister( pastix_starpu_interface_ops.interfaceid ); */
        /* #endif                                                                               */
        pastix_starpu_interface_ops.interfaceid = STARPU_UNKNOWN_INTERFACE_ID;
    }
}

/**
 *
 *******************************************************************************
 *
 * @brief Register a cblk at the StarPU level
 *
 *******************************************************************************
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
 *
 ********************************************************************************/
void
pastix_starpu_register( starpu_data_handle_t *handleptr,
                        const SolverCblk     *cblk,
                        pastix_coefside_t     side,
                        pastix_coeftype_t     flttype )
{
    pastix_starpu_interface_t interf = {
        .id        = PASTIX_STARPU_INTERFACE_ID,
        .flttype   = flttype,
        .offset    = -1,
        .nbblok    =  0,
        .allocsize = -1,
        .cblk      = cblk,
        .dataptr   = NULL,
    };
    SolverBlok *fblok     = cblk[0].fblokptr;
    SolverBlok *lblok     = cblk[1].fblokptr;
    size_t      size      = 0;
    int         home_node = -1;

    assert( side != PastixLUCoef );

    pastix_starpu_logger;

    /*
     * Get the right dataptr
     * coeftab if the cblk is not compressed
     * else lrblock
     */
    if ( !( cblk->cblktype & CBLK_COMPRESSED ) ) {
        size           = cblk->stride * cblk_colnbr( cblk ) * pastix_size_of( flttype );
        interf.dataptr = side == PastixLCoef ? cblk->lcoeftab : cblk->ucoeftab;
        if ( interf.dataptr != NULL ) {
            home_node = STARPU_MAIN_RAM;
        }
    }
    else {
        interf.dataptr = cblk->fblokptr->LRblock[side];
        if ( interf.dataptr != NULL )
        {
            size      = ( lblok - fblok ) * sizeof( pastix_lrblock_t );
            home_node = STARPU_MAIN_RAM;
        }
        /* FANIN / 1D */
        if ( (cblk->cblktype & CBLK_FANIN ) &&
             !(cblk->cblktype & CBLK_TASKS_2D) )
        {
            size = ( lblok - fblok ) * sizeof( pastix_lrblock_t );
        }
    }

    interf.nbblok    = lblok - fblok;
    interf.allocsize = size;

#if defined(PASTIX_STARPU_INTERFACE_DEBUG)
    fprintf( stderr,
             "cblk (%9s, size=%8zu, nbblok=%2ld )\n",
             cblk->cblktype & CBLK_COMPRESSED ? "Low-rank" : "Full-rank",
             interf.allocsize, (long)(interf.nbblok) );
#endif

    starpu_data_register( handleptr, home_node, &interf, &pastix_starpu_interface_ops );
}

void
pastix_starpu_register_ws( starpu_data_handle_t *handleptr,
                           const SolverCblk     *cblk,
                           pastix_coeftype_t     flttype )
{
    pastix_starpu_interface_t interf = {
        .id        = PASTIX_STARPU_INTERFACE_ID,
        .flttype   = flttype,
        .offset    = -1,
        .nbblok    =  0,
        .allocsize = -1,
        .cblk      = cblk,
        .dataptr   = NULL,
    };
    SolverBlok  *fblok = cblk[0].fblokptr;
    SolverBlok  *lblok = cblk[1].fblokptr;
    size_t       size  = 0;
    pastix_int_t M     = cblk->stride;
    pastix_int_t N     = cblk_colnbr( cblk );

    pastix_starpu_logger;

    /*
     * Get the correct size to allocate
     */
    interf.nbblok = lblok - fblok;
    if ( (M - N) > 0 ) /* Check for useless cases where we have only a diagonal block */
    {
        if ( (cblk->cblktype & CBLK_COMPRESSED) &&
             (cblk_getdataL( cblk ) != NULL) )
        {
            size = M * N * pastix_size_of( flttype );
            size += interf.nbblok * sizeof( pastix_lrblock_t );
        }

        if ( !(cblk->cblktype & CBLK_COMPRESSED) )
        {
            size = M * N * pastix_size_of( flttype );
        }
    }
    interf.allocsize = size;

#if defined(PASTIX_STARPU_INTERFACE_DEBUG)
    fprintf( stderr,
             "cblk (%9s, size=%8zu, nbblok=%2ld )\n",
             cblk->cblktype & CBLK_COMPRESSED ? "Low-rank" : "Full-rank",
             interf.allocsize, (long)(interf.nbblok) );
#endif

    starpu_data_register( handleptr, -1, &interf, &pastix_starpu_interface_ops );
}

void
pastix_starpu_register_blok( starpu_data_handle_t *handleptr,
                             const SolverCblk     *cblk,
                             const SolverBlok     *blok,
                             pastix_coeftype_t     flttype )
{
    pastix_starpu_interface_t interf = {
        .id        = PASTIX_STARPU_INTERFACE_ID,
        .flttype   = flttype,
        .offset    = -1,
        .nbblok    =  1,
        .allocsize = -1,
        .cblk      = cblk,
        .dataptr   = NULL,
    };
    SolverBlok  *fblok = cblk[0].fblokptr;
    SolverBlok  *lblok = cblk[1].fblokptr;
    size_t       size  = 0;
    pastix_int_t M     = blok_rownbr( blok );
    pastix_int_t N     = cblk_colnbr( cblk );

    pastix_starpu_logger;

    /* Count the number of blocks and number of rows */
    {
        const SolverBlok *cblok = blok;

        while( (cblok[0].lcblknm == cblok[1].lcblknm) &&
               (cblok[0].fcblknm == cblok[1].fcblknm) )
        {
            cblok++;
            interf.nbblok++;
            M += blok_rownbr( cblok );
        }
    }

    /*
     * Get the right dataptr
     * coeftab if the cblk is not compressed
     * else lrblock
     */
    if ( !( cblk->cblktype & CBLK_COMPRESSED ) )
    {
        size = M * N * pastix_size_of( flttype );
    }

    interf.nbblok    = lblok - fblok;
    interf.allocsize = size;

#if defined(PASTIX_STARPU_INTERFACE_DEBUG)
    fprintf( stderr,
             "cblk (%9s, size=%8zu, nbblok=%2ld )\n",
             cblk->cblktype & CBLK_COMPRESSED ? "Low-rank" : "Full-rank",
             interf.allocsize, (long)(interf.nbblok) );
#endif

    starpu_data_register( handleptr, -1, &interf, &pastix_starpu_interface_ops );
}
/**
 * @}
 */
