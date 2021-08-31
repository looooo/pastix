/**
 *
 * @file pastix_starpu_interface.c
 *
 * Interface used by StarPU to handle factorization.
 *
 * @copyright 2021-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.1
 * @author Nolan Bredel
 * @author Mathieu Faverge
 * @date 2021-07-20
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

static inline void
psi_register_data_handle( starpu_data_handle_t handle, unsigned home_node, void *data_interface )
{
    pastix_starpu_interface_t *interface = (pastix_starpu_interface_t *)data_interface;
    unsigned                   node;

    for ( node = 0; node < STARPU_MAXNODES; node++ ) {
        pastix_starpu_interface_t *local_interface =
            (pastix_starpu_interface_t *)starpu_data_get_interface_on_node( handle, node );

        memcpy( local_interface, interface, sizeof( pastix_starpu_interface_t ) );

        if ( node != home_node ) {
            local_interface->dataptr = NULL;
        }
    }
}

static inline starpu_ssize_t
psi_allocate_data_on_node( void *data_interface, unsigned node )
{
    pastix_starpu_interface_t *interface = (pastix_starpu_interface_t *)data_interface;
    starpu_ssize_t             allocated_memory;
    uintptr_t                  addr = 0;
    uintptr_t                  handle;

    allocated_memory = interface->allocsize;
    if ( allocated_memory <= 0 ) {
        return 0;
    }

    handle = starpu_malloc_on_node( node, allocated_memory );
    if ( !handle ) {
        return -ENOMEM;
    }

    if ( starpu_node_get_kind( node ) != STARPU_OPENCL_RAM ) {
        addr = handle;
    }

    /* update the data properly */
    interface->dataptr = (void *)addr;

    /* Allocate the workspace for the low-rank blocks */
    if ( interface->cblk->cblktype & CBLK_COMPRESSED )
    {
        SolverBlok       *blok    = interface->cblk->fblokptr;
        pastix_lrblock_t *LRblock = interface->dataptr;
        int               offset  = pastix_imax( 0, interface->offset );
        int               i, ncols, M;

        assert( node == STARPU_MAIN_RAM );

        ncols = cblk_colnbr( interface->cblk );
        blok += offset;
        for ( i = 0; i < interface->nbblok; i++, blok++, LRblock++ ) {
            M = blok_rownbr( blok );

            /* Allocate the LR block to its tight space */
            switch ( interface->flttype ) {
                case PastixComplex64:
                    core_zlralloc( M, ncols, -1, LRblock );
                    break;
                case PastixComplex32:
                    core_clralloc( M, ncols, -1, LRblock );
                    break;
                case PastixDouble:
                    core_dlralloc( M, ncols, -1, LRblock );
                    break;
                case PastixFloat:
                    core_slralloc( M, ncols, -1, LRblock );
                    break;
                default:
                    assert( 0 );
            }
        }
    }

    return allocated_memory;
}

static inline void
psi_free_data_on_node( void *data_interface, unsigned node )
{
    pastix_starpu_interface_t *interface = (pastix_starpu_interface_t *)data_interface;

    SolverCblk *cblk = interface->cblk;

    if ( cblk->cblktype & CBLK_COMPRESSED )
    {
        pastix_lrblock_t *LRblock = interface->dataptr;
        int               i;

        for ( i = 0; i < interface->nbblok; i++, LRblock++ ) {
            core_zlrfree( LRblock );
        }
    }

    starpu_free_on_node( node, (uintptr_t)interface->dataptr, interface->allocsize );

    interface->dataptr = NULL;
}

static inline void
psi_init( void *data_interface )
{
    pastix_starpu_interface_t *interface = data_interface;
    interface->id                        = PASTIX_STARPU_INTERFACE_ID;
    interface->allocsize                 = -1;
}

static inline void *
psi_to_pointer( void *data_interface, unsigned node )
{
    (void)node;
    pastix_starpu_interface_t *interface = (pastix_starpu_interface_t *)data_interface;

    return interface->dataptr;
}

static inline int
psi_pointer_is_inside( void *data_interface, unsigned node, void *ptr )
{
    (void)node;
    (void)data_interface;
    (void)ptr;
    return 0;
}

static inline size_t
psi_get_size( starpu_data_handle_t handle )
{
    pastix_starpu_interface_t *interface =
        starpu_data_get_interface_on_node( handle, STARPU_MAIN_RAM );
    size_t      size;
    SolverCblk *cblk = interface->cblk;

    if ( !( cblk->cblktype & CBLK_COMPRESSED ) ) {
        if ( interface->offset == -1 ) {
            pastix_int_t ncols = cblk_colnbr( cblk );
            size               = cblk->stride * ncols * pastix_size_of( interface->flttype );
        }
        else {
            size = interface->allocsize;
        }
    }
    else {
        size = interface->nbblok * sizeof( pastix_lrblock_t );
    }

#ifdef STARPU_DEBUG
    STARPU_ASSERT_MSG( interface->id == PASTIX_STARPU_INTERFACE_ID,
                       "Error. The given data is not a pastix interface for starpu." );
#endif

    return size;
}

static inline size_t
psi_get_alloc_size( starpu_data_handle_t handle )
{
    pastix_starpu_interface_t *interface =
        starpu_data_get_interface_on_node( handle, STARPU_MAIN_RAM );

#ifdef STARPU_DEBUG
    STARPU_ASSERT_MSG( interface->id == PASTIX_STARPU_INTERFACE_ID,
                       "Error. The given data is not a ." );
#endif

    return interface->allocsize;
}

static inline uint32_t
psi_footprint( starpu_data_handle_t handle )
{
    pastix_starpu_interface_t *interface =
        starpu_data_get_interface_on_node( handle, STARPU_MAIN_RAM );
    SolverCblk *cblk = interface->cblk;

    return starpu_hash_crc32c_be( cblk->gcblknum, interface->offset + 1 );
}

static inline uint32_t
psi_alloc_footprint( starpu_data_handle_t handle )
{
    pastix_starpu_interface_t *interface =
        starpu_data_get_interface_on_node( handle, STARPU_MAIN_RAM );

    return starpu_hash_crc32c_be( interface->allocsize, 0 );
}

static inline int
psi_compare( void *data_interface_a, void *data_interface_b )
{
    pastix_starpu_interface_t *pastix_interface_a = (pastix_starpu_interface_t *)data_interface_a;
    pastix_starpu_interface_t *pastix_interface_b = (pastix_starpu_interface_t *)data_interface_b;

    void *solva = ( pastix_interface_a->cblk );
    void *solvb = ( pastix_interface_b->cblk );

    /* Two interfaces are considered compatible if they point to the same solver */
    return ( solva == solvb );
}

static inline int
psi_alloc_compare( void *data_interface_a, void *data_interface_b )
{
    pastix_starpu_interface_t *pastix_interface_a = (pastix_starpu_interface_t *)data_interface_a;
    pastix_starpu_interface_t *pastix_interface_b = (pastix_starpu_interface_t *)data_interface_b;

    /* Two matrices are considered compatible if they have the same allocated size */
    return ( pastix_interface_a->allocsize == pastix_interface_b->allocsize );
}

static inline void
psi_display( starpu_data_handle_t handle, FILE *f )
{
    pastix_starpu_interface_t *interface =
        (pastix_starpu_interface_t *)starpu_data_get_interface_on_node( handle, STARPU_MAIN_RAM );

    SolverCblk *cblk = interface->cblk;

    if ( interface->offset == -1 ) {
        fprintf( f, "Cblk%ld", (long)( cblk->gcblknum ) );
    }
    else {
        fprintf( f, "Cblk%ldBlok%ld", (long)( cblk->gcblknum ), (long)( interface->offset ) );
    }
}

static inline size_t
psi_compute_size_lr( pastix_starpu_interface_t *interface )
{
    assert( interface->cblk->cblktype & CBLK_COMPRESSED );

    size_t            elemsize = pastix_size_of( interface->flttype );
    pastix_lrblock_t *LRblock  = interface->dataptr;
    pastix_int_t      N        = cblk_colnbr( interface->cblk );
    pastix_int_t      suv, rk, M;

    SolverBlok *blok = interface->cblk->fblokptr + pastix_imax( 0, interface->offset );

    size_t size = 0;
    int    i;
    for ( i = 0; i < interface->nbblok; i++, blok++, LRblock++ ) {
        M  = blok_rownbr( blok );
        rk = LRblock->rk;
        if ( rk != -1 ) {
            suv = rk * ( M + N );
        }
        else {
            suv = M * N;
        }
        size += sizeof( int ) + suv * elemsize;
    }
    return size;
}

static inline size_t
psi_compute_size_fr( pastix_starpu_interface_t *interface )
{
    assert( !( interface->cblk->cblktype & CBLK_COMPRESSED ) );

    return interface->allocsize;
}

static inline void
psi_pack_lr( pastix_starpu_interface_t *interface, void **ptr )
{
    assert( interface->cblk->cblktype & CBLK_COMPRESSED );

    size_t            elemsize = pastix_size_of( interface->flttype );
    char             *tmp      = *ptr;
    pastix_lrblock_t *LRblock  = interface->dataptr;
    int               N        = cblk_colnbr( interface->cblk );
    int               j        = 0;
    int               M, rk;

    SolverBlok *blok = interface->cblk->fblokptr + pastix_imax( 0, interface->offset );

    for ( ; j < interface->nbblok; j++, blok++, LRblock++ ) {
        M  = blok_rownbr( blok );
        rk = LRblock->rk;

        /* Get the rank */
        assert( rk >= -1 );
        memcpy( tmp, &rk, sizeof( int ) );
        tmp += sizeof( int );

        /* Copy the data */
        if ( rk != -1 ) {
            /* Pack U */
            memcpy( tmp, LRblock->u, M * rk * elemsize );
            tmp += M * rk * elemsize;

            /* Pack V */
            memcpy( tmp, LRblock->v, N * rk * elemsize );
            tmp += N * rk * elemsize;
        }
        else {
            /* Pack the full block */
            memcpy( tmp, LRblock->u, M * N * elemsize );
            tmp += M * N * elemsize;
        }
    }
}

static inline void
psi_pack_fr( pastix_starpu_interface_t *interface, void **ptr, starpu_ssize_t *count )
{
    assert( !( interface->cblk->cblktype & CBLK_COMPRESSED ) );

    memcpy( *ptr, interface->dataptr, *count );
}

static inline int
psi_pack_data( starpu_data_handle_t handle, unsigned node, void **ptr, starpu_ssize_t *count )
{
    STARPU_ASSERT( starpu_data_test_if_allocated_on_node( handle, node ) );

    pastix_starpu_interface_t *interface =
        (pastix_starpu_interface_t *)starpu_data_get_interface_on_node( handle, node );

    SolverCblk *cblk = interface->cblk;

    if ( cblk->cblktype & CBLK_COMPRESSED ) {
        *count = psi_compute_size_lr( interface );
    }
    else {
        *count = psi_compute_size_fr( interface );
    }

    if ( ptr != NULL ) {
        *ptr = (void *)starpu_malloc_on_node_flags( node, *count, 0 );

        if ( cblk->cblktype & CBLK_COMPRESSED ) {
            psi_pack_lr( interface, ptr );
        }
        else {
            psi_pack_fr( interface, ptr, count );
        }
    }

    return 0;
}

static inline void
psi_unpack_lr( pastix_starpu_interface_t *interface, void *ptr )
{
    assert( interface->cblk->cblktype & CBLK_COMPRESSED );

    size_t            elemsize = pastix_size_of( interface->flttype );
    SolverBlok       *blok     = interface->cblk->fblokptr + pastix_imax( 0, interface->offset );
    pastix_lrblock_t *LRblock  = interface->dataptr;
    char             *tmp      = ptr;
    int               N        = cblk_colnbr( interface->cblk );
    int               rk;

    int i = 0;
    for ( ; i < interface->nbblok; i++, blok++, LRblock++ ) {
        int M = blok_rownbr( blok );

        /* Get the rank */
        memcpy( &rk, tmp, sizeof( int ) );
        tmp += sizeof( int );

        /* Allocate the LR block to its tight space */
        switch ( interface->flttype ) {
            case PastixComplex64:
                core_zlralloc( M, N, rk, LRblock );
                break;
            case PastixComplex32:
                core_clralloc( M, N, rk, LRblock );
                break;
            case PastixDouble:
                core_dlralloc( M, N, rk, LRblock );
                break;
            case PastixFloat:
                core_slralloc( M, N, rk, LRblock );
                break;
            default:
                assert( 0 );
        }

        /* Copy the data */
        if ( rk != -1 ) {
            /* Unpack U */
            memcpy( LRblock->u, tmp, M * rk * elemsize );
            tmp += M * rk * elemsize;

            /* Unpack V */
            memcpy( LRblock->v, tmp, N * rk * elemsize );
            tmp += N * rk * elemsize;
        }
        else {
            /* Unpack the full block */
            memcpy( LRblock->u, tmp, M * N * elemsize );
            tmp += M * N * elemsize;
        }
    }
}

static inline void
psi_unpack_fr( pastix_starpu_interface_t *interface, void *ptr, size_t count )
{
    memcpy( interface->dataptr, ptr, count );
}

static inline int
psi_peek_data( starpu_data_handle_t handle, unsigned node, void *ptr, size_t count )
{
    pastix_starpu_interface_t *interface =
        (pastix_starpu_interface_t *)starpu_data_get_interface_on_node( handle, node );

    STARPU_ASSERT( starpu_data_test_if_allocated_on_node( handle, node ) );

    if ( interface->cblk->cblktype & CBLK_COMPRESSED ) {
        psi_unpack_lr( interface, ptr );
    }
    else {
        psi_unpack_fr( interface, ptr, count );
    }

    return 0;
}

static inline int
psi_unpack_data( starpu_data_handle_t handle, unsigned node, void *ptr, size_t count )
{
    psi_peek_data( handle, node, ptr, count );

    /* Free the received information */
    starpu_free_on_node_flags( node, (uintptr_t)ptr, count, 0 );

    return 0;
}

static inline starpu_ssize_t
psi_describe( void *data_interface, char *buf, size_t size )
{
    pastix_starpu_interface_t *interface = (pastix_starpu_interface_t *)data_interface;
    SolverCblk                *cblk      = interface->cblk;

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

    SolverCblk *cblk = pastix_src->cblk;
    int         ret  = 0;
    size_t      size = 0;

    if ( cblk->cblktype & CBLK_COMPRESSED ) {
        size = pastix_src->allocsize;
    }
    else {
        pastix_int_t ncols   = cblk_colnbr( cblk );
        size_t       coefnbr = cblk->stride * ncols;

        size = coefnbr * pastix_size_of( pastix_src->flttype );
    }

    if ( starpu_interface_copy( (uintptr_t)pastix_src->dataptr, 0, src_node,
                                (uintptr_t)pastix_dst->dataptr, 0, dst_node,
                                size,
                                async_data ) ) {
        ret = -EAGAIN;
    }
    starpu_interface_data_copy( src_node, dst_node, size );

    return ret;
}

static const struct starpu_data_copy_methods psi_copy_methods = {
    .any_to_any = psi_copy_any_to_any,
};

struct starpu_data_interface_ops pastix_starpu_interface_ops = {
    .init                  = psi_init,
    .register_data_handle  = psi_register_data_handle,
    .allocate_data_on_node = psi_allocate_data_on_node,
    .free_data_on_node     = psi_free_data_on_node,
    .to_pointer            = psi_to_pointer,
    .pointer_is_inside     = psi_pointer_is_inside,
    .get_size              = psi_get_size,
    .get_alloc_size        = psi_get_alloc_size,
    .footprint             = psi_footprint,
    .alloc_footprint       = psi_alloc_footprint,
    .compare               = psi_compare,
    .alloc_compare         = psi_alloc_compare,
    .display               = psi_display,
    .pack_data             = psi_pack_data,
#if defined( HAVE_STARPU_DATA_PEEK )
    .peek_data             = psi_peek_data,
#endif
    .unpack_data           = psi_unpack_data,
    .describe              = psi_describe,
    .copy_methods          = &psi_copy_methods,
    .interfaceid           = STARPU_UNKNOWN_INTERFACE_ID,
    .interface_size        = sizeof( pastix_starpu_interface_t ),
    .name                  = "PASTIX_STARPU_INTERFACE"
};

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

void
pastix_starpu_register( starpu_data_handle_t *handleptr,
                        int                   home_node,
                        SolverCblk           *cblk,
                        pastix_coefside_t     side,
                        pastix_coeftype_t     flttype )
{
    pastix_starpu_interface_t interface = {
        .id        = PASTIX_STARPU_INTERFACE_ID,
        .flttype   = flttype,
        .offset    = -1,
        .nbblok    =  1,
        .allocsize = -1,
        .cblk      = cblk,
        .dataptr   = NULL,
    };
    SolverBlok *fblok = cblk[0].fblokptr;
    SolverBlok *lblok = cblk[1].fblokptr;
    size_t      size;

    assert( side != PastixLUCoef );

    /*
     * Get the right dataptr
     * coeftab if the cblk is not compressed
     * else lrblock
     */
    if ( !( cblk->cblktype & CBLK_COMPRESSED ) ) {
        size              = cblk->stride * cblk_colnbr( cblk ) * pastix_size_of( flttype );
        interface.dataptr = side == PastixLCoef ? cblk->lcoeftab : cblk->ucoeftab;
    }
    else {
        size              = ( lblok - fblok ) * sizeof( pastix_lrblock_t );
        interface.dataptr = cblk->fblokptr->LRblock[side];
    }

    interface.nbblok    = lblok - fblok;
    interface.allocsize = size;

    starpu_data_register( handleptr, home_node, &interface, &pastix_starpu_interface_ops );
}
