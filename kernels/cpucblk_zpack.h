/**
 *
 * @file cpucblk_zpack.h
 *
 * Precision dependent routines to pack and unpack cblks.
 *
 * @copyright 2021-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.2
 * @author Nolan Bredel
 * @date 2023-07-21
 *
 * @precisions normal z -> s d c
 *
 **/

#ifndef _cpucblk_zpack_h_
#define _cpucblk_zpack_h_

size_t cpublok_zcompute_size_lr( pastix_coefside_t  side,
                                 pastix_int_t       N,
                                 const SolverBlok  *blok );

pastix_uint_t cpucblk_zcompute_size_lr( pastix_coefside_t  side,
                                        const SolverCblk  *cblk );

size_t cpucblk_zcompute_size( pastix_coefside_t  side,
                              const SolverCblk  *cblk );

char *cpublok_zpack_lr( pastix_coefside_t  side,
                        pastix_uint_t      N,
                        const SolverBlok  *blok,
                        char *             buffer );

void *cpucblk_zpack_lr( pastix_coefside_t  side,
                        SolverCblk        *cblk,
                        size_t             size );

char *cpublok_zunpack_lr( pastix_coefside_t  side,
                          pastix_int_t       N,
                          SolverBlok        *blok,
                          char              *buffer );

void cpucblk_zunpack_lr( pastix_coefside_t  side,
                         SolverCblk        *cblk,
                         void              *buffer );

void *cpucblk_zpack_fr( pastix_coefside_t  side,
                        const SolverCblk  *cblk );

void cpucblk_zunpack_fr( pastix_coefside_t   side,
                         SolverCblk         *cblk,
                         pastix_complex64_t *buffer );

void *cpucblk_zpack( pastix_coefside_t  side,
                     SolverCblk        *cblk,
                     size_t             size );

void cpucblk_zunpack( pastix_coefside_t  side,
                      SolverCblk        *cblk,
                      void              *buffer );

#endif /* _cpucblk_zpack_h_ */
