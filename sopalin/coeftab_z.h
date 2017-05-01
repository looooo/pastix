/**
 * @file coeftab_z.h
 *
 * Precision dependent coeficient array header.
 *
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author David Goudin
 * @author Pascal Henon
 * @author Francois Pellegrini
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 * @precisions normal z -> s d c
 *
 * @addtogroup coeftab
 * @{
 *
 **/
#ifndef _coeftab_z_h_
#define _coeftab_z_h_

/**
 *    @name PastixComplex64 compression/uncompression routines
 *    @{
 */
pastix_int_t coeftab_zcompress_one  ( SolverCblk *cblk, pastix_lr_t lowrank );
void         coeftab_zalloc_one     ( SolverCblk *cblk );
void         coeftab_zuncompress_one( SolverCblk *cblk, int factoLU );
pastix_int_t coeftab_zmemory_one    ( const SolverCblk *cblk, int factoLU );

pastix_int_t coeftab_zcompress  ( SolverMatrix *solvmtx );
void         coeftab_zuncompress( SolverMatrix *solvmtx );
pastix_int_t coeftab_zmemory    ( const SolverMatrix *solvmtx );

/**
 *    @}
 *    @name PastixComplex64 initialization routines
 *    @{
 */
void coeftab_zffbcsc  ( const SolverMatrix  *solvmtx,
                        const pastix_bcsc_t *bcsc,
                        pastix_int_t         itercblk );
void coeftab_zinitcblk( const SolverMatrix  *solvmtx,
                        const pastix_bcsc_t *bcsc,
                        pastix_int_t itercblk,
                        int fakefillin, int factoLU );

/**
 *    @}
 *    @name PastixComplex64 Schur routines
 *    @{
 */
void coeftab_zgetschur_one_fullrank( const SolverCblk *cblk, int upper_part,
                                     pastix_complex64_t *S, pastix_int_t lds );
void coeftab_zgetschur_one_lowrank ( const SolverCblk *cblk, int upper_part,
                                     pastix_complex64_t *S, pastix_int_t lds );
void coeftab_zgetschur             ( const SolverMatrix *solvmtx,
                                     pastix_complex64_t *S, pastix_int_t lds );

/**
 *    @}
 *    @name PastixComplex64 debug routines
 *    @{
 */
void coeftab_zdumpcblk( const SolverCblk   *cblk,
                        const void         *array,
                        FILE               *stream );
void coeftab_zdump    ( const SolverMatrix *solvmtx,
                        const char         *filename );
int  coeftab_zdiff    ( const SolverMatrix *solvA,
                        SolverMatrix       *solvB );

/**
 *    @}
 */
#endif /* _coeftab_z_h_ */

/**
 * @}
 */
