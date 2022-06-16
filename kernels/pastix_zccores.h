/**
 * @file pastix_zccores.h
 *
 * PaStiX kernel header.
 *
 * @copyright 2011-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Esragul Korkmaz
 * @author Gregoire Pichon
 * @author Tony Delarue
 * @author Brieuc Nicolas
 * @date 2022-07-07
 * @precisions mixed zc -> ds
 *
 */
#ifndef _pastix_zccores_h_
#define _pastix_zccores_h_

int cpucblk_zcfillin( pastix_coefside_t    side,
                      const SolverMatrix  *solvmtx,
                      const pastix_bcsc_t *bcsc,
                      pastix_int_t         itercblk );
void cpucblk_zcinit( pastix_coefside_t    side,
                     const SolverMatrix  *solvmtx,
                     const pastix_bcsc_t *bcsc,
                     pastix_int_t         itercblk,
                     const char          *directory );

#endif /* _pastix_zccores_h_ */
