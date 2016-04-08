/**
 * @file coeftab.h
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author David Goudin
 * @author Pascal Henon
 * @author Francois Pellegrini
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 **/
#ifndef _COEFTAB_H_
#define _COEFTAB_H_

void coeftab_zcompress_one( SolverCblk *cblk, double tol);
void coeftab_zuncompress_one( SolverCblk *cblk );
void coeftab_zuncompress( SolverMatrix *solvmtx );
void coeftab_zffbcsc( const SolverMatrix  *solvmtx,
                      const pastix_bcsc_t *bcsc,
                      pastix_int_t         itercblk );
void coeftab_zinitcblk( const SolverMatrix  *solvmtx,
                        const pastix_bcsc_t *bcsc,
                        pastix_int_t itercblk,
                        int fakefillin, int factoLU );
void coeftab_zdumpcblk( const SolverCblk *cblk,
                        void *array,
                        FILE *stream );
void coeftab_zdump( const SolverMatrix *solvmtx,
                    const char   *filename );

void coeftab_ccompress_one( SolverCblk *cblk, double tol);
void coeftab_cuncompress_one( SolverCblk *cblk );
void coeftab_cuncompress( SolverMatrix *solvmtx );
void coeftab_cffbcsc( const SolverMatrix  *solvmtx,
                      const pastix_bcsc_t *bcsc,
                      pastix_int_t         itercblk );
void coeftab_cinitcblk( const SolverMatrix  *solvmtx,
                        const pastix_bcsc_t *bcsc,
                        pastix_int_t itercblk,
                        int fakefillin, int factoLU );
void coeftab_cdumpcblk( const SolverCblk *cblk,
                        void *array,
                        FILE *stream );
void coeftab_cdump( const SolverMatrix *solvmtx,
                    const char   *filename );

void coeftab_dcompress_one( SolverCblk *cblk, double tol);
void coeftab_duncompress_one( SolverCblk *cblk );
void coeftab_duncompress( SolverMatrix *solvmtx );
void coeftab_dffbcsc( const SolverMatrix  *solvmtx,
                      const pastix_bcsc_t *bcsc,
                      pastix_int_t         itercblk );
void coeftab_dinitcblk( const SolverMatrix  *solvmtx,
                        const pastix_bcsc_t *bcsc,
                        pastix_int_t itercblk,
                        int fakefillin, int factoLU );
void coeftab_ddumpcblk( const SolverCblk *cblk,
                        void *array,
                        FILE *stream );
void coeftab_ddump( const SolverMatrix *solvmtx,
                    const char   *filename );

void coeftab_scompress_one( SolverCblk *cblk, double tol);
void coeftab_suncompress_one( SolverCblk *cblk );
void coeftab_suncompress( SolverMatrix *solvmtx );
void coeftab_sffbcsc( const SolverMatrix  *solvmtx,
                      const pastix_bcsc_t *bcsc,
                      pastix_int_t         itercblk );
void coeftab_sinitcblk( const SolverMatrix  *solvmtx,
                        const pastix_bcsc_t *bcsc,
                        pastix_int_t itercblk,
                        int fakefillin, int factoLU );
void coeftab_sdumpcblk( const SolverCblk *cblk,
                        void *array,
                        FILE *stream );
void coeftab_sdump( const SolverMatrix *solvmtx,
                    const char   *filename );

#endif /* _COEFTAB_H_ */
