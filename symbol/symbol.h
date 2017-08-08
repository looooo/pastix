/**
 *
 * @file symbol.h
 *
 * PaStiX symbol structure routines
 *
 * @copyright 2004-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author David Goudin
 * @author Francois Pellegrini
 * @author Mathieu Faverge
 * @author Pascal Henon
 * @author Pierre Ramet
 * @date 2013-06-24
 *
 * @addtogroup pastix_symbol
 * @{
 *   @brief Functions to generate and manipulate the symbolic factorization
 *   structure
 *
 *   This module provides the set of function to generate the symbolic
 *   factorization structure based on a given graph, and an associated
 *   ordering. The symbolic structure is described in the symbol_matrix_t
 *   structure, and it can be generated through two different algorithms: Fax or
 *   Kass. The first one is used when no amalgamation is required. This is the
 *   case when the ordering comes from Scotch for example. The second one, is
 *   used when the elementary elimination tree has been rediscovered, and
 *   amalgamation needs to be performed to improve solver efficiency.
 *
 **/
#ifndef _SYMBOL_H_
#define _SYMBOL_H_

/**
 * @brief Symbol column block structure.
 */
typedef struct symbol_cblk_s {
    pastix_int_t fcolnum;    /**< First column index               */
    pastix_int_t lcolnum;    /**< Last column index (inclusive)    */
    pastix_int_t bloknum;    /**< First block in column (diagonal) */
    pastix_int_t brownum;    /**< First block in row facing the diagonal block in browtab, 0-based */
#if defined(PASTIX_SYMBOL_DUMP_SYMBMTX)
    pastix_int_t split_cblk;
#endif
} SymbolCblk;

/**
 * @brief Symbol block structure.
 */
typedef struct symbol_blok_s {
    pastix_int_t frownum; /**< First row index            */
    pastix_int_t lrownum; /**< Last row index (inclusive) */
    pastix_int_t lcblknm; /**< Local column block         */
    pastix_int_t fcblknm; /**< Facing column block        */
} SymbolBlok;

/**
 * @brief Symbol matrix structure.
 *
 * This structure describes the symbolic block structure of the factorized
 * matrix L, U is never stored as it is a symmetry of L. This structure is
 * global and replicated on all processes. The default way to number the block
 * is the CSC format where block are continuously number per column, the browtab
 * array stores the CSR representation of the L structure to have a faster
 * access to the list of blocks updating a column block.
 *
 */
typedef struct symbol_matrix_s {
    pastix_int_t            baseval;  /**< Base value for numbering                */
    pastix_int_t            dof;      /**< Degrees of freedom per node
                                          (constant if > 0, unconstant if 0 (not implemented)) */
    pastix_int_t            cblknbr;  /**< Number of column blocks                 */
    pastix_int_t            bloknbr;  /**< Number of blocks                        */
    pastix_int_t            nodenbr;  /**< Number of node in the compressed symbol */
    pastix_int_t            schurfcol;/**< First column of the schur complement    */
    SymbolCblk   * restrict cblktab;  /**< Array of column blocks [+1,based]       */
    SymbolBlok   * restrict bloktab;  /**< Array of blocks in CSC format [based]   */
    pastix_int_t * restrict browtab;  /**< Array of blocks in CSR format [based]   */
} symbol_matrix_t;

/**
 * @name Symbol basic subroutines
 * @{
 */
void symbolInit       (      symbol_matrix_t *symbptr);
void symbolExit       (      symbol_matrix_t *symbptr);
void symbolBase       (      symbol_matrix_t *symbptr, const pastix_int_t baseval);
void symbolRealloc    (      symbol_matrix_t *symbptr);
int  symbolCheck      (const symbol_matrix_t *symbptr);

/**
 * @}
 * @name Symbol IO subroutines
 * @{
 */
int  symbolSave       (const symbol_matrix_t *symbptr, FILE *stream);
int  symbolLoad       (      symbol_matrix_t *symbptr, FILE *stream);
int  symbolDraw       (const symbol_matrix_t *symbptr, FILE *stream);

/**
 * @}
 * @name Symbol statistical information subroutines
 * @{
 */
void         symbolPrintStats ( const symbol_matrix_t *symbptr );
pastix_int_t symbolGetNNZ     ( const symbol_matrix_t *symbptr );
void         symbolGetFlops   ( const symbol_matrix_t *symbmtx,
                                pastix_coeftype_t flttype, pastix_factotype_t factotype,
                                double *thflops, double *rlflops );
void         symbolGetTimes   ( const symbol_matrix_t *symbmtx,
                                pastix_coeftype_t flttype, pastix_factotype_t factotype,
                                double *cblkcost, double *blokcost );

/**
 * @}
 * @name Symbol reordering subroutines
 * @{
 */
void         symbolReordering( const symbol_matrix_t *, pastix_order_t *, pastix_int_t, int );
void         symbolReorderingPrintComplexity( const symbol_matrix_t *symbptr );

/**
 * @}
 * @name Symbol construction subroutines
 * @{
 */
int          symbolFaxGraph  ( symbol_matrix_t         *symbptr,
                               const pastix_int_t    vertnbr,
                               const pastix_int_t   *verttab,
                               const pastix_int_t   *edgetab,
                               const pastix_order_t *ordeptr );
int          symbolKass      ( int             verbose,
                               int             ilu,
                               int             levelk,
                               int             rat_cblk,
                               int             rat_blas,
                               symbol_matrix_t   *symbmtx,
                               pastix_graph_t *graph,
                               pastix_order_t *orderptr,
                               MPI_Comm        pastix_comm );
void         symbolRustine   ( symbol_matrix_t *symbptr, symbol_matrix_t *symbptr2 );
void         symbolBuildRowtab( symbol_matrix_t *symbptr );
pastix_int_t symbolGetFacingBloknum( const symbol_matrix_t *symbptr,
                                     pastix_int_t bloksrc,
                                     pastix_int_t bloknum,
                                     pastix_int_t startsearch,
                                     int ricar );
/**
 * @}
 * @} End of pastix_symbol group
 */

#endif /* SYMBOL_H */
