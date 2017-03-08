/**
 *
 * @file symbol.h
 *
 * PaStiX symbol structure routines
 *
 * @copyright (c) 2004-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                          Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author David Goudin
 * @author Francois Pellegrini
 * @author Mathieu Faverge
 * @author Pascal Henon
 * @author Pierre Ramet
 * @date 2013-06-24
 *
 * @defgroup pastix_symbol Symbolic Factorization
 * @ingroup  pastix_analyze
 * @brief Functions to generate and manipulate the symbolic factorization
 * structure
 *
 * This module provides the set of function to generate the symbolic
 * factorization structure based on a given graph, and an associated
 * ordering. The symbolic structure is described in the SymbolMatrix structure,
 * and it can be generated through two different algorithms: Fax or Kass. The
 * first one is used when no amalgamation is required. This is the case when the
 * ordering comes from Scotch for example. The second one, is used when the
 * elementary elimination tree has been rediscovered, and amalgamation needs to
 * be performed to improve solver efficiency.
 *
 * @{
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
} SymbolMatrix;

/**
 * @name Symbol basic subroutines
 * @{
 */
void symbolInit       (      SymbolMatrix *symbptr);
void symbolExit       (      SymbolMatrix *symbptr);
void symbolBase       (      SymbolMatrix *symbptr, const pastix_int_t baseval);
void symbolRealloc    (      SymbolMatrix *symbptr);
int  symbolCheck      (const SymbolMatrix *symbptr);

/**
 * @}
 * @name Symbol IO subroutines
 * @{
 */
int  symbolSave       (const SymbolMatrix *symbptr, FILE *stream);
int  symbolLoad       (      SymbolMatrix *symbptr, FILE *stream);
int  symbolDraw       (const SymbolMatrix *symbptr, FILE *stream);

/**
 * @}
 * @name Symbol statistical information subroutines
 * @{
 */
void         symbolPrintStats ( const SymbolMatrix *symbptr );
pastix_int_t symbolGetNNZ     ( const SymbolMatrix *symbptr );
void         symbolGetFlops   ( const SymbolMatrix *symbmtx,
                                pastix_coeftype_t flttype, pastix_factotype_t factotype,
                                double *thflops, double *rlflops );
void         symbolGetTimes   ( const SymbolMatrix *symbmtx,
                                pastix_coeftype_t flttype, pastix_factotype_t factotype,
                                double *cblkcost, double *blokcost );

/**
 * @}
 * @name Symbol reordering subroutines
 * @{
 */
void         symbolReordering( const SymbolMatrix *, Order *, pastix_int_t, int );
void         symbolReorderingPrintComplexity( const SymbolMatrix *symbptr );

/**
 * @}
 * @name Symbol construction subroutines
 * @{
 */
int          symbolFaxGraph  ( SymbolMatrix * const symbptr,
                               const pastix_int_t   vertnbr,
                               const pastix_int_t * verttab,
                               const pastix_int_t * edgetab,
                               const Order  * const ordeptr );
int          symbolKass      ( int             verbose,
                               int             ilu,
                               int             levelk,
                               int             rat_cblk,
                               int             rat_blas,
                               SymbolMatrix   *symbmtx,
                               pastix_graph_t *graph,
                               Order          *orderptr,
                               MPI_Comm        pastix_comm );
void         symbolRustine   ( SymbolMatrix *symbptr, SymbolMatrix *symbptr2 );
void         symbolBuildRowtab( SymbolMatrix *symbptr );
pastix_int_t symbolGetFacingBloknum( const SymbolMatrix *symbptr,
                                     pastix_int_t bloksrc,
                                     pastix_int_t bloknum,
                                     pastix_int_t startsearch,
                                     int ricar );
/**
 * @}
 * @} End of pastix_symbol group
 */

#endif /* SYMBOL_H */
