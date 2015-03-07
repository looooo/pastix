/**
 *
 * @file symbol.h
 *
 *  PaStiX symbol structure routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author David Goudin
 * @author Francois Pelegrin
 * @author Mathieu Faverge
 * @author Pascal Henon
 * @author Pierre Ramet
 * @date 2013-06-24
 *
 **/
#ifndef _SYMBOL_H_
#define _SYMBOL_H_

/**
 * @ingroup pastix_symbol
 * @struct symbolcblk_s - Symbol column block structure.
 */
typedef struct SymbolCblk_ {
  pastix_int_t fcolnum;  /*< First column index               */
  pastix_int_t lcolnum;  /*< Last column index (inclusive)    */
  pastix_int_t bloknum;  /*< First block in column (diagonal) */
  pastix_int_t brownum;  /*< First block in row facing the diagonal block in browtab/crowtab */
} SymbolCblk;

/**
 * @ingroup pastix_symbol
 * @struct symbolblok_s - Symbol block structure.
 */
typedef struct SymbolBlok_ {
  pastix_int_t frownum;  /*< First row index            */
  pastix_int_t lrownum;  /*< Last row index (inclusive) */
  pastix_int_t cblknum;  /*< Facing column block        */
  pastix_int_t levfval;  /*< Level-of-fill value        */
} SymbolBlok;

/**
 * @ingroup pastix_symbol
 * @struct symbolmtx_s - Symbol matrix structure.
 */
typedef struct SymbolMatrix_ {
  pastix_int_t            baseval;  /*< Base value for numberings         */
  pastix_int_t            cblknbr;  /*< Number of column blocks           */
  pastix_int_t            bloknbr;  /*< Number of blocks                  */
  SymbolCblk   * restrict cblktab;  /*< Array of column blocks [+1,based] */
  SymbolBlok   * restrict bloktab;  /*< Array of blocks [based]           */
  pastix_int_t * restrict crowtab;  /*< Array of column blocks [based]    */
  pastix_int_t * restrict browtab;  /*< Array of blocks [based]           */
  pastix_int_t            nodenbr;  /*< Number of nodes in matrix         */
#ifdef STARPU_GET_TASK_CTX
  pastix_int_t            starpu_subtree_nbr;
#endif
} SymbolMatrix;

/*
**  The function prototypes.
*/

int  symbolInit          (SymbolMatrix * const symbptr);
void symbolExit          (SymbolMatrix * const symbptr);
void symbolBase          (SymbolMatrix * const symbptr, const pastix_int_t baseval);
void symbolRealloc       (SymbolMatrix * const symbptr);
int  symbolLoad          (SymbolMatrix * const symbptr, FILE * const stream);
int  symbolSave          (const SymbolMatrix * const symbptr, FILE * const stream);
int  symbolCheck         (const SymbolMatrix * const symbptr);
int  symbolDraw          (const SymbolMatrix * const symbptr, FILE * const stream);
#ifdef DOF_H
int  symbolCost          (const SymbolMatrix * const symbptr, const Dof * const deofptr,
                          const SymbolCostType typeval, double * const nnzptr, double * const opcptr);
int  symbolCosti         (const SymbolMatrix * const symbptr, const Dof * const deofptr,
                          const SymbolCostType typeval, const pastix_int_t levfval, double * const nnzptr, double * const opcptr);
int  symbolLevf          (const SymbolMatrix * const symbptr, pastix_int_t * const levfmax, pastix_int_t ** const levftab);
int  symbolTree          (const SymbolMatrix * const symbptr, const Dof * const deofptr,
                          pastix_int_t * const leafnbr, pastix_int_t * const heigmin, pastix_int_t * const heigmax,
                          double * const heigavg, double * const heigdlt);
int  symbolNonzeros      (const SymbolMatrix * const symbptr, FILE * const stream);
#endif /* DOF_H */

void symbolRustine (SymbolMatrix *       matrsymb,
                    SymbolMatrix * const matrsymb2);


pastix_int_t
symbolGetFacingBloknum(const SymbolMatrix *symbptr,
                       pastix_int_t bloksrc,
                       pastix_int_t bloknum,
                       pastix_int_t startsearch,
                       int ricar);

pastix_int_t
symbolGetNNZ(const SymbolMatrix *symbptr);

void
symbolPrintStats( const SymbolMatrix * );

#endif /* SYMBOL_H */
