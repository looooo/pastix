/************************************************************/
/**                                                        **/
/**   NAME       : solver.h                                **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the solver matrix.                  **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     28 oct 1998     **/
/**                # Version 1.0  : from : 06 jun 2002     **/
/**                                 to     06 jun 2002     **/
/**                                                        **/
/************************************************************/

#ifndef SOLVER_H
#define SOLVER_H

/*
**  The type and structure definitions.
*/

#define COMP_1D                     0
#define DIAG                        1
#define E1                          2
#define E2                          3
#define DRUNK                       4

typedef struct Task_ {
  pastix_int_t                       taskid;               /*+ COMP_1D DIAG E1 E2                                     +*/
  pastix_int_t                       prionum;              /*+ Priority value for the factorization                   +*/
  pastix_int_t                       cblknum;              /*+ Attached column block                                  +*/
  pastix_int_t                       bloknum;              /*+ Attached block                                         +*/
  pastix_int_t volatile              ftgtcnt;              /*+ Number of fan-in targets                               +*/
  pastix_int_t volatile              ctrbcnt;              /*+ Total number of contributions                          +*/
  pastix_int_t                       indnum;               /*+ For E2 (COMP_1D), index of ftgt (>0) else if local = -taskdest
                                                            For DIAG and E1 , index of btag (>0) if there is a
                                                            local one it must be the first of the chain of local E1   +*/
#if (defined PASTIX_DYNSCHED) || (defined TRACE_SOPALIN)
  pastix_int_t                       threadid;             /*+ Index of the bubble which contains the task +*/
#endif
#ifdef TRACE_SOPALIN
  pastix_int_t                       fcandnum;             /*+ First thread candidate                      +*/
  pastix_int_t                       lcandnum;		  /*+ Last thread candidate                       +*/
  pastix_int_t                       id;                   /*+ Global cblknum of the attached column block +*/
#endif
} Task;


/*+ Solver block structure. +*/

typedef struct SolverBlok_ {
  pastix_int_t                       frownum;              /*+ First row index            +*/
  pastix_int_t                       lrownum;              /*+ Last row index (inclusive) +*/
  pastix_int_t                       cblknum;              /*+ Facing column block        +*/
  pastix_int_t                       levfval;              /*+ Level-of-fill value        +*/
  pastix_int_t                       coefind;              /*+ Index in coeftab           +*/
} SolverBlok;

/*+ Solver column block structure. +*/

typedef struct SolverCblk_  {
  pastix_int_t                       fcolnum;              /*+ First column index                     +*/
  pastix_int_t                       lcolnum;              /*+ Last column index (inclusive)          +*/
  SolverBlok    *                    fblokptr;            /*+ First block in column (diagonal)       +*/
  pastix_int_t                       stride;               /*+ Column block stride                    +*/
  pastix_int_t                       procdiag;             /*+ Processor owner of diagonal block      +*/
  pastix_int_t                       gcblknum;             /*+ Global column block index              +*/
  pastix_float_t * restrict          coeftab;              /*+ Coefficients access vector             +*/
  pastix_float_t * restrict          ucoeftab;             /*+ Coefficients access vector             +*/
} SolverCblk;

/*+ Solver matrix structure. +*/

/* All data are local to one cluster */
typedef struct SolverMatrix_ {
  pastix_int_t              baseval;              /*+ Base value for numberings                         +*/

  pastix_int_t              nodenbr;              /*+ Number of nodes before dof extension              +*/
  pastix_int_t              coefnbr;              /*+ Number of coefficients (node after dof extension) +*/
  pastix_int_t              gcblknbr;             /*+ Global number of column blocks                    +*/
  pastix_int_t              cblknbr;              /*+ Number of column blocks                   +*/
  pastix_int_t              bloknbr;              /*+ Number of blocks                          +*/
  SolverCblk * restrict     cblktab;              /*+ Array of solver column blocks             +*/
  SolverBlok * restrict     bloktab;              /*+ Array of solver blocks                    +*/

#ifdef PASTIX_WITH_STARPU
  /* All this part concern halo of the local matrix
   * ie: column blocks which will:
   *  - be updated by local column blocks
   *  - update local column blocks
   */
  pastix_int_t              hcblknbr;             /*+ Number of column block in the halo        +*/
  pastix_int_t *            gcblk2halo;           /*+ Indicate the local number corresponding
                                                   *  global column block.
                                                   *  gcblk2halo[gcblk] == 0 : gcblk not local nor in halo
                                                   *                    >  0 : local cblk number + 1
                                                   *                    <  0 : - (halo cblk number + 1)
                                                   *                                            +*/
  SolverCblk   * restrict hcblktab;             /*+ Array of halo column blocks               +*/
  SolverBlok   * restrict hbloktab;             /*+ Array of halo blocks                      +*/
  pastix_int_t *          fcblknbr;               /*+ Number of fanin buffer to send or recv    +*/
  SolverCblk  ** restrict fcblktab;               /*+ Fanin column block array                  +*/
  SolverBlok  ** restrict fbloktab;               /*+ Fanin block array                         +*/
#endif

  pastix_int_t              ftgtnbr;              /*+ Number of fanintargets                    +*/
  pastix_int_t              ftgtcnt;              /*+ Number of fanintargets to receive         +*/
  FanInTarget * restrict    ftgttab;              /*+ Fanintarget access vector                 +*/

  pastix_int_t              coefmax;              /*+ Working block max size (cblk coeff 1D)    +*/
  pastix_int_t              nbftmax;              /*+ Maximum block number in ftgt              +*/
  pastix_int_t              arftmax;              /*+ Maximum block area in ftgt                +*/

  pastix_int_t              clustnum;             /*+ current processor number                  +*/
  pastix_int_t              clustnbr;             /*+ number of processors                      +*/
  pastix_int_t              procnbr;              /*+ Number of physical processor used         +*/
  pastix_int_t              thrdnbr;              /*+ Number of local computation threads       +*/
  pastix_int_t              bublnbr;              /*+ Number of local computation threads       +*/
  BubbleTree  * restrict    btree;                /*+ Bubbles tree                              +*/

  pastix_int_t              indnbr;
  pastix_int_t * restrict   indtab;
  Task * restrict           tasktab;              /*+ Task access vector                        +*/
  pastix_int_t              tasknbr;              /*+ Number of Tasks                           +*/
  pastix_int_t **           ttsktab;              /*+ Task access vector by thread              +*/
  pastix_int_t *            ttsknbr;              /*+ Number of tasks by thread                 +*/

  pastix_int_t *            proc2clust;           /*+ proc -> cluster                           +*/
  pastix_int_t              gridldim;             /*+ Dimensions of the virtual processors      +*/
  pastix_int_t              gridcdim;             /*+ grid if dense end block                   +*/
  UpDownVector              updovct;              /*+ UpDown vector                             +*/
} SolverMatrix;

pastix_int_t
sizeofsolver(const SolverMatrix *solvptr,
             pastix_int_t *iparm );

/**
 * Indicates whether a column block is in halo.
 *
 * @param datacode SolverMatrix structure.
 * @param column block SolverCblk structure to test.
 *
 * @retval API_YES if the column block is in halo.
 * @retval API_NO  if the column block is not in halo.
 */
static inline
int cblk_ishalo( SolverMatrix * datacode,
                 SolverCblk   * cblk ) {
    if ((size_t)cblk >= (size_t)datacode->hcblktab &&
        (size_t)cblk < (size_t)(datacode->hcblktab+datacode->hcblknbr)) {
        return API_YES;
    }
    return API_NO;
}

/**
 * Indicates whether a column block is a fanin column block.
 *
 * @param datacode SolverMatrix structure.
 * @param column block SolverCblk structure to test.
 *
 * @retval API_YES if the column block is a fanin column block.
 * @retval API_NO  if the column block is not a fanin column block.
 */
static inline
int cblk_isfanin( SolverMatrix * datacode,
                  SolverCblk   * cblk ) {
    pastix_int_t clustnum;
    for (clustnum = 0; clustnum < datacode->clustnbr; clustnum++) {
        if ((size_t)cblk >= (size_t)datacode->fcblktab[clustnum] &&
            (size_t)cblk < (size_t)(datacode->fcblktab[clustnum] +
                                    datacode->fcblknbr[clustnum])) {
            return API_YES;
        }
    }
    return API_NO;
}


/**
 * Indicates whether a column block is a local column block.
 *
 * @param datacode SolverMatrix structure.
 * @param column block SolverCblk structure to test.
 *
 * @retval API_YES if the column block is a local column block.
 * @retval API_NO  if the column block is not a local column block.
 */
static inline
int cblk_islocal( SolverMatrix * datacode,
                  SolverCblk   * cblk ) {
    if ((size_t)cblk >= (size_t)datacode->cblktab &&
        (size_t)cblk < (size_t)(datacode->cblktab+datacode->cblknbr)) {
        return API_YES;
    }
    return API_NO;
}


/**
 * Get the index of a local column block.
 *
 * @param datacode SolverMatrix structure.
 * @param column block SolverCblk structure to test.
 *
 * @returns the index of the column block.
 */
static inline
pastix_int_t cblk_getnum( SolverMatrix * datacode,
                          SolverCblk   * cblk ) {
    assert(cblk_islocal(datacode, cblk) == API_YES);
    return cblk - datacode->cblktab;
}


/**
 * Get the index of a halo column block.
 *
 * @param datacode SolverMatrix structure.
 * @param column block SolverCblk structure to test.
 *
 * @returns the index of the column block.
 */
static inline
pastix_int_t hcblk_getnum( SolverMatrix * datacode,
                           SolverCblk   * cblk ) {
    assert(cblk_ishalo(datacode, cblk) == API_YES);
    return cblk - datacode->hcblktab;
}


/**
 * Get the index of a fanin column block.
 *
 * @param datacode SolverMatrix structure.
 * @param column block SolverCblk structure to test.
 *
 * @returns the index of the column block.
 */
static inline
pastix_int_t fcblk_getnum( SolverMatrix * datacode,
                           SolverCblk   * cblk,
                           pastix_int_t   procnum ) {
    assert(cblk_isfanin(datacode, cblk) == API_YES);
    return cblk - datacode->fcblktab[procnum];
}

/**
 * Get the number of columns of a column block.
 *
 * @param column block SolverCblk structure.
 *
 * @returns the number of columns in the column block.
 */
static inline
pastix_int_t cblk_colnbr( SolverCblk * cblk ) {
    return cblk->lcolnum - cblk->fcolnum + 1;
}
/**
 * Indicate if a blok is included inside an other block.
 * i.e. indicate if the row range of the first block is included in the
 * one of the second.
 *
 * @param first block SolverBlok structure to test.
 * @param second block SolverBlok structure to test.
 *
 * @retval true   if the first block is     included in the second one.
 * @retval false  if the first block is not included in the second one.
 */
#  if defined(NAPA_SOPALIN)
static inline int is_block_inside_fblock( SolverBlok *blok, SolverBlok *fblok ) {
    return (((blok->frownum >= fblok->frownum) &&
             (blok->lrownum <= fblok->lrownum)) ||
            ((blok->frownum <= fblok->frownum) &&
             (blok->lrownum >= fblok->lrownum)) ||
            ((blok->frownum <= fblok->frownum) &&
             (blok->lrownum >= fblok->frownum)) ||
            ((blok->frownum <= fblok->lrownum) &&
             (blok->lrownum >= fblok->lrownum)));
}
#  else
static inline int is_block_inside_fblock( SolverBlok *blok, SolverBlok *fblok ) {
    return ((blok->frownum >= fblok->frownum) &&
            (blok->lrownum <= fblok->lrownum));
}
#  endif /* defined(NAPA_SOPALIN) */

#endif /* SOLVER_H */
