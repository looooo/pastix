/**
 * @file solver.h
 *
 * PaStiX solver structure header.
 *
 * @copyright 2004-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author David Goudin
 * @author Pascal Henon
 * @author Francois Pellegrini
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Xavier Lacoste
 * @author Esragul Korkmaz
 * @author Gregoire Pichon
 * @author Tony Delarue
 * @author Alycia Lisito
 * @author Brieuc Nicolas
 * @author Nolan Bredel
 * @date 2024-07-05
 *
 **/
#ifndef _solver_h_
#define _solver_h_

#ifndef DOXYGEN_SHOULD_SKIP_THIS
struct blendctrl_s;
typedef struct blendctrl_s BlendCtrl;

struct simuctrl_s;
typedef struct simuctrl_s SimuCtrl;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#include "pastix_lowrank.h"

/**
 * @name Cblk properties
 * @{
 *  The type and structure definitions.
 *  Define the mask for the cblks in the cblktype field:
 *   - 1st bit: The cblk is a local fanin accumulating the contributions for a remote cblk
 *   - 2nd bit: The cblk is stored in a 2D layout fashion as in a tiled matrix, otherwise the standard 1D lapack layout is used
 *   - 3rd bit: The cblk generates 2D granularity tasks, instead of a single 1D tasks that perform factorization, solves and updates
 *   - 4th bit: The cblk is compressed in Low-Rank (implies CBLK_LAYOUT_2D), otherwise it is stored in dense
 *   - 5th bit: The cblk is part of the Schur complement if set
 *   - 6th bit: The cblk is a local copy of a remote fanin to be received.
 */
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define CBLK_FANIN      (1 << 0)
#define CBLK_LAYOUT_2D  (1 << 1)
#define CBLK_TASKS_2D   (1 << 2)
#define CBLK_COMPRESSED (1 << 3)
#define CBLK_IN_SCHUR   (1 << 4)
#define CBLK_IN_LAST    (1 << 5)
#define CBLK_RECV       (1 << 6)
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *@}
 */

/*
 * The type and structure definitions.
 */
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define COMP_1D                     0
#define DIAG                        1
#define E1                          2
#define E2                          3
#define DRUNK                       4
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 * @brief Tags used in MPI communications.
 */
typedef enum solve_step_e {
    PastixSolveForward,
    PastixSolveBackward,
    PastixFacto,
} solve_step_t;

/**
 * @brief Arguments for the solve.
 */
typedef struct args_solve_s
{
    pastix_scheduler_t  sched;
    solve_step_t        solve_step;
    pastix_solv_mode_t  mode;
    pastix_side_t       side;
    pastix_uplo_t       uplo;
    pastix_trans_t      trans;
    pastix_diag_t       diag;
} args_solve_t;

/**
 * @brief Solver recv block structure.
 */
typedef struct solver_blok_recv_s {
    pastix_int_t frownum; /**< First row contributed in the block (lrownum+1 of the original block if none) */
    pastix_int_t lrownum; /**< Last row contributed in the block (frownum-1 of the original block if none)  */
} solver_blok_recv_t;

/**
 * @brief Solver recv column block structure.
 *
 * The structure works as a list of cblk for each remote owner to compute the exact contributions.
 */
typedef struct solver_cblk_recv_s {
    struct solver_cblk_recv_s *next;
    pastix_int_t ownerid;
    pastix_int_t fcolnum; /**< First column contributed in the block (lcolnum+1 of the original cblk if none) */
    pastix_int_t lcolnum; /**< Last column contributed in the block (fcolnum-1 of the original cblk if none)  */
    solver_blok_recv_t bloktab[1]; /**< Array of reception blocks */
} solver_cblk_recv_t;

/**
 * @brief The task structure for the numerical factorization
 */
typedef struct task_s {
    pastix_int_t          taskid;  /**< COMP_1D DIAG E1 E2                                        */
    pastix_int_t          prionum; /**< Priority value for the factorization                      */
    pastix_int_t          cblknum; /**< Attached column block                                     */
    pastix_int_t          bloknum; /**< Attached block                                            */
    pastix_int_t volatile ctrbcnt; /**< Total number of contributions                             */
#if defined(PASTIX_DYNSCHED)
    int                   threadid;/**< Index of the bubble which contains the task               */
#endif
} Task;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define GPUID_UNDEFINED -2 /**< GPU still undefined       */
#define GPUID_NONE      -1 /**< Block not computed on GPU */
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 * @brief Solver block structure.
 */
typedef struct solver_blok_s {
    void        *handler[2]; /**< Runtime data handler                     */
    pastix_int_t lcblknm;    /**< Local column block                       */
    pastix_int_t fcblknm;    /**< Facing column block                      */
    pastix_int_t gbloknm;    /**< Index in global bloktab (UNUSED)         */
    pastix_int_t gfaninnm;   /**< Global fanin block index, -1 if cblk not (FANIN | RECV) */
    pastix_int_t frownum;    /**< First row index                          */
    pastix_int_t lrownum;    /**< Last row index (inclusive)               */
    pastix_int_t coefind;    /**< Index in coeftab                         */
    pastix_int_t browind;    /**< Index in browtab                         */
    int8_t       inlast;     /**< Index of the block among last separator (2), coupling with last separator (1) or other blocks (0) */
    int          iluklvl;    /**< The block ILU(k) level */

    /* LR structures */
    pastix_lrblock_t *LRblock[2]; /**< Store the blok (L/U) in LR format. Allocated for the cblk. */
} SolverBlok;

/**
 * @brief Solver column block structure.
 */
typedef struct solver_cblk_s  {
    pastix_atomic_lock_t lock;        /**< Lock to protect computation on the cblk         */
    volatile int32_t     ctrbcnt;     /**< Number of contribution to receive               */
    int8_t               cblktype;    /**< Type of cblk                                    */
    int8_t               partitioned; /**< Bitmask to know if one of the side has been partitioned or not */
    pastix_int_t         fcolnum;     /**< First column index (Global numbering)           */
    pastix_int_t         lcolnum;     /**< Last column index (Global numbering, inclusive) */
    SolverBlok          *fblokptr;    /**< First block in column (diagonal)                */
    pastix_int_t         stride;      /**< Column block stride                             */
    pastix_int_t         lcolidx;     /**< First column index (Local numbering), used for the rhs vectors      */
    pastix_int_t         brownum;     /**< First block in row facing the diagonal block in browtab, 0-based    */
    pastix_int_t         brown2d;     /**< First 2D-block in row facing the diagonal block in browtab, 0-based */
    pastix_int_t         sndeidx;     /**< Global index of the original supernode the cblk belongs to          */
    pastix_int_t         gcblknum;    /**< Global column block index                                           */
    pastix_int_t         bcscnum;     /**< Index in the bcsctab if local cblk, -1 otherwise (FANIN | RECV)     */
    pastix_int_t         gfaninnum;   /**< Global fanin column block index, -1 if not (FANIN | RECV)           */
    void                *lcoeftab;    /**< Coefficients access vector, lower part  */
    void                *ucoeftab;    /**< Coefficients access vector, upper part  */
    void                *handler[2];  /**< Runtime data handler                    */
    pastix_int_t         selevtx;     /**< Index to identify selected cblk for which intra-separator contributions are not compressed */
    int                  ownerid;     /**< Rank of the owner                       */
    int                  threadid;    /**< Rank of the accessing thread            */
    pastix_int_t         priority;    /**< Priority of the cblk                    */
} SolverCblk;


#ifndef DOXYGEN_SHOULD_SKIP_THIS
struct parsec_sparse_matrix_desc_s;
typedef struct parsec_sparse_matrix_desc_s parsec_sparse_matrix_desc_t;

struct starpu_sparse_matrix_desc_s;
typedef struct starpu_sparse_matrix_desc_s starpu_sparse_matrix_desc_t;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 * @brief Solver column block structure.
 *
 * This structure stores all the numerical information about the factorization,
 * as well as the structure of the problem. Only local information to each
 * process is stored in this structure.
 *
 */
struct solver_matrix_s {
    int restore; /**< Flag to indicate if it is require to restore data with
                      solverBackupRestore: 0: No need, 1:After solve,
                      2:After Factorization */
    pastix_int_t            baseval;       /**< Base value for numberings                                    */
    pastix_int_t            nodenbr;       /**< Number of local columns (after dof extension)                */
    pastix_int_t            coefnbr;       /**< Number of locally stored values (after dof extension)        */
    pastix_int_t            gcblknbr;      /**< Global number of column blocks                               */
    pastix_int_t            cblknbr;       /**< Local number of column blocks                                */
    pastix_int_t            gfanincblknbr; /**< Global number of fanin column blocks                         */
    pastix_int_t            faninnbr;      /**< Local number of fanin cblk (included in cblknbr)             */
    pastix_int_t            fanincnt;      /**< Number of sends to realize                                   */
    pastix_int_t            maxrecv;       /**< Maximum blok size for a cblk_recv                            */
    pastix_int_t            recvnbr;       /**< Local number of recv cblk (included in cblknbr)              */
    pastix_int_t            recvcnt;       /**< Number of receptions to realize                              */
    pastix_int_t            cblkmax1d;     /**< Rank of the last cblk not beeing enabled for 2D computations */
    pastix_int_t            cblkmin2d;     /**< Rank of the first cblk beeing enabled for 2D computations    */
    pastix_int_t            cblkmaxblk;    /**< Maximum number of blocks per cblk                            */
    pastix_int_t            cblkschur;     /**< Index of the first local cblk in Schur                       */
    pastix_int_t            nb2dcblk;      /**< Number of 2D cblks                                           */
    pastix_int_t            nb2dblok;      /**< Number of 2D blocks                                          */
    pastix_int_t            bloknbr;       /**< Number of blocks                                             */
    pastix_int_t            gbloknbr;      /**< Global number of blocks                                      */
    pastix_int_t            gfaninbloknbr; /**< Global number of fanin blocks                                */
    pastix_int_t            brownbr;       /**< Size of the browtab array                                    */
    SolverCblk   * restrict cblktab;       /**< Array of solver column blocks [+1]                           */
    SolverBlok   * restrict bloktab;       /**< Array of solver blocks        [+1]                           */
    pastix_int_t * restrict browtab;       /**< Array of blocks                                              */
    pastix_coeftype_t       flttype;       /**< valtab datatype: PastixFloat, PastixDouble, PastixComplex32 or PastixComplex64 */
    int                     globalalloc;   /**< Boolean for global allocation of coeftab  */

    pastix_int_t           *gcbl2loc;      /**< Array of local cblknum corresponding to gcblknum */

    pastix_lr_t             lowrank;       /**< Low-rank parameters                       */
    pastix_factotype_t      factotype;     /**< General or symmetric factorization?       */
    double                  diagthreshold; /**< Diagonal threshold for pivoting           */
    volatile int32_t        nbpivots;      /**< Number of pivots during the factorization */

#if defined(PASTIX_WITH_PARSEC)
    parsec_sparse_matrix_desc_t *parsec_desc;
#endif
#if defined(PASTIX_WITH_STARPU)
    starpu_sparse_matrix_desc_t *starpu_desc;
#endif

    pastix_int_t              offdmax;              /*+ Maximum size of the off-diagonal blocks for hetrf/sytrf temporary buffers +*/
    pastix_int_t              gemmmax;              /*+ Maximum size of the GEMM update for 1d GEMM computations                  +*/
    pastix_int_t              blokmax;              /*+ Maximum size of 2D blocks                 +*/
    pastix_int_t              colmax;               /*+ Maximum column width in solvmtx           +*/

    int                       clustnum;             /*+ current processor number                  +*/
    int                       clustnbr;             /*+ number of processors                      +*/
    pastix_int_t              procnbr;              /*+ Number of physical processor used         +*/
    pastix_int_t              thrdnbr;              /*+ Number of local computation threads       +*/
    pastix_int_t              bublnbr;              /*+ Number of local computation threads       +*/
    /* BubbleTree   * restrict   btree;                /\*+ Bubbles tree                              +*\/ */

    Task         * restrict   tasktab;              /*+ Task access vector                        +*/
    pastix_int_t              tasknbr;              /*+ Number of Tasks in 1d mode                +*/
    pastix_int_t              tasknbr_1dp;          /*+ Number of Tasks in 1d plus mode           +*/
    pastix_int_t **           ttsktab;              /*+ Task access vector by thread              +*/
    pastix_int_t *            ttsknbr;              /*+ Number of tasks by thread                 +*/
    pastix_queue_t **         computeQueue;         /*+ Queue of task to compute by thread        +*/

    pastix_int_t             *selevtx;              /*+ Array to identify which cblk are pre-selected +*/

    MPI_Request              *reqtab;               /**< Array of requests for MPI asynchronous messages      */
    pastix_int_t             *reqidx;               /**< Array of local cblknum index corresponding to reqtab */
    pastix_int_t              reqnbr;               /**< Length of the reqtab/reqidx arrays                   */
    volatile int32_t          reqnum;               /**< Current amount of active requests (packed)           */
    pastix_atomic_lock_t      reqlock;              /**< Lock to access the request arrays                    */
    void                     *rcoeftab;             /**< Reception buffer for the communication               */

    size_t                   *com_vector;           /**< Matrix of communications between nodes.              */

    PASTIX_Comm               solv_comm;            /*+ Copy of the pastix_data->inter_node_comm       +*/
};

/**
 *******************************************************************************
 *
 * @brief Computes the current solve step.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specify whether the off-diagonal blocks appear on the left or right
 *          in the equation. It has to be either PastixLeft or PastixRight.
 *
 * @param[in] uplo
 *          Specify whether the off-diagonal blocks are upper or lower
 *          triangular. It has to be either PastixUpper or PastixLower.
 *
 * @param[in] trans
 *          Specify the transposition used for the off-diagonal blocks. It has
 *          to be either PastixTrans or PastixConjTrans.
 *
 *******************************************************************************
 *
 * @return PastixSolveForward or PastixSolveBackward
 *
 *******************************************************************************/
static inline solve_step_t
compute_solve_step( pastix_side_t  side,
                    pastix_uplo_t  uplo,
                    pastix_trans_t trans )
{
    if ( ( (side == PastixLeft)  && (uplo == PastixUpper) && (trans == PastixNoTrans) ) ||
         ( (side == PastixLeft)  && (uplo == PastixLower) && (trans != PastixNoTrans) ) ||
         ( (side == PastixRight) && (uplo == PastixUpper) && (trans != PastixNoTrans) ) ||
         ( (side == PastixRight) && (uplo == PastixLower) && (trans == PastixNoTrans) ) )
    {
        return PastixSolveBackward;
    }
    else {
        return PastixSolveForward;
    }
}

/**
 * @brief     Compute the number of columns in a column block.
 * @param[in] cblk
 *            The pointer to the column block.
 * @return    The number of columns in the cblk.
 */
static inline pastix_int_t
cblk_colnbr( const SolverCblk *cblk )
{
    return cblk->lcolnum - cblk->fcolnum + 1;
}

/**
 * @brief     Get the pointer to the data associated to the lower part of the cblk.
 * @param[in] cblk
 *            The pointer to the column block.
 * @return    lcoeftab if the cblk is not compressed.
 *            fblokptr->LRblock[0] if the cblk is compressed.
 */
static inline void *
cblk_getdataL( const SolverCblk *cblk ) {
    return (cblk->cblktype & CBLK_COMPRESSED) ? cblk->fblokptr->LRblock[0] : cblk->lcoeftab;
}

/**
 * @brief     Get the pointer to the data associated to the upper part of the cblk.
 * @param[in] cblk
 *            The pointer to the column block.
 * @return    lcoeftab if the cblk is not compressed.
 *            fblokptr->LRblock[1] if the cblk is compressed.
 */
static inline void *
cblk_getdataU( const SolverCblk *cblk ) {
    return (cblk->cblktype & CBLK_COMPRESSED) ? cblk->fblokptr->LRblock[1] : cblk->ucoeftab;
}

/**
 * @brief     Get the pointer to the data associated to the side part of the cblk.
 * @param[in] cblk
 *            The pointer to the column block.
 * @param[in] side
 *            side of the data needed.
 *            PastixLCoef for the lower part PastixUCoef for the upper part.
 * @return    [lu]coeftab if the cblk is not compressed.
 *            fblokptr->LRblock[side] if the cblk is compressed.
 */
static inline void *
cblk_getdata( const SolverCblk *cblk, pastix_coefside_t side ) {
    return cblk->cblktype & CBLK_COMPRESSED
               ? cblk->fblokptr->LRblock[side]
           : side == PastixUCoef ? cblk->ucoeftab
                                 : cblk->lcoeftab;
}

/**
 * @brief     Compute the number of blocks in a column block.
 * @param[in] cblk
 *            The pointer to the column block.
 * @return    The number of blocks in the cblk including the diagonal block.
 */
static inline pastix_int_t
cblk_bloknbr( const SolverCblk *cblk )
{
    return (cblk+1)->fblokptr - cblk->fblokptr + 1;
}

/**
 * @brief     Compute the number of rows of a block.
 * @param[in] blok
 *            The pointer to the block.
 * @return    The number of rows in the block.
 */
static inline pastix_int_t
blok_rownbr( const SolverBlok *blok )
{
    return blok->lrownum - blok->frownum + 1;
}

/**
 * @brief     Compute the number of rows of a contiguous block in front of the same cblk.
 * @param[in] blok
 *            The pointer to the block.
 * @return    The number of rows in the block.
 */
static inline pastix_int_t
blok_rownbr_ext( const SolverBlok *blok )
{
    pastix_int_t rownbr = blok_rownbr( blok );

    while( (blok[0].lcblknm == blok[1].lcblknm) &&
           (blok[0].fcblknm == blok[1].fcblknm) )
    {
        blok++;
        rownbr += blok_rownbr( blok );
    }
    return rownbr;
}

/**
 * @brief Return if a block is preselected as either part of the projection, or
 *        as a sub-diagonal block.
 * @param[in] cblk
 *            The pointer to the cblk to which belong the block.
 * @param[in] blok
 *            The pointer to the block.
 * @param[in] fcbk
 *            The pointer to the facing cblk of the block.
 * @return True is the block is preselected, false otherwise.
 */
static inline int
blok_is_preselected( const SolverCblk *cblk,
                     const SolverBlok *blok,
                     const SolverCblk *fcbk )
{
    int is_preselected = ( fcbk->selevtx );
    int is_firstoffd   = ( blok == (cblk->fblokptr + 1) );

    return ( fcbk->sndeidx == cblk->sndeidx ) && ( is_preselected | is_firstoffd );
}

/**
 * @brief     Compute the number of rows of a column block.
 * @param[in] cblk
 *            The pointer to the column block.
 * @return    The number of rows in the column block.
 */
static inline pastix_int_t
cblk_rownbr( const SolverCblk *cblk )
{
    pastix_int_t rownbr = 0;
    SolverBlok * blok;
    for (blok = cblk->fblokptr; blok < cblk[1].fblokptr; blok++) {
        rownbr += blok_rownbr(blok);
    }
    return rownbr;
}

/**
 * @brief    Task stealing method.
 *
 * @param[inout] solvmtx
 *            The pointer to the solverMatrix.
 * @param[in] rank
 *            Rank of the computeQueue.
 * @param[in] nbthreads
 *            Total amount of threads in the node.
 * @return    The concerned cblk if it exists. -1 otherwhise.
 */
static inline pastix_int_t
stealQueue( SolverMatrix *solvmtx,
            int           rank,
            int           nbthreads )
{
    int rk = (rank + 1)%nbthreads;
    pastix_queue_t *stoleQueue;
    pastix_int_t    cblknum = -1;
    while( rk != rank )
    {
        assert( solvmtx->computeQueue[ rk ] );
        stoleQueue = solvmtx->computeQueue[ rk ];
        if( (cblknum = pqueuePop(stoleQueue)) != -1 ){
            return cblknum;
        }
        rk = (rk + 1)%nbthreads;
    }
    return cblknum;
}

/**
 * @brief Check if a block is included inside another one.
 *
 * Indicate if a blok is included inside an other block.
 * i.e. indicate if the row range of the first block is included in the
 * one of the second.
 *
 * @param[in] blok  The block that is tested for inclusion.
 * @param[in] fblok The block that is suppose to include the first one.
 *
 * @retval true   if the first block is     included in the second one.
 * @retval false  if the first block is not included in the second one.
 */
static inline int
is_block_inside_fblock( const SolverBlok *blok,
                        const SolverBlok *fblok )
{
#  if defined(NAPA_SOPALIN)
    return (((blok->frownum >= fblok->frownum) &&
             (blok->lrownum <= fblok->lrownum)) ||
            ((blok->frownum <= fblok->frownum) &&
             (blok->lrownum >= fblok->lrownum)) ||
            ((blok->frownum <= fblok->frownum) &&
             (blok->lrownum >= fblok->frownum)) ||
            ((blok->frownum <= fblok->lrownum) &&
             (blok->lrownum >= fblok->lrownum)));
#  else
    return ((blok->frownum >= fblok->frownum) &&
            (blok->lrownum <= fblok->lrownum));
#  endif /* defined(NAPA_SOPALIN) */
}

void solverInit( SolverMatrix *solvmtx );
void solverExit( SolverMatrix *solvmtx );

int  solverMatrixGen( SolverMatrix          *solvmtx,
                      const symbol_matrix_t *symbmtx,
                      const pastix_order_t  *ordeptr,
                      const SimuCtrl        *simuctl,
                      const BlendCtrl       *ctrl,
                      PASTIX_Comm            comm,
                      isched_t              *isched );

int  solverMatrixGenSeq( SolverMatrix          *solvmtx,
                         const symbol_matrix_t *symbmtx,
                         const pastix_order_t  *ordeptr,
                         const SimuCtrl        *simuctl,
                         const BlendCtrl       *ctrl,
                         PASTIX_Comm            comm,
                         isched_t              *isched,
                         pastix_int_t           is_dbg );

int  solverLoad( SolverMatrix       *solvptr,
                 FILE               *stream );
int  solverSave( const SolverMatrix *solvptr,
                 FILE               *stream );

void          solverRealloc( SolverMatrix       *solvptr);
SolverMatrix *solverCopy   ( const SolverMatrix *solvptr,
                             pastix_coeftype_t   flttype );

int           solverCheck     ( const SolverMatrix *solvmtx );
int           solverDraw      ( const SolverMatrix *solvptr,
                                FILE               *stream,
                                int                 verbose,
                                const char         *directory );
void          solverPrintStats( const SolverMatrix *solvptr );

void solverRequestInit( solve_step_t  solve_step,
                        SolverMatrix *solvmtx );
void solverRequestExit( SolverMatrix *solvmtx );

void solverRecvInit( pastix_coefside_t  side,
                     SolverMatrix      *solvmtx,
                     pastix_coeftype_t  flttype );
void solverRecvExit( SolverMatrix      *solvmtx );

void solverRhsRecvInit( solve_step_t       solve_step,
                        SolverMatrix      *solvmtx,
                        pastix_coeftype_t  flttype,
                        pastix_rhs_t       rhsb   );
void solverRhsRecvExit( SolverMatrix      *solvmtx );

/*
 * Solver backup
 */
#ifndef DOXYGEN_SHOULD_SKIP_THIS
struct SolverBackup_s;
typedef struct SolverBackup_s SolverBackup_t;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

SolverBackup_t *solverBackupInit   ( const SolverMatrix *solvmtx                          );
int             solverBackupRestore(       SolverMatrix *solvmtx, const SolverBackup_t *b );
void            solverBackupExit   (                                    SolverBackup_t *b );

#endif /* _solver_h_ */
