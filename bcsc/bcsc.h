/**
 *
 * @file bcsc.h
 *
 * @copyright 2004-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Alycia Lisito
 * @date 2022-10-11
 *
 * @addtogroup bcsc
 * @{
 *   @brief Describe all the internals routines to manipulate the internal block csc.
 *
 *   These functions provide a set of subroutines to manipulate the permuted
 *   sparse matrix stored in block of columns following the partition.
 *
 **/
#ifndef _bcsc_h_
#define _bcsc_h_

typedef enum bcsc_tag_ {
    PastixTagCount,
    PastixTagIndexesA,
    PastixTagIndexesAt,
    PastixTagValuesA,
    PastixTagValuesAt,
} bcsc_tag_e;

/**
 * @brief Information about the amount of data exchanged.
 */
typedef struct bcsc_data_amount_s
{
    pastix_int_t  idx_A;  /**< Amount of indexes of A which will be exchanged.  */
    pastix_int_t  val_A;  /**< Amount of values of A which will be exchanged.   */
    pastix_int_t  idx_At; /**< Amount of indexes of At which will be exchanged. */
    pastix_int_t  val_At; /**< Amount of values of At which will be exchanged.  */
} bcsc_data_amount_t;

/**
 * @brief Informations of the data exchanged with other processors.
 */
typedef struct bcsc_proc_comm_s
{
    bcsc_data_amount_t nsends;       /**< Number of indexes and values to send to clustnum.                  */
    bcsc_data_amount_t nrecvs;       /**< Number of indexes and values to receive from clustnum.             */
    pastix_int_t      *indexes_A;    /**< Array of indexes of A to send to clustnum.                         */
                                     /* indexes_A[2*k]   = kth column index to send to proc clustnum.        */
                                     /* indexes_A[2*k+1] = kth row index to send to proc clustnum.           */
    void              *values_A;     /**< Array of values of A to send to clustnum.                          */
                                     /* values_A is sorted the same way as the indexes_A, for each indexes   */
                                     /* there are dofi*dofj values.                                          */
    pastix_int_t      *indexes_At;   /**< Array of indexes of At to send to clustnum.                        */
                                     /* indexes_At[2*k]   = kth column index to send to proc clustnum.       */
                                     /* indexes_At[2*k+1] = kth row index to send to proc clustnum.          */
    void              *values_At;    /**< Array of values of At to send to clustnum.                         */
                                     /* values_At is sorted the same way as the indexes_At, for each indexes */
                                     /* there are dofi*dofj values.                                          */
} bcsc_proc_comm_t;

/**
 * @brief Structure to manage communications with distributed spm.
 */
typedef struct bcsc_handle_comm_s
{
    pastix_int_t      clustnbr;     /**< Number of processes in the cluster.                         */
    pastix_int_t      clustnum;     /**< ID of the current process in the cluster.                   */
    PASTIX_Comm       comm;         /**< PaStiX MPI communicator used for the ordering step.         */
    pastix_coeftype_t flttype;      /**< valtab datatype: PastixFloat, PastixDouble, PastixComplex32 or PastixComplex64 */
    bcsc_proc_comm_t  data_comm[1]; /**< Array of size clustnbr.                                     */
                                    /* data_comm[c]: contains the data clustnum has to send to c     */
                                    /*               and the amount of data clustnum will receive    */
                                    /*               from c.                                         */
                                    /* data_comm[clustnum]: contains the data clustnum will recevied */
                                    /*                     from the other processors.                */
} bcsc_handle_comm_t;

/**
 * @brief Compressed colptr format for the bcsc
 */
typedef struct bcsc_cblk_s {
    pastix_int_t  colnbr;  /**< Number of columns in the block column.                                    */
    pastix_int_t  cblknum; /**< Index of the corresponding cblk in the local solver matrix                */
    pastix_int_t *coltab;  /**< Array of indexes of the start of each column in the row and value arrays. */
} bcsc_cblk_t;

/**
 * @brief Internal column block distributed CSC matrix.
 */
struct pastix_bcsc_s {
    pastix_int_t        gN;        /**< Global number of vertices                                                          */
    pastix_int_t        n;         /**< Local number of vertices                                                           */
    pastix_mtxtype_t    mtxtype;   /**< Matrix structure: PastixGeneral, PastixSymmetric or PastixHermitian.               */
    pastix_coeftype_t   flttype;   /**< valtab datatype: PastixFloat, PastixDouble, PastixComplex32 or PastixComplex64     */
    pastix_int_t        cscfnbr;   /**< Number of column blocks.                                                           */
    bcsc_cblk_t        *cscftab;   /**< Array of Block column structures of size cscfnbr. (<pastix_bcscFormat_t>)          */
    pastix_int_t       *rowtab;    /**< Array of rows in the matrix.                                                       */
    void               *Lvalues;   /**< Array of values of the matrix A                                                    */
    void               *Uvalues;   /**< Array of values of the matrix A^t                                                  */
    pastix_int_t       *col2cblk;  /**< Array which gives the repartition of the solvmtx columns into the block structure. */
                                   /*   If the column k belongs to a remote block then col2cblk[k] = -(owner_proc + 1).    */
                                   /*   If the column k belongs to a local block then col2cblk[k] = block_num.             */
    bcsc_handle_comm_t *bcsc_comm; /**< Structure which handles the MPI communication (= NULL if PASTIX_WITH_MPI=OFF).     */
};

double bcscInit( const spmatrix_t     *spm,
                 const pastix_order_t *ord,
                 const SolverMatrix   *solvmtx,
                 pastix_int_t          initAt,
                 pastix_bcsc_t        *bcsc );

void   bcscExit( pastix_bcsc_t *bcsc );

/**
 * @}
 *
 * @addtogroup bcsc_internal
 * @{
 */
pastix_int_t bcsc_init_coltab( const spmatrix_t     *spm,
                               const pastix_order_t *ord,
                               const SolverMatrix   *solvmtx,
                               pastix_bcsc_t        *bcsc );
void bcsc_restore_coltab( pastix_bcsc_t *bcsc );

#if defined(PASTIX_WITH_MPI)
void bcsc_handle_comm_init( const SolverMatrix *solvmtx,
                            bcsc_handle_comm_t *bcsc_comm );
void bcsc_handle_comm_exit( bcsc_handle_comm_t *bcsc_comm );
void bcsc_exchange_amount_of_data( bcsc_handle_comm_t *bcsc_comm );
void bcsc_exchange_indexes( bcsc_handle_comm_t   *bcsc_comm );
#endif

/**
 * @}
 */
#endif /* _bcsc_h_ */
