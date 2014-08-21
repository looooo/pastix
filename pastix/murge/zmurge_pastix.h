/**
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/
#ifndef ZMURGE_PASTIX_H
#define ZMURGE_PASTIX_H
#ifdef FORCE_NOMPI
#  include "nompi.h"
#else /* not FORCE_NOMPI */
#  include <mpi.h>
#endif /* not FORCE_NOMPI */
#include "zmurge.h"
/******************************************************************************
 * Title: PaStiX specific addons to Murge                                     *
 *                                                                            *
 * Functions declaration for Murge function defined only in PaStiX.           *
 *                                                                            *
 * Authors:                                                                   *
 *   Xavier Lacoste - xavier.lacoste@inria.fr                                 *
 *                                                                            *
 * More informations can be found in <murge.h> documentation.                 *
 *                                                                            *
 ******************************************************************************/


/******************************************************************************
 * Function: ZMURGE_Analyze                                                   *
 *                                                                            *
 * Perform matrix analyze.                                                    *
 *                                                                            *
 * Follows several steps :                                                    *
 *   - Compute a new ordering of the unknows                                  *
 *   - Compute the symbolic factorisation of the matrix                       *
 *   - Distribute column blocks and computation on processors                 *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *                                                                            *
 * Returns:                                                                   *
 *   MURGE_SUCCESS       - If function runned succesfuly.                     *
 *   MURGE_ERR_ORDER     - If function the graph is not built.                *
 *   MURGE_ERR_PARAMETER - If *murge_id* is not a valid ID.                   *
 *                                                                            *
 * Fortran interface:                                                         *
 >
 > SUBROUTINE ZMURGE_Analyze(ID, IERROR)
 >   INTS, INTENT(IN)  :: ID
 >   INTS, INTENT(OUT) :: IERROR
 > END SUBROUTINE ZMURGE_Analyze
 ******************************************************************************/

INTS ZMURGE_Analyze(INTS id);

/******************************************************************************
 * Function: ZMURGE_Factorize                                                  *
 *                                                                            *
 * Perform matrix factorization.                                              *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *                                                                            *
 * Returns:                                                                   *
 *   MURGE_SUCCESS       - If function runned succesfuly.                     *
 *   MURGE_ERR_ORDER     - If function the graph is not built.                *
 *   MURGE_ERR_PARAMETER - If *murge_id* is not a valid ID.                   *
 *                                                                            *
 * Fortran interface:                                                         *
 >
 > SUBROUTINE ZMURGE_Factorize(ID, IERROR)
 >   INTS, INTENT(IN)  :: ID
 >   INTS, INTENT(OUT) :: IERROR
 > END SUBROUTINE ZMURGE_Factorize
 ******************************************************************************/

INTS ZMURGE_Factorize(INTS id);

/******************************************************************************
 * Function: ZMURGE_SetOrdering                                                *
 *                                                                            *
 * Set a personal ordering to perform factorization.                          *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *   permutation - Permutation to perform factorization.                      *
 *                                                                            *
 * Fortran interface:                                                         *
 >
 > SUBROUTINE ZMURGE_SetOrdering(ID, PERMUTATION, IERROR)
 >   INTS, INTENT(IN)               :: ID
 >   INTS, DIMENSION(0), INTENT(IN) :: PERMUTATION
 >   INTS, INTENT(OUT)              :: IERROR
 > END SUBROUTINE ZMURGE_SetOrdering

 ******************************************************************************/
INTS ZMURGE_SetOrdering(INTS   id,
                       INTS * permutation);

/******************************************************************************
 * Function: ZMURGE_ProductSetLocalNodeNbr                                     *
 *                                                                            *
 * Set local node number if users only perform product and doesn't need       *
 * to perform first steps.                                                    *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *   n  - Number of local nodes.                                              *
 *                                                                            *
 Fortran intrerface:
 >
 > SUBROUTINE ZMURGE_SetLocalNodeNbr(ID, N, IERROR)
 >   INTS, INTENT(IN)               :: ID, N
 >   INTS, INTENT(OUT)              :: IERROR
 > END SUBROUTINE ZMURGE_SetLocalNodeNbr
 ******************************************************************************/
INTS ZMURGE_ProductSetLocalNodeNbr (INTS id, INTS n);

/******************************************************************************
 * Function: ZMURGE_ProductSetGlobalNodeNbr                                    *
 *                                                                            *
 * Set global node number if users only perform product and doesn't need      *
 * to perform first steps.                                                    *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *   N  - Number of global nodes.                                             *
 *                                                                            *
 Fortran intrerface:
 >
 > SUBROUTINE ZMURGE_SetGlobalNodeNbr(ID, N, IERROR)
 >   INTS, INTENT(IN)               :: ID, N
 >   INTS, INTENT(OUT)              :: IERROR
 > END SUBROUTINE ZMURGE_SetGlobalNodeNbr
 ******************************************************************************/
INTS ZMURGE_ProductSetGlobalNodeNbr (INTS id, INTS N);

/******************************************************************************
 * Function: ZMURGE_ProductSetLocalNodeList                                    *
 *                                                                            *
 * Set local node list if users only perform product and doesn't need         *
 * to perform first steps.                                                    *
 *                                                                            *
 * Parameters:                                                                *
 *   id  - Solver instance identification number.                             *
 *   l2g - Local to global node numbers.                                      *
 *                                                                            *
 Fortran interface:
 >
 > SUBROUTINE ZMURGE_ProductSetLocalNodeList(ID, L2G, IERROR)
 >   INTS, INTENT(IN)               :: ID
 >   INTS, DIMENSION(0), INTENT(IN) :: L2G
 >   INTS, INTENT(OUT)              :: IERROR
 > END SUBROUTINE ZMURGE_ProductSetLocalNodeList
 ******************************************************************************/
INTS ZMURGE_ProductSetLocalNodeList (INTS id, INTS * l2g);

/******************************************************************************
 * Function: ZMURGE_GetLocalProduct                                            *
 *                                                                            *
 * Perform the product A * X.                                                 *
 *                                                                            *
 * The vector must have been given trough <ZMURGE_SetLocalRHS> or              *
 * <ZMURGE_SetGlobalRHS>.                                                      *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *   x  - Array in which the local part of the product will be stored.        *
 *                                                                            *
 Fortran interface:
 >
 > SUBROUTINE ZMURGE_GetLocalProduct(ID, X, IERROR)
 >   INTS, INTENT(IN)               :: ID
 >   COEF, DIMENSION(0), INTENT(INOUT) :: X
 >   INTS, INTENT(OUT)              :: IERROR
 > END SUBROUTINE ZMURGE_GetLocalProduct
 ******************************************************************************/
INTS ZMURGE_GetLocalProduct (INTS id, COEF *x);

/******************************************************************************
 * Function: ZMURGE_GetGlobalProduct                                           *
 *                                                                            *
 * Perform the product A * X.                                                 *
 *                                                                            *
 * The vector must have been given trough <ZMURGE_SetLocalRHS> or              *
 * <ZMURGE_SetGlobalRHS>.                                                      *
 *                                                                            *
 * Parameters:                                                                *
 *   id   - Solver instance identification number.                            *
 *   x    - Array in which the product will be stored.                        *
 *   root - Rank of the process which will own the product at end of call,    *
 *          use -1 for all processes.                                         *
 * Returns:                                                                   *
 *   ZMURGE_ERR_ORDER  - If values have not been set.                          *
 *                                                                            *
 *                                                                            *
 Fortran interface:
 >
 > SUBROUTINE ZMURGE_GetGlobalProduct(ID, X, ROOT, IERROR)
 >   INTS, INTENT(IN)               :: ID, ROOT
 >   COEF, DIMENSION(0), INTENT(INOUT) :: X
 >   INTS, INTENT(OUT)              :: IERROR
 > END SUBROUTINE ZMURGE_GetGlobalProduct
 ******************************************************************************/
INTS ZMURGE_GetGlobalProduct (INTS id, COEF *x, INTS root);

/******************************************************************************
 * Function: ZMURGE_ForceNoFacto                                               *
 *                                                                            *
 * Prevent Zmurge from running factorisation even if matrix has changed.       *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 * Returns:                                                                   *
 *   ZMURGE_SUCCESS                                                            *
 *                                                                            *
 Fortran interface:
 >  SUBROUTINE ZMURGE_FORCENOFACTO(ID, IERROR)
 >    INTS,               INTENT(IN)  :: ID
 >    INTS,               INTENT(OUT) :: IERROR
 >  END SUBROUTINE ZMURGE_FORCENOFACTO
 ******************************************************************************/
INTS ZMURGE_ForceNoFacto(INTS id);

/******************************************************************************
 * Function: ZMURGE_SetLocalNodeList                                           *
 *                                                                            *
 * Deprecated, need to be checked                                             *
 *                                                                            *
 >  SUBROUTINE ZMURGE_SETLOCALNODELIST(ID, N, LIST, IERROR)
 >    INTS,               INTENT(IN)  :: ID, N
 >    INTS, DIMENSION(N), INTENT(IN)  :: LIST
 >    INTS,               INTENT(OUT) :: IERROR
 >  END SUBROUTINE ZMURGE_SETLOCALNODELIST
 ******************************************************************************/
INTS ZMURGE_SetLocalNodeList(INTS id, INTS n, INTS * list);

/******************************************************************************
 * Function: ZMURGE_AssemblySetSequence                                        *
 *                                                                            *
 * Create a sequence of entries to build a matrix and store it for being      *
 * reused.                                                                    *
 *                                                                            *
 * Parameters:                                                                *
 *   id      - Solver instance identification number.                         *
 *   coefnbr - Number of entries.                                             *
 *   ROWs    - List of rows in the sequence.                                  *
 *   COLs    - List of columns in the sequence.                               *
 *   op      - Operation to perform for coefficient which appear              *
 *             several tim (see <MURGE_ASSEMBLY_OP>).                         *
 *   op2     - Operation to perform when a coefficient is set by              *
 *             two different processors (see <MURGE_ASSEMBLY_OP>).            *
 *   mode    - Indicates if user ensure he will respect solvers distribution  *
 *             (see <MURGE_ASSEMBLY_MODE>).                                   *
 *   nodes   - 0 entries are entered value by value,                          *
 *             1 entries are entries node by node.                            *
 *   id_seq  - Sequence ID.                                                   *
 *                                                                            *
 * Returns:                                                                   *
 *   MURGE_SUCCESS       - If function runned successfully.                   *
 *   MURGE_ERR_ORDER     - If graph hasn't been built before.                 *
 *   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.                       *
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or          *
 *                         *op*, *mode*, *sym*, or *coefnbr* are not valid.   *
 >  SUBROUTINE ZMURGE_ASSEMBLYSETSEQUENCE (ID, NNZ, ROWs, COLs, OP, OP2, MODE, NODES, ID_SEQ, IERROR)
 >    INTS,                     INTENT(IN)  :: ID, OP, OP2, MODE, NODES
 >    INTL,                     INTENT(IN)  :: NNZ
 >    INTS, DIMENSION(NNZ),     INTENT(IN)  :: ROWs, COLs
 >    INTS,                     INTENT(OUT) :: ID_SEQ, IERROR
 >  END SUBROUTINE ZMURGE_ASSEMBLYSETSEQUENCE
 ******************************************************************************/
INTS ZMURGE_AssemblySetSequence (INTS id, INTL coefnbr, INTS * ROWs, INTS * COLs, INTS op, INTS op2, INTS mode, INTS nodes,INTS * id_seq);

/******************************************************************************
 * ZMURGE_AssemblySetSequence                                                  *
 *                                                                            *
 * Assembly the matrix using a stored sequence.                               *
 *                                                                            *
 * Parameters:                                                                *
 *   id      - Solver instance identification number.                         *
 *   id_seq  - Sequence ID.                                                   *
 *   values  - Values to insert in the CSC.                                   *
 *                                                                            *
 * Returns:                                                                   *
 *   MURGE_SUCCESS       - If function runned successfully.                   *
 *   MURGE_ERR_ORDER     - If graph hasn't been built before.                 *
 *   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.                       *
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or          *
 *                         *id_seq* or *values* are not valid.                *
 Fortran interface:
 >  SUBROUTINE ZMURGE_ASSEMBLYUSESEQUENCE(ID, ID_SEQ, VALUES, IERROR)
 >    INTS,               INTENT(IN)  :: ID, ID_SEQ 
 >    COEF, DIMENSION(0), INTENT(IN)  :: VALUES
 >    INTS,               INTENT(OUT) :: IERROR
 >  END SUBROUTINE ZMURGE_ASSEMBLYUSESEQUENCE
 ******************************************************************************/
INTS ZMURGE_AssemblyUseSequence(INTS id, INTS id_seq, COEF * values);

/******************************************************************************
 * Function: ZMURGE_AssemblyDeleteSequence                                     *
 *                                                                            *
 * Destroy an assembly sequence                                               *
 *                                                                            *
 *   id      - Solver instance identification number.                         *
 *   id_seq  - Sequence ID.                                                   *
 *                                                                            *
 * Returns:                                                                   *
 *   MURGE_SUCCESS       - If function runned successfully.                   *
 *   MURGE_ERR_ORDER     - If graph hasn't been built before.                 *
 *   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.                       *
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or          *
 *                         *id_seq* is not valid.                             *
 Fortran interface:
 >  SUBROUTINE ZMURGE_ASSEMBLYDELETESEQUENCE(ID, ID_SEQ, IERROR)
 >    INTS,               INTENT(IN)  :: ID, ID_SEQ 
 >    INTS,               INTENT(OUT) :: IERROR
 >  END SUBROUTINE ZMURGE_ASSEMBLYDELETESEQUENCE
 ******************************************************************************/
INTS ZMURGE_AssemblyDeleteSequence(INTS id, INTS id_seq);



INTS ZMURGE_GetCommRank(INTS id, int * rank);
INTS ZMURGE_GetCommSize(INTS id, int * size);
INTS ZMURGE_GetOptionINT(INTS id, INTS index, INTS * value);
INTS ZMURGE_GetComm(INTS id, MPI_Comm * comm);
/* INTS ZMURGE_GetLocalElementNbr(INTS id, INTS N, INTS globalElementNbr, INTS * localElementNbr, INTS mode, void * d); */
/* INTS ZMURGE_GetLocalElementList(INTS id, INTS * element_list); */
INTS ZMURGE_GraphSetEdge (INTS id, INTS ROW, INTS COL);
INTS ZMURGE_GraphSetBlockEdges(INTS id, INTS nROW, INTS *ROWlist, INTS nCOL, INTS *COLlist);

/*
 >  SUBROUTINE ZMURGE_SETDROPNODES(ID, NODENBR, DROPMASK, IERROR)
 >    INTS,               INTENT(IN) :: ID, NODENBR
 >    INTS, DIMENSION(0), INTENT(IN) :: DROPMASK
 >    INTS,               INTENT(OUT) :: IERROR
 >  END SUBROUTINE ZMURGE_SETDROPNODES
 */
INTS ZMURGE_SetDropNodes(INTS id, INTS nodenbr, INTS * dropmask);

/*
 >  SUBROUTINE ZMURGE_SETDROPCOLS(ID, NODENBR, DROPMASK, IERROR)
 >    INTS,               INTENT(IN) :: ID, NODENBR
 >    INTS, DIMENSION(0), INTENT(IN) :: DROPMASK
 >    INTS,               INTENT(OUT) :: IERROR
 >  END SUBROUTINE ZMURGE_SETDROPCOLS
 */
INTS ZMURGE_SetDropCols(INTS id, INTS nodenbr, INTS * dropmask);
/*
 >  SUBROUTINE ZMURGE_SETDROPROWS(ID, NODENBR, DROPMASK, IERROR)
 >    INTS,               INTENT(IN) :: ID, NODENBR
 >    INTS, DIMENSION(0), INTENT(IN) :: DROPMASK
 >    INTS,               INTENT(OUT) :: IERROR
 >  END SUBROUTINE ZMURGE_SETDROPROWS
 */
INTS ZMURGE_SetDropRows(INTS id, INTS nodenbr, INTS * dropmask);
/*
 >  SUBROUTINE ZMURGE_COLGETNONZEROSNBR(ID, ROW, NNZNBR, IERROR)
 >    INTS,               INTENT(IN) :: ID, ROW
 >    INTS,               INTENT(INOUT) :: NNZNBR
 >    INTS,               INTENT(OUT) :: IERROR
 >  END SUBROUTINE ZMURGE_COLGETNONZEROSNBR
 */
INTS ZMURGE_ColGetNonZerosNbr(INTS id, INTS COL, INTS * nnzNbr);
/*
 >  SUBROUTINE ZMURGE_COLGETNONZEROSIDX(ID, ROW, NNZIDX, IERROR)
 >    INTS,               INTENT(IN) :: ID, ROW
 >    INTS, DIMENSION(0), INTENT(INOUT) :: NNZIDX
 >    INTS,               INTENT(OUT) :: IERROR
 >  END SUBROUTINE ZMURGE_COLGETNONZEROSIDX
 */
INTS ZMURGE_ColGetNonZerosIdx(INTS id, INTS COL, INTS * indexes);

INTS ZMURGE_AssemblySetListOfBlockValues(INTS id, INTS nBlocks, INTS nROW, INTS *ROWlist, INTS nCOL, INTS *COLlist, COEF *values);


/* This comment must remain for generation purpose */
/*
 Murge's constants
 */
enum MURGE_PASTIX_ERR {
  MURGE_ERR_MPI      = 1024,
  MURGE_ERR_INTERNAL = 1025
};

enum MURGE_ELEMENT_DIST {
  MURGE_DUPLICATE_ELEMENTS,
  MURGE_DISTRIBUTE_ELEMENTS
};

enum PASTIX_MURGE_ASSEMBLY_MODE {
  MURGE_ASSEMBLY_DROPNONLOCAL = 3
};
#endif /* ZMURGE_PASTIX_H */
