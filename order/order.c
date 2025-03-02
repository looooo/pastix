/**
 *
 * @file order.c
 *
 * PaStiX order structure routines
 *
 * @copyright 2004-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Francois Pellegrini
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Vincent Bridonneau
 * @date 2024-07-05
 *
 **/
#include "common.h"
#include "order/order_internal.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief Allocate the order structure.
 *
 * The base value is set to 0 by default.
 *
 *******************************************************************************
 *
 * @param[inout] ordeptr
 *          The data structure is set to 0 and then initialize.
 *          Need to call pastixOrderExit() to release the memory first if required to
 *          prevent memory leak.
 *
 * @param[in] vertnbr
 *          The number of nodes, this is the size of the internal permtab and
 *          peritab arrays.
 *          If vertnbr == 0, permtab and peritab are not allocated.
 *
 * @param[in] cblknbr
 *          The number of supernodes. The internal rangtab array is of size
 *          cblknbr+1, and treetab of size cblknbr.
 *          If cblknbr == 0, rangtab and treetab are not allocated.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect,
 * @retval PASTIX_ERR_OUTOFMEMORY if one allocation failed.
 *
 *******************************************************************************/
int
pastixOrderAlloc( pastix_order_t * const ordeptr,
                  pastix_int_t           vertnbr,
                  pastix_int_t           cblknbr )
{
    /* Parameter checks */
    if ( ordeptr == NULL ) {
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( vertnbr < 0 ) {
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( cblknbr < 0 ) {
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( cblknbr > vertnbr ) {
        return PASTIX_ERR_BADPARAMETER;
    }

    memset( ordeptr, 0, sizeof(pastix_order_t) );

    ordeptr->vertnbr = vertnbr;
    ordeptr->cblknbr = cblknbr;
    ordeptr->sndenbr = cblknbr;
    ordeptr->sndetab = NULL;

    if ( vertnbr > 0 ) {
        MALLOC_INTERN( ordeptr->permtab, vertnbr, pastix_int_t );
        MALLOC_INTERN( ordeptr->peritab, vertnbr, pastix_int_t );
    }

    if ( cblknbr > 0 ) {
        MALLOC_INTERN( ordeptr->rangtab, cblknbr+1, pastix_int_t );
        MALLOC_INTERN( ordeptr->treetab, cblknbr,   pastix_int_t );
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief Allocate the order structure for a given number of vertices with no
 * cblk, and id permutation.
 *
 * The base value is set to 0 by default.
 *
 *******************************************************************************
 *
 * @param[inout] ordeptr
 *          The data structure is set to 0 and then initialize.
 *          Need to call pastixOrderExit() to release the memory first if required to
 *          prevent memory leak.
 *
 * @param[in] vertnbr
 *          The number of nodes, this is the size of the internal permtab and
 *          peritab arrays.
 *          If vertnbr == 0, permtab and peritab are not allocated.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect,
 * @retval PASTIX_ERR_OUTOFMEMORY if one allocation failed.
 *
 *******************************************************************************/
int
pastixOrderAllocId( pastix_order_t * const ordeptr,
                    pastix_int_t           vertnbr )
{
    pastix_int_t *perm;
    pastix_int_t *invp;
    pastix_int_t  i;
    int rc;

    rc = pastixOrderAlloc( ordeptr, vertnbr, 1 );
    if ( rc != PASTIX_SUCCESS ) {
        return rc;
    }

    perm = ordeptr->permtab;
    invp = ordeptr->peritab;

    for(i=0; i<vertnbr; i++, perm++, invp++) {
        *perm = i;
        *invp = i;
    }

    ordeptr->rangtab[0] = 0;
    ordeptr->rangtab[1] = vertnbr;
    ordeptr->treetab[0] = -1;

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief Initialize the order structure with the given values.
 *
 * The base value is set to 0 by default. This is useful to give a personal
 * ordering to the pastix_task_order() function.
 *
 *******************************************************************************
 *
 * @param[inout] ordeptr
 *          The data structure is set to 0 and then initialize.
 *          Need to call orderExit to release the memory first if required to
 *          prevent memory leak.
 *
 * @param[in] baseval
 *          The base value used in the given arrays. Usually 0 for C, 1 for
 *          Fortran. Must be >= 0.
 *
 * @param[in] vertnbr
 *          The number of nodes, this is the size of the internal permtab and
 *          peritab arrays.
 *
 * @param[in] cblknbr
 *          The number of supernodes. The internal rangtab and treetab arrays
 *          are of size cblknbr+1.
 *
 * @param[in] permtab
 *          The permutation array which must be of size vertnbr, and based on
 *          baseval value.
 *          If NULL, the permtab field is not initialized.
 *
 * @param[in] peritab
 *          The inverse permutation array which must be of size vertnbr, and
 *          based on baseval value.
 *          If NULL, the peritab field is not initialized.
 *
 * @param[in] rangtab
 *          The rangtab array that describes the supernodes in the graph. This
 *          array must be of size cblknbr+1, and based on baseval value.
 *          If NULL, the rangtab field is not initialized.
 *
 * @param[in] treetab
 *          The treetab array that describes the elimination tree connecting the
 *          supernodes. This array must be defined as follow:
 *              - of size cblknbr;
 *              - based on baseval value;
 *              - each treetab[i] > i, unless i is a root and treetab[i] == -1
 *              - all roots of the tree must have -1 as father
 *          If NULL, the treetab field is not initialized.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 * @retval PASTIX_ERR_OUTOFMEMORY if one allocation failed.
 *
 *******************************************************************************/
int
pastixOrderInit( pastix_order_t * const ordeptr,
                 pastix_int_t           baseval,
                 pastix_int_t           vertnbr,
                 pastix_int_t           cblknbr,
                 pastix_int_t   * const permtab,
                 pastix_int_t   * const peritab,
                 pastix_int_t   * const rangtab,
                 pastix_int_t   * const treetab )
{
    /* Parameter checks */
    if ( ordeptr == NULL ) {
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( vertnbr < 0 ) {
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( cblknbr < 0 ) {
        return PASTIX_ERR_BADPARAMETER;
    }

    memset(ordeptr, 0, sizeof(pastix_order_t));

    ordeptr->baseval = baseval;
    ordeptr->vertnbr = vertnbr;
    ordeptr->cblknbr = cblknbr;
    ordeptr->sndenbr = cblknbr;
    ordeptr->sndetab = NULL;

    if ( permtab ) {
        ordeptr->permtab = permtab;
    }
    if ( peritab ) {
        ordeptr->peritab = peritab;
    }
    if ( rangtab ) {
        ordeptr->rangtab = rangtab;
        ordeptr->sndetab = malloc( (ordeptr->sndenbr+1) * sizeof(pastix_int_t) );
        memcpy( ordeptr->sndetab, ordeptr->rangtab, (ordeptr->sndenbr+1) * sizeof(pastix_int_t) );
    }
    if ( treetab ) {
        ordeptr->treetab = treetab;
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief Free the arrays initialized in the order structure.
 *
 *******************************************************************************
 *
 * @param[inout] ordeptr
 *          The data structure to clean. All arrays of the structure are freed
 *          and the structure is set to 0.
 *
 *******************************************************************************/
void
pastixOrderExit( pastix_order_t * const ordeptr )
{
    /* Parameter check */
    if ( ordeptr == NULL ) {
        return;
    }

    if (ordeptr->permtab != NULL) {
        memFree_null (ordeptr->permtab);
    }
    if (ordeptr->peritab != NULL) {
        memFree_null (ordeptr->peritab);
    }
    if (ordeptr->rangtab != NULL) {
        memFree_null (ordeptr->rangtab);
    }
    if (ordeptr->treetab != NULL) {
        memFree_null (ordeptr->treetab);
    }
    if (ordeptr->selevtx != NULL) {
        memFree_null (ordeptr->selevtx);
    }
    if (ordeptr->sndetab != NULL) {
        memFree_null (ordeptr->sndetab);
    }
    if (ordeptr->peritab_exp != NULL) {
        memFree_null (ordeptr->peritab_exp);
    }
    memset(ordeptr, 0, sizeof(pastix_order_t) );
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief This routine sets the base of the given ordering structure to the
 * given base value.
 *
 *******************************************************************************
 *
 * @param[inout] ordeptr
 *          The ordering to rebase.
 *
 * @param[in] baseval
 *          The base value to be used (needs to be 0 or 1).
 *
 *******************************************************************************/
void
pastixOrderBase( pastix_order_t * const ordeptr,
                 pastix_int_t           baseval )
{
    pastix_int_t baseadj;                    /* Base adjust */
    pastix_int_t cblknum;
    pastix_int_t vertnum;

    /* Parameter checks */
    if ( ordeptr == NULL ) {
        pastix_print_error( "pastixOrderBase: ordeptr pointer is NULL" );
        return;
    }
    if ( (baseval != 0) &&
         (baseval != 1) )
    {
        pastix_print_error( "pastixOrderBase: baseval is incorrect, must be 0 or 1" );
        return;
    }

    baseadj = baseval - ordeptr->baseval; /* Set base adjust     */
    if (baseadj == 0)                     /* If base already set */
	return;

    if (ordeptr->permtab != NULL) {
	for (vertnum = 0; vertnum < ordeptr->vertnbr; vertnum ++) {
	    ordeptr->permtab[vertnum] += baseadj;
        }
    }
    if (ordeptr->peritab != NULL) {
	for (vertnum = 0; vertnum < ordeptr->vertnbr; vertnum ++) {
	    ordeptr->peritab[vertnum] += baseadj;
        }
    }

    if (ordeptr->rangtab != NULL) {
        for (cblknum = 0; cblknum <= ordeptr->cblknbr; cblknum ++) {
            ordeptr->rangtab[cblknum] += baseadj;
        }
    }
    if (ordeptr->treetab != NULL) {
        for (cblknum = 0; cblknum < ordeptr->cblknbr; cblknum ++) {
            ordeptr->treetab[cblknum] += baseadj;
        }
    }
    if (ordeptr->sndetab != NULL) {
        pastix_int_t sndenum;
        for (sndenum = 0; sndenum <= ordeptr->sndenbr; sndenum ++) {
            ordeptr->sndetab[sndenum] += baseadj;
        }
    }

    ordeptr->baseval = baseval;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief This routine expand the permutation arrays and the rangtab when the
 * spm is using multiple dof per unknown.
 *
 *******************************************************************************
 *
 * @param[inout] ordeptr
 *          The ordering to expand. On entry, the order of the compressed
 *          matrix/graph. On exit, the ordering is 0-based and contains the
 *          permutation for the expanded matrix. peritab_ext remains untouched
 *          (=NULL) because peritab will already be expanded.
 *
 * @param[in] spm
 *          The sparse matrix structure providing dof information.
 *
 *******************************************************************************/
void pastixOrderExpand( pastix_order_t   *ordeptr,
                        const spmatrix_t *spm )
{
    pastix_int_t       *peritab;
    pastix_int_t        i, j, n;
    pastix_int_t        begin, end, sum_rang, sum_snde;
    pastix_int_t       *newperi;
    pastix_int_t       *rangtab;
    pastix_int_t       *sndetab;
    pastix_int_t        baseval;
    const pastix_int_t *dofs;

    baseval = spm->baseval;
    pastixOrderBase( ordeptr, 0 );

    n = spm->gNexp;

    /*
     * Initialize inverse permutation and rangtab
     */
    peritab = ordeptr->peritab;
    rangtab = ordeptr->rangtab;
    sndetab = ordeptr->sndetab;

    MALLOC_INTERN( ordeptr->peritab, n, pastix_int_t );
    newperi = ordeptr->peritab;

    dofs = spm->dofs;

    sum_rang = 0;
    sum_snde = 0;
    for (i=0; i<ordeptr->vertnbr; i++)
    {
        if ( spm->dof <= 0 ) {
            begin = dofs[ peritab[i] ] - baseval;
            end   = dofs[ peritab[i] + 1 ] - baseval;
        }
        else {
            begin = peritab[i] * spm->dof;
            end   = begin + spm->dof;
        }

        if ( i == rangtab[1] ) {
            rangtab[1] = rangtab[0] + sum_rang;
            rangtab++;
            sum_rang = 0;
        }
        if ( i == sndetab[1] ) {
            sndetab[1] = sndetab[0] + sum_snde;
            sndetab++;
            sum_snde = 0;
        }

        sum_rang += (end - begin);
        sum_snde += (end - begin);

        for ( j=begin; j<end; ++j, ++newperi ) {
            *newperi = j;
        }
    }
    rangtab[1] = rangtab[0] + sum_rang;
    sndetab[1] = sndetab[0] + sum_snde;

    ordeptr->vertnbr = n;
    memFree_null( peritab );

    /*
     * Update permtab
     */
    memFree_null( ordeptr->permtab );
    MALLOC_INTERN( ordeptr->permtab, n, pastix_int_t );
    for( i=0; i<n; i++ ) {
        j = ordeptr->peritab[i];
        ordeptr->permtab[j] = i;
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief This routine expand the peritab array for multi-dof matrices.
 *
 *******************************************************************************
 *
 * @param[inout] ordeptr
 *          TODO
 *
 * @param[in] spm
 *          The sparse matrix structure providing the dof information.
 *
 *******************************************************************************
 *
 * @return The expanded peritab array.
 *
 *******************************************************************************/
pastix_int_t *
orderGetExpandedPeritab( pastix_order_t   *ordeptr,
                         const spmatrix_t *spm )
{
    pastix_int_t       *peritab, *peritab_exp;
    pastix_int_t        i, j;
    pastix_int_t        begin, end;
    pastix_int_t        baseval_spm;
    pastix_int_t        baseval_ord;
    const pastix_int_t *dofs;

    assert( ordeptr != NULL );
    assert( spm     != NULL );
    assert( spm->gN == ordeptr->vertnbr );

    if ( spm->dof == 1 ) {
        return ordeptr->peritab;
    }

    if ( ordeptr->peritab_exp != NULL ) {
        return ordeptr->peritab_exp;
    }

    MALLOC_INTERN( peritab_exp, spm->gNexp, pastix_int_t );
    ordeptr->peritab_exp = peritab_exp;

    baseval_spm = spm->baseval;
    baseval_ord = ordeptr->baseval;

    /* Shift the dof ptr to take into account the baseval of peritab */
    dofs    = spm->dofs - baseval_ord;
    peritab = ordeptr->peritab;

    for ( i = 0; i < ordeptr->vertnbr; i++, peritab++ ) {
        if ( spm->dof <= 0 ) {
            begin = dofs[ peritab[0]   ] - baseval_spm;
            end   = dofs[ peritab[0]+1 ] - baseval_spm;
        }
        else {
            begin = (peritab[0] - baseval_ord) * spm->dof;
            end   = begin + spm->dof;
        }
        for ( j = begin; j < end; ++j, ++peritab_exp ) {
            *peritab_exp = j;
        }
    }

    return ordeptr->peritab_exp;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief This routine copy a given ordering in a new one.
 *
 * This function copies an order structure into another one. If all subpointers
 * are NULL, then they are all allocated and conatins the original ordesrc
 * values on exit. If one or more array pointers are not NULL, then, only those
 * are copied to the ordedst structure.
 *
 *******************************************************************************
 *
 * @param[inout] ordedst
 *          The destination ordering
 *
 * @param[in] ordesrc
 *          The source ordering
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
pastixOrderCopy( pastix_order_t       * const ordedst,
                 const pastix_order_t * const ordesrc )
{
    /* Parameter checks */
    if ( ordedst == NULL ) {
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( ordesrc == NULL ) {
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( ordesrc == ordedst ) {
        return PASTIX_ERR_BADPARAMETER;
    }

    ordedst->baseval = ordesrc->baseval;
    ordedst->vertnbr = ordesrc->vertnbr;
    ordedst->cblknbr = ordesrc->cblknbr;
    ordedst->sndenbr = ordesrc->sndenbr;

    if ( (ordedst->permtab == NULL) &&
         (ordedst->peritab == NULL) &&
         (ordedst->rangtab == NULL) &&
         (ordedst->treetab == NULL) )
    {
        pastixOrderAlloc( ordedst, ordedst->vertnbr, ordedst->cblknbr );
    }

    if ( (ordesrc->permtab != NULL) && (ordedst->permtab != NULL) )
    {
        memcpy( ordedst->permtab, ordesrc->permtab, ordesrc->vertnbr * sizeof(pastix_int_t) );
    }

    if ( (ordesrc->peritab != NULL) && (ordedst->peritab != NULL) )
    {
        memcpy( ordedst->peritab, ordesrc->peritab, ordesrc->vertnbr * sizeof(pastix_int_t) );
    }

    if ( (ordesrc->rangtab != NULL) && (ordedst->rangtab != NULL) )
    {
        memcpy( ordedst->rangtab, ordesrc->rangtab, (ordesrc->cblknbr+1) * sizeof(pastix_int_t) );
    }

    if ( (ordesrc->treetab != NULL) && (ordedst->treetab != NULL) )
    {
        memcpy( ordedst->treetab, ordesrc->treetab, ordesrc->cblknbr * sizeof(pastix_int_t) );
    }

    if ( (ordesrc->sndetab != NULL) && (ordedst->sndetab != NULL) )
    {
        memcpy( ordedst->sndetab, ordesrc->sndetab, (ordesrc->sndenbr+1) * sizeof(pastix_int_t) );
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief This routine returns the pointer to the internal order structure to
 * access permutation information.
 *
 * @warning The data returned by the routines must not be freed or modified by
 * the user, and are valid as long as no operation on the ordering is performed
 * (pastix_subtask_order(), pastix_subtask_reordering(), or pastixFinalize()).
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure of the problem
 *
 *******************************************************************************
 *
 * @return The pointer to the internal ordering structure with permutation
 *         information, or NULL if pastix_data is invalid.
 *
 *******************************************************************************/
pastix_order_t *
pastixOrderGet( const pastix_data_t * const pastix_data )
{
    /* Parameter checks */
    if ( pastix_data == NULL ) {
        return NULL;
    }
    return pastix_data->ordemesh;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief This routine broadcast the ordemesh structure from node root to all
 *        the other nodes.
 *
 *******************************************************************************
 *
 * @param[inout] ordemesh
 *          The ordering structure of the problem.
 *
 * @param[in] root
 *          The node that will broadcast its ordemesh.
 *
 * @param[in] pastix_comm
 *          The MPI communicator of the problem.
 *
 *******************************************************************************/
void
pastixOrderBcast( pastix_order_t *ordemesh,
                  int             root,
                  PASTIX_Comm     pastix_comm )
{
    pastix_int_t    vertnbr, cblknbr, sndenbr;
    int             clustnum;

    MPI_Comm_rank( pastix_comm, &clustnum );

    /* Free previous order structure */
    if ( clustnum != root ) {
        pastixOrderExit( ordemesh );
    }

    /* Copy ordemesh datas from node 0 */
    MPI_Bcast( ordemesh, sizeof(pastix_order_t), MPI_BYTE,
               root, pastix_comm );

    vertnbr = ordemesh->vertnbr;
    cblknbr = ordemesh->cblknbr;
    sndenbr = ordemesh->sndenbr;

    if ( ordemesh->permtab ) {
        if ( clustnum != root ) {
            MALLOC_INTERN( ordemesh->permtab, vertnbr, pastix_int_t );
        }
        MPI_Bcast( ordemesh->permtab, vertnbr, PASTIX_MPI_INT,
                   root, pastix_comm );
    }

    if ( ordemesh->peritab ) {
        if ( clustnum != root ) {
            MALLOC_INTERN( ordemesh->peritab, vertnbr, pastix_int_t );
        }
        MPI_Bcast( ordemesh->peritab, vertnbr, PASTIX_MPI_INT,
                   root, pastix_comm );
    }

    if ( ordemesh->rangtab ) {
        if ( clustnum != root ) {
            MALLOC_INTERN( ordemesh->rangtab, cblknbr+1, pastix_int_t );
        }
        MPI_Bcast( ordemesh->rangtab, cblknbr+1, PASTIX_MPI_INT,
                   root, pastix_comm );
    }

    if ( ordemesh->treetab ) {
        if ( clustnum != root ) {
            MALLOC_INTERN( ordemesh->treetab, cblknbr, pastix_int_t );
        }
        MPI_Bcast( ordemesh->treetab, cblknbr, PASTIX_MPI_INT,
                   root, pastix_comm );
    }

    if ( ordemesh->selevtx ) {
        if ( clustnum != root ) {
            MALLOC_INTERN( ordemesh->selevtx, cblknbr, int8_t );
        }
        MPI_Bcast( ordemesh->selevtx, cblknbr * sizeof(int8_t), MPI_BYTE,
                   root, pastix_comm );
    }

    if ( ordemesh->sndetab ) {
        if ( clustnum != root ) {
            MALLOC_INTERN( ordemesh->sndetab, sndenbr+1, pastix_int_t );
        }
        MPI_Bcast( ordemesh->sndetab, sndenbr+1, PASTIX_MPI_INT,
                   root, pastix_comm );
    }
}
