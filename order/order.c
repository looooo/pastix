/**
 *
 * @file order.c
 *
 *  PaStiX order routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * Contains basic functions to manipulate the order structure.
 *
 * @version 5.1.0
 * @author Francois Pellegrini
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#include "common.h"
#include "order.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_ordering
 *
 * orderAlloc - Allocate the order structure. The base value is set to 0 by
 * default.
 *
 *******************************************************************************
 *
 * @param[in,out] ordeptr
 *          The data structure is set to 0 and then initialize.
 *          Need to call orderExit to release the memory first if required to
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
 * @return
 *          \retval PASTIX_SUCCESS on successful exit
 *          \retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *          \retval PASTIX_ERR_OUTOFMEMORY if one allocation failed.
 *
 *******************************************************************************/
int
orderAlloc( Order * const ordeptr,
            pastix_int_t vertnbr,
            pastix_int_t cblknbr)
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

    memset(ordeptr, 0, sizeof(Order));

    ordeptr->vertnbr = vertnbr;
    ordeptr->cblknbr = cblknbr;

    if (vertnbr != 0) {
        MALLOC_INTERN(ordeptr->permtab, vertnbr, pastix_int_t);
        MALLOC_INTERN(ordeptr->peritab, vertnbr, pastix_int_t);
    }

     if (cblknbr != 0) {
        MALLOC_INTERN(ordeptr->rangtab, cblknbr+1, pastix_int_t);
        MALLOC_INTERN(ordeptr->treetab, cblknbr,   pastix_int_t);
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_ordering
 *
 * orderInit - Initialize the order structure with the given values. The base
 * value is set to 0 by default. This is useful to give a personal ordering to
 * the pastix_task_order() function.
 *
 *******************************************************************************
 *
 * @param[in,out] ordeptr
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
 *              - each treetab[i] > i
 *              - all roots of the tree must have -1 as father
 *          If NULL, the rangtab field is not initialized.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PASTIX_SUCCESS on successful exit
 *          \retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *          \retval PASTIX_ERR_OUTOFMEMORY if one allocation failed.
 *
 *******************************************************************************/
int
orderInit ( Order * const ordeptr,
            pastix_int_t baseval,
            pastix_int_t vertnbr,
            pastix_int_t cblknbr,
            pastix_int_t *permtab,
            pastix_int_t *peritab,
            pastix_int_t *rangtab,
            pastix_int_t *treetab )
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

    memset(ordeptr, 0, sizeof(Order));

    ordeptr->baseval = baseval;
    ordeptr->vertnbr = vertnbr;
    ordeptr->cblknbr = cblknbr;

    if ( permtab ) {
        ordeptr->permtab = permtab;
    }
    if ( peritab ) {
        ordeptr->peritab = peritab;
    }
    if ( rangtab ) {
        ordeptr->rangtab = rangtab;
    }
    if ( treetab ) {
        ordeptr->treetab = treetab;
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_ordering
 *
 * orderExit - Free the arrays initialized in the order structure.
 *
 *******************************************************************************
 *
 * @param[in,out] ordeptr
 *          The data structure to clean. All arrays of the structure are freed
 *          and the structure is set to 0.
 *
 *******************************************************************************/
void
orderExit (Order * const ordeptr)
{
    /* Parameter check */
    if ( ordeptr == NULL ) {
        return;
    }

    if (ordeptr->permtab != NULL)
        memFree_null (ordeptr->permtab);
    if (ordeptr->peritab != NULL)
        memFree_null (ordeptr->peritab);
    if (ordeptr->rangtab != NULL)
        memFree_null (ordeptr->rangtab);
    if (ordeptr->treetab != NULL)
        memFree_null (ordeptr->treetab);

    memset(ordeptr, 0, sizeof(Order) );
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_ordering
 *
 * orderBase - This routine sets the base of the given ordering structure to the
 * given base value.
 *
 *******************************************************************************
 *
 * @param[in,out] ordeptr
 *          The ordering to rebase.
 *
 * @param[in] baseval
 *          The base value to be used (needs to be 0 or 1).
 *
 *******************************************************************************/
void
orderBase (Order * const ordeptr,
	   pastix_int_t  baseval)
{
    pastix_int_t baseadj;                    /* Base adjust */
    pastix_int_t cblknum;
    pastix_int_t vertnum;

    /* Parameter checks */
    if ( ordeptr == NULL ) {
        errorPrint("orderBase: ordeptr pointer is NULL");
        return;
    }
    if ( (baseval != 0) &&
         (baseval != 1) )
    {
        errorPrint("orderBase: baseval is incorrect, must be 0 or 1");
        return;
    }

    baseadj = baseval - ordeptr->baseval; /* Set base adjust     */
    if (baseadj == 0)                     /* If base already set */
	return;

    if (ordeptr->permtab != NULL) {
	for (vertnum = 0; vertnum < ordeptr->vertnbr; vertnum ++)
	    ordeptr->permtab[vertnum] += baseadj;
    }
    if (ordeptr->peritab != NULL) {
	for (vertnum = 0; vertnum < ordeptr->vertnbr; vertnum ++)
	    ordeptr->peritab[vertnum] += baseadj;
    }

    if (ordeptr->rangtab != NULL) {
        for (cblknum = 0; cblknum <= ordeptr->cblknbr; cblknum ++)
            ordeptr->rangtab[cblknum] += baseadj;
    }
    if (ordeptr->treetab != NULL) {
        for (cblknum = 0; cblknum < ordeptr->cblknbr; cblknum ++)
            ordeptr->treetab[cblknum] += baseadj;
    }

    ordeptr->baseval = baseval;
}
