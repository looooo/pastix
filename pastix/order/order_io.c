/**
 *
 * @file order_io.c
 *
 *  PaStiX order routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * Contains functions to read/write the order structure from/to the disk.
 *
 * @version 5.1.0
 * @author Francois Pellegrini
 * @date 2013-06-24
 *
 **/
#include "common.h"
#include "order.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_ordering_internal
 *
 * orderLoad - This routine reads the given ordering structure from the given
 * stream.
 *
 *******************************************************************************
 *
 * @param[in,out] ordeptr
 *          The ordering structure to fill in.
 *
 * @param[in] stream
 *          The stream where to read the informations.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PASTIX_SUCESS on successful exit.
 *          \retval PASTIX_ERR_FILE if a problem occurs during the read.
 *
 *******************************************************************************/
int
orderLoad (Order * const ordeptr,
           FILE  * const stream)
{
    pastix_int_t  versval;                      /* Version number */
    pastix_int_t  cblknbr;
    pastix_int_t  cblknum;
    pastix_int_t  vertnbr;
    pastix_int_t  vertnum;
    pastix_int_t  vertnnd;
    pastix_int_t *permtax;
    pastix_int_t *peritax;
    int           i;

    if ((intLoad (stream, &versval) +
         intLoad (stream, &cblknbr) +
         intLoad (stream, &vertnbr) != 3) ||
        (versval != 0)                    ||
        (cblknbr > vertnbr)) {
        errorPrint ("orderLoad: bad input (1)");
        return PASTIX_ERR_FILE;
    }

    orderInit(ordeptr, vertnbr, cblknbr);
    ordeptr->vertnbr = vertnbr;
    ordeptr->cblknbr = cblknbr;

    for (cblknum = 0, i = 1; (i == 1) && (cblknum <= cblknbr); cblknum ++) /* Read column-block data */
        i = intLoad (stream, &ordeptr->rangtab[cblknum]);

    for (vertnum = 0; (i == 1) && (vertnum < vertnbr); vertnum ++) /* Read direct permutation */
        i = intLoad (stream, &ordeptr->permtab[vertnum]);

    if (i != 1) {
        errorPrint ("orderLoad: bad input (2)");
        orderExit  (ordeptr);
        return PASTIX_ERR_FILE;
    }

    permtax = ordeptr->permtab - ordeptr->rangtab[0];
    peritax = ordeptr->peritab - ordeptr->rangtab[0];
    for (vertnum = ordeptr->rangtab[0], vertnnd = vertnum + vertnbr; /* Build inverse permutation */
         vertnum < vertnnd; vertnum ++)
        peritax[permtax[vertnum]] = vertnum;

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_ordering_internal
 *
 * orderSave - This routine writes the given ordering structure to the given
 * stream.
 *
 *******************************************************************************
 *
 * @param[in,out] ordeptr
 *          The ordering structure to dump to disk.
 *
 * @param[in] stream
 *          The stream where to write the ordering.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PASTIX_SUCESS on successful exit.
 *          \retval PASTIX_ERR_BADPARAMETER if the ordeptr structure is incorrect.
 *          \retval PASTIX_ERR_FILE if a problem occurs during the write.
 *
 *******************************************************************************/
int
orderSave (const Order * const ordeptr,
           FILE * const        stream)
{
    pastix_int_t vertnbr;
    pastix_int_t vertnum;
    pastix_int_t cblknum;
    int          o;

    if (ordeptr->rangtab == NULL) {
        errorPrint ("orderSave: cannot save ordering without column block data");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (ordeptr->permtab == NULL) {
        errorPrint ("orderSave: cannot save ordering without direct permutation data");
        return PASTIX_ERR_BADPARAMETER;
    }

    vertnbr = ordeptr->rangtab[ordeptr->cblknbr] -  /* Get number of nodes */
        ordeptr->rangtab[0];

    if (fprintf (stream, "0\n%ld\t%ld\n",
                 (long) ordeptr->cblknbr,
                 (long) vertnbr) == EOF) {
        errorPrint ("orderSave: bad output (1)");
        return PASTIX_ERR_FILE;
    }

    for (cblknum = 0, o = 1; (o == 1) && (cblknum < ordeptr->cblknbr); cblknum ++) { /* Save column-block range array */
        o = intSave (stream, ordeptr->rangtab[cblknum]);
        putc (((cblknum & 7) == 7) ? '\n' : '\t', stream);
    }
    o = intSave (stream, ordeptr->rangtab[cblknum]);
    putc ('\n', stream);

    for (vertnum = 0; (o == 1) && (vertnum < (vertnbr - 1)); vertnum ++) { /* Save direct permutation */
        o = intSave (stream, ordeptr->permtab[vertnum]);
        putc (((vertnum & 7) == 7) ? '\n' : '\t', stream);
    }
    o = intSave (stream, ordeptr->permtab[vertnum]);
    putc ('\n', stream);

    if (o != 1) {
        errorPrint ("orderSave: bad output (2)");
        return PASTIX_ERR_FILE;
    }

    return PASTIX_SUCESS;
}
