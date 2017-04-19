/**
 *
 * @file order_io.c
 *
 * PaStiX order functions to read/write the order structure from/to the disk.
 *
 * @copyright 2004-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Francois Pellegrini
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#include "common.h"
#include "order.h"

/**
 *******************************************************************************
 *
 * @ingroup order_dev
 *
 * @brief Load an ordering from a file.
 *
 *******************************************************************************
 *
 * @param[inout] ordeptr
 *          The ordering structure to fill in.
 *
 * @param[in] stream
 *          The stream where to read the informations.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCESS on successful exit,
 * @retval PASTIX_ERR_FILE if a problem occurs during the read.
 *
 *******************************************************************************/
static inline int
ordering_load(Order * ordeptr,
              FILE  * stream)
{
    pastix_int_t  versval;
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
         intLoad (stream, &vertnbr) != 3)  ||
        ((versval != 0) && (versval != 1)) ||
        (cblknbr > vertnbr)) {
        errorPrint ("orderLoad: bad input (1)");
        return PASTIX_ERR_FILE;
    }

    orderAlloc(ordeptr, vertnbr, cblknbr);
    ordeptr->vertnbr = vertnbr;
    ordeptr->cblknbr = cblknbr;

    /* Read column-block data */
    for (cblknum = 0, i = 1; (i == 1) && (cblknum <= cblknbr); cblknum ++)
        i = intLoad (stream, &ordeptr->rangtab[cblknum]);

    /* Read direct permutation */
    for (vertnum = 0; (i == 1) && (vertnum < vertnbr); vertnum ++)
        i = intLoad (stream, &ordeptr->permtab[vertnum]);

    /* Read treetab data */
    if ( versval == 1 ) {
        for (cblknum = 0, i = 1; (i == 1) && (cblknum < cblknbr); cblknum ++)
            i = intLoad (stream, &ordeptr->treetab[cblknum]);
    }
    else {
        free( ordeptr->treetab );
        ordeptr->treetab = NULL;
    }

    if (i != 1) {
        errorPrint ("orderLoad: bad input (2)");
        orderExit  (ordeptr);
        return PASTIX_ERR_FILE;
    }

    /* Build inverse permutation */
    permtax = ordeptr->permtab - ordeptr->rangtab[0];
    peritax = ordeptr->peritab - ordeptr->rangtab[0];
    for (vertnum = ordeptr->rangtab[0], vertnnd = vertnum + vertnbr;
         vertnum < vertnnd; vertnum ++)
        peritax[permtax[vertnum]] = vertnum;

    ordeptr->baseval = ordeptr->rangtab[0];

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief Load an ordering from a file.
 *
 *******************************************************************************
 *
 * @param[inout] ordemesh
 *          The initialized ordering structure to fill in.
 *
 * @param[in] filename
 *          The filename where to read the ordering. If filename == NULL, the
 *          environment variable PASTIX_FILE_ORDER is used, and if
 *          PASTIX_FILE_ORDER is not defined, the default filename "ordername" in
 *          the current directory is used.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect,
 * @retval PASTIX_ERR_FILE if a problem occurs during the read.
 *
 *******************************************************************************/
int orderLoad( Order *ordemesh,
               char  *filename )
{
    FILE  *stream;
    int rc = PASTIX_SUCCESS;
    int env = 0;

    /* Parameter checks */
    if ( ordemesh == NULL ) {
        return PASTIX_ERR_BADPARAMETER;
    }

    /*
     * Get the environment variable as second option
     */
    if ( filename == NULL ) {
        filename = pastix_getenv( "PASTIX_FILE_ORDER" );
        env = 1;
    }

    /*
     * Get the default name as third option
     */
    if ( filename == NULL ) {
        filename = "ordername";
        env = 0;
    }

    PASTIX_FOPEN(stream, filename, "r");
    rc = ordering_load(ordemesh, stream);
    if (rc != PASTIX_SUCCESS)
    {
        errorPrint("test: cannot load order");
        EXIT(MOD_SOPALIN, PASTIX_ERR_INTERNAL);
    }
    fclose(stream);

    if (env) {
        pastix_cleanenv( filename );
    }

    return rc;
}

/**
 *******************************************************************************
 *
 * @ingroup order_dev
 *
 * @brief Save an ordering to a file.
 *
 *******************************************************************************
 *
 * @param[inout] ordeptr
 *          The ordering structure to dump to disk.
 *
 * @param[in] stream
 *          The stream where to write the ordering.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if the ordeptr structure is incorrect,
 * @retval PASTIX_ERR_FILE if a problem occurs during the write.
 *
 *******************************************************************************/
static inline int
ordering_save(const Order * const ordeptr,
              FILE * const        stream)
{
    pastix_int_t vertnbr;
    pastix_int_t vertnum;
    pastix_int_t cblknum;
    int          o;

    if (ordeptr->permtab == NULL) {
        errorPrint ("orderSave: cannot save ordering without direct permutation data");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (ordeptr->rangtab == NULL) {
        errorPrint ("orderSave: cannot save ordering without rangtab array");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (ordeptr->treetab == NULL) {
        errorPrint ("orderSave: cannot save ordering without treetab array");
        return PASTIX_ERR_BADPARAMETER;
    }

    vertnbr = ordeptr->rangtab[ordeptr->cblknbr] -  /* Get number of nodes */
        ordeptr->rangtab[0];

    assert( vertnbr == ordeptr->vertnbr );
    assert( ordeptr->rangtab[0] == ordeptr->baseval );

    if (fprintf (stream, "1\n%ld\t%ld\n",
                 (long) ordeptr->cblknbr,
                 (long) vertnbr) == EOF) {
        errorPrint ("orderSave: bad output (1)");
        return PASTIX_ERR_FILE;
    }

    /* Save column-block range array */
    for (cblknum = 0, o = 1; (o == 1) && (cblknum < ordeptr->cblknbr); cblknum ++) {
        o = intSave (stream, ordeptr->rangtab[cblknum]);
        putc (((cblknum & 7) == 7) ? '\n' : '\t', stream);
    }
    o = intSave (stream, ordeptr->rangtab[cblknum]);
    putc ('\n', stream);

    /* Save direct permutation */
    for (vertnum = 0; (o == 1) && (vertnum < (vertnbr - 1)); vertnum ++) {
        o = intSave (stream, ordeptr->permtab[vertnum]);
        putc (((vertnum & 7) == 7) ? '\n' : '\t', stream);
    }
    o = intSave (stream, ordeptr->permtab[vertnum]);
    putc ('\n', stream);

    /* Save treetab */
    for (cblknum = 0, o = 1; (o == 1) && (cblknum < ordeptr->cblknbr - 1); cblknum ++) {
        o = intSave (stream, ordeptr->treetab[cblknum]);
        putc (((cblknum & 7) == 7) ? '\n' : '\t', stream);
    }
    o = intSave (stream, ordeptr->treetab[cblknum]);
    putc ('\n', stream);

    if (o != 1) {
        errorPrint ("orderSave: bad output (2)");
        return PASTIX_ERR_FILE;
    }

    return PASTIX_SUCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief Save an ordering to a file.
 *
 *******************************************************************************
 *
 * @param[in] ordemesh
 *          The initialized ordering structure to save.
 *
 * @param[in] filename
 *          The filename where to save the ordering. If filename == NULL, the
 *          environment variable PASTIX_FILE_ORDER is used, and if
 *          PASTIX_FILE_ORDER is not defined, the default filename "ordergen" in
 *          the current directory is used.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect,
 * @retval PASTIX_ERR_FILE if a problem occurs during the write.
 *
 *******************************************************************************/
int orderSave( const Order * const ordemesh,
               char  *filename )
{
    FILE *stream;
    int rc = PASTIX_SUCCESS;
    int env = 0;

    /* Parameter checks */
    if ( ordemesh == NULL ) {
        return PASTIX_ERR_BADPARAMETER;
    }

    /*
     * Get the environment variable as second option
     */
    if ( filename == NULL ) {
        filename = pastix_getenv( "PASTIX_FILE_ORDER" );
        env = 1;
    }

    /*
     * Get the default name as third option
     */
    if ( filename == NULL ) {
        filename = "ordergen";
        env = 0;
    }

    PASTIX_FOPEN(stream, filename, "w");
    rc = ordering_save(ordemesh, stream);
    if (rc != PASTIX_SUCCESS )
    {
        errorPrint ("cannot save order");
    }
    fclose(stream);

    if (env) {
        pastix_cleanenv( filename );
    }

    return rc;
}
