/*
 *  File: order_compute_io.c
 *
 *  Function to get the ordering from the files, and it's opposite to save ordering.
 *
 *  Authors:
 *    Mathieu  Faverge    - faverge@labri.fr
 *    Xavier   LACOSTE    - lacoste@labri.fr
 *    Pierre   RAMET      - ramet@labri.fr
 *
 *  Dates:
 *    Version 0.0 - from 08 may 1998
 *                  to   08 jan 2001
 *    Version 1.0 - from 06 jun 2002
 *                  to   06 jun 2002
 */
#include "common.h"

int orderLoadFiles( pastix_data_t *pastix_data )
{
    Order        *ordemesh = pastix_data->ordemesh;
    FILE *stream;
    int   retval = PASTIX_SUCCESS;

    PASTIX_FOPEN(stream, "ordername", "r");
    if (orderLoad(ordemesh, stream) != 0)
    {
        errorPrint("test: cannot load order");
        EXIT(MOD_SOPALIN, PASTIX_ERR_INTERNAL);
    }
    fclose(stream);

    return retval;
}

int orderSaveFiles( pastix_data_t *pastix_data )
{
    Order        *ordemesh = pastix_data->ordemesh;
    pastix_int_t  procnum  = pastix_data->procnum;
    FILE *stream;
    int   retval = PASTIX_SUCCESS;

    if (procnum == 0) {
        PASTIX_FOPEN(stream, "ordergen", "w");
        if (orderSave (ordemesh, stream) != 0) {
            errorPrint ("cannot save order");
            retval = PASTIX_ERR_INTERNAL;
        }
        fclose(stream);
    }

    return retval;
}
