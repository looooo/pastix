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
#if defined(HAVE_SCOTCH)
#include <scotch.h>
#endif
#include "csc_utils.h"

int orderLoadFiles( pastix_data_t *pastix_data, pastix_csc_t *csc )
{
    pastix_int_t *iparm    = pastix_data->iparm;
    Order        *ordemesh = pastix_data->ordemesh;
    pastix_int_t  procnum  = pastix_data->procnum;
    FILE *stream;
    int   retval = PASTIX_SUCCESS;
    int   strategy = iparm[IPARM_IO_STRATEGY];

    PASTIX_FOPEN(stream, "ordername", "r");
    if (orderLoad(ordemesh, stream) != 0)
    {
        errorPrint("test: cannot load order");
        EXIT(MOD_SOPALIN, PASTIX_ERR_INTERNAL);
    }
    fclose(stream);

    /* Load graph if required */
#if defined(HAVE_SCOTCH)
    if (PASTIX_MASK_ISTRUE(strategy, API_IO_LOAD_GRAPH))
    {
        PASTIX_FOPEN(stream, "graphname", "r");
        if (SCOTCH_graphLoad(&(ordemesh->grafmesh), stream, 0, 0) != 0) {
            errorPrint ("test: cannot load mesh");
            EXIT(MOD_SOPALIN, PASTIX_ERR_INTERNAL);
        }
        fclose (stream);
    }
#endif

    if (PASTIX_MASK_ISTRUE(strategy, API_IO_LOAD_CSC))
    {
        pastix_int_t ncol = 0;
        pastix_int_t *colptr, *rows;
        int dof = 1;

        if (csc->colptr   != NULL) memFree_null(csc->colptr);
        if (csc->rows     != NULL) memFree_null(csc->rows);
        if (csc->loc2glob != NULL) memFree_null(csc->loc2glob);

        if (procnum == 0) {
            //TODO
            PASTIX_FOPEN(stream, "cscname","r");
            retval = csc_load(&ncol, &colptr, &rows, NULL, &dof, stream);
            fclose(stream);
        }

        /* Number of columns */
        MPI_Bcast(&ncol, 1, PASTIX_MPI_INT, 0, pastix_data->pastix_comm);

        /* Colptr */
        if (procnum != 0) {
            MALLOC_INTERN(colptr, ncol+1, pastix_int_t);
        }
        MPI_Bcast(colptr, ncol+1, PASTIX_MPI_INT,
                  0, pastix_data->pastix_comm);

        /* Rows */
        if  (procnum != 0) {
            MALLOC_INTERN(rows, colptr[ncol]-colptr[0], pastix_int_t);
        }
        MPI_Bcast(rows, colptr[ncol]-colptr[0], PASTIX_MPI_INT,
                  0, pastix_data->pastix_comm);

        csc->n        = ncol;
        csc->gN       = ncol;
        csc->colptr   = colptr;
        csc->rows     = rows;
        csc->loc2glob = NULL;
    }
    return retval;
}

int orderSaveFiles( pastix_data_t *pastix_data )
{
    pastix_int_t *iparm    = pastix_data->iparm;
    Order        *ordemesh = pastix_data->ordemesh;
    pastix_int_t  procnum  = pastix_data->procnum;
    FILE *stream;
    int   retval = PASTIX_SUCCESS;
    int   strategy = iparm[IPARM_IO_STRATEGY];

    if (procnum == 0) {
        PASTIX_FOPEN(stream, "ordergen", "w");
        if (orderSave (ordemesh, stream) != 0) {
            errorPrint ("cannot save order");
            retval = PASTIX_ERR_INTERNAL;
        }
        fclose(stream);

#if defined(HAVE_SCOTCH)
        if (PASTIX_MASK_ISTRUE(strategy, API_IO_SAVE_GRAPH))
        {
            PASTIX_FOPEN(stream, "graphgen", "w");
            if (SCOTCH_graphSave (&(ordemesh->grafmesh), stream) != 0) {
                errorPrint ("cannot save graph");
                retval = PASTIX_ERR_INTERNAL;
            }
            fclose(stream);
        }
#endif
        if (PASTIX_MASK_ISTRUE(strategy, API_IO_SAVE_CSC)) {
            PASTIX_FOPEN(stream, "cscgen", "w");
            // TODO
            retval = csc_save( pastix_data->n2,
                               pastix_data->col2,
                               pastix_data->row2,
                               NULL, 1, stream);
            fclose(stream);
        }
    }

    return retval;
}
