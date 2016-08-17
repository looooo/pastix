/**
 *
 * @file pastix_task_sopalin.c
 *
 *  PaStiX factorization routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Pascal Henon
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 * @precisions normal z -> s d c
 *
 **/
#include "cblas.h"
#include "common.h"
#include "solver.h"
#include "bcsc.h"
#include "sopalin_data.h"
#include "pastix_zcores.h"

void
sequential_ztrsm( pastix_data_t *pastix_data, int side, int uplo, int trans, int diag,
                  sopalin_data_t *sopalin_data,
                  int nrhs, pastix_complex64_t *b, int ldb )
{
    SolverMatrix *datacode = sopalin_data->solvmtx;
    SolverCblk *cblk;
    pastix_int_t i;
    (void)pastix_data;

    /*
     *  Left / Upper / NoTrans
     */
    if (side == PastixLeft) {
        if (uplo == PastixUpper) {
            /*  We store U^t, so we swap uplo and trans */
            if (trans == PastixNoTrans) {
                cblk = datacode->cblktab + datacode->cblknbr - 1;
                for (i=0; i<datacode->cblknbr; i++, cblk--){
                    solve_ztrsmsp( side, uplo, trans, diag,
                                   datacode, cblk, nrhs, b, ldb );
                }
            }
        }
        else {
            /*
             *  Left / Lower / NoTrans
             */
            if (trans == PastixNoTrans) {
                cblk = datacode->cblktab;
                for (i=0; i<datacode->cblknbr; i++, cblk++){
                    solve_ztrsmsp( side, uplo, trans, diag,
                                   datacode, cblk, nrhs, b, ldb );
                }
            }
            /*
             *  Left / Lower / [Conj]Trans
             */
            else {
                cblk = datacode->cblktab + datacode->cblknbr - 1;
                for (i=0; i<datacode->cblknbr; i++, cblk--){
                    solve_ztrsmsp( side, uplo, trans, diag,
                                   datacode, cblk, nrhs, b, ldb );
                }
            }
        }
    }
    /**
     * Right
     */
    else {
    }
}
