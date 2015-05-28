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
#include <cblas.h>
#include "common.h"
#include "csc.h"
#include "bcsc.h"
#include "sopalin_data.h"

static pastix_complex64_t zone  =  1.0;
static pastix_complex64_t mzone = -1.0;

void
pastix_static_ztrsm( int side, int uplo, int trans, int diag,
                     sopalin_data_t *sopalin_data,
                     int nrhs, pastix_complex64_t *b, int ldb )
{
    SolverMatrix *datacode = sopalin_data->solvmtx;
    SolverCblk *cblk;
    SolverBlok *blok;
    pastix_int_t i, ii, tempm, tempn;
    pastix_int_t ttsknbr, *ttsktab;
    Task *t;

    ttsknbr = datacode->ttsknbr[0];
    ttsktab = datacode->ttsktab[0];

    /*
     *  Left / Upper / NoTrans
     */
    if (side == PastixLeft) {
        if (uplo == PastixUpper) {
            /*  We store U^t, so we swap uplo and trans */
            if (trans == PastixNoTrans) {
            }
        }
        else {
            /*
             *  Left / Lower / NoTrans
             */
            if (trans == PastixNoTrans) {
                for (ii=0; ii<ttsknbr; ii++){
                    i = ttsktab[ii];
                    t = datacode->tasktab + i;
                    cblk = datacode->cblktab + t->cblknum;

                    if ( t->taskid != COMP_1D )
                        continue;

                    tempn = cblk->lcolnum - cblk->fcolnum + 1;

                    /* Solve the diagonal block */
                    cblas_ztrsm(
                        CblasColMajor, CblasLeft, CblasLower,
                        CblasNoTrans, (enum CBLAS_DIAG)diag,
                        tempn, nrhs, CBLAS_SADDR(zone),
                        cblk->lcoeftab, cblk->stride,
                        b + cblk->lcolidx, ldb );

                    /* Apply the update */
                    for (blok = cblk->fblokptr; blok < cblk[1].fblokptr; blok++ ) {
                        SolverCblk *fcbk = datacode->cblktab + blok->cblknum;
                        tempm = blok->lrownum - blok->frownum + 1;

                        cblas_zgemm(
                            CblasColMajor, CblasNoTrans, CblasNoTrans,
                            tempm, nrhs, tempn,
                            CBLAS_SADDR(mzone), cblk->lcoeftab, cblk->stride,
                                                b + cblk->lcolidx, ldb,
                            CBLAS_SADDR(zone),  b + fcbk->lcolidx, ldb );
                    }
                }
            }
            /*
             *  Left / Lower / [Conj]Trans
             */
            else {
                for (ii=ttsknbr-1; ii>=0; ii--){
                    i = ttsktab[ii];
                    t = datacode->tasktab + i;
                    cblk = datacode->cblktab + t->cblknum;

                    if ( t->taskid != COMP_1D )
                        continue;

                    tempn = cblk->lcolnum - cblk->fcolnum + 1;

                    /* Solve the diagonal block */
                    cblas_ztrsm(
                        CblasColMajor, CblasLeft, CblasLower,
                        (enum CBLAS_TRANSPOSE)trans, (enum CBLAS_DIAG)diag,
                        tempn, nrhs, CBLAS_SADDR(zone),
                        cblk->lcoeftab,    cblk->stride,
                        b + cblk->lcolidx, ldb );

                    /* Apply the update */
                    /* for (blok = cblk->fblokptr; blok < cblk[1].fcblkptr; blok++ ) { */
                    /*     SolverCblk *fcbk = datacode->cblktab + blok->cblknum; */
                    /*     tempm = blok->lrownum - fcblk->stride - tempn; */

                    /*     cblas_zgemm( */
                    /*         CblasColMajor, CblasNoTrans, CblasNoTrans, */
                    /*         tempm, nrhs, tempn, */
                    /*         CBLAS_SADDR(mzone), cblk->lcoeftab, cblk->stride, */
                    /*                             b + cblk->lcolidx, ldb, */
                    /*         CBLAS_SADDR(zone),  b + fcbk->lcolidx, ldb ); */
                    /* } */
                }
            }
        }
    }
}
