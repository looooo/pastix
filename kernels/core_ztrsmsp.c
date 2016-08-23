/**
 *
 * @file core_ztrsmsp.c
 *
 *  PaStiX kernel routines
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
#include "common.h"
#include "cblas.h"
#include "blend/solver.h"
#include "pastix_zcores.h"

static pastix_complex64_t  zone =  1.;
static pastix_complex64_t zzero =  0.;
static pastix_complex64_t mzone = -1.;

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_ztrsmsp_1d1d - Apply all the trsm updates on a panel stored in 1D
 * layout.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          If uplo == PastixLower, the contribution of:
 *          (block .. (cblk[1].fblokptr-1)) -by- block is computed and added to
 *          C, otherwise the contribution:
 *          (block+1 .. (cblk[1].fblokptr-1)) -by- block is computed and added
 *          to C.
 *          The pointer to the data structure that describes the panel from
 *          which we compute the contributions. Next column blok must be
 *          accessible through cblk[1].
 *
 * @param[in] trans
 *          Specify the transposition used for the B matrix. It has to be either
 *          PastixTrans or PastixConjTrans.
 *
 * @param[in] cblk
 *          The cblk structure to which block belongs to. The A and B pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[in] blok
 *          The block from which we compute the contributions.
 *
 * @param[in] fcblk
 *          The pointer to the data structure that describes the panel on which
 *          we compute the contributions. The C pointer must be one of the
 *          oceftab from this fcblk. Next column blok must be accessible through
 *          fcblk[1].
 *
 * @param[in] A
 *          The pointer to the coeftab of the cblk.lcoeftab matrix storing the
 *          coefficients of the panel when the Lower part is computed,
 *          cblk.ucoeftab otherwise. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in] B The pointer to the coeftab of the cblk.lcoeftab matrix storing
 *          the coefficients of the panel, if Symmetric/Hermitian cases or if
 *          upper part is computed; cblk.ucoeftab otherwise. Must be of size
 *          cblk.stride -by- cblk.width
 *
 * @param[in,out] C
 *          The pointer to the fcblk.lcoeftab if the lower part is computed,
 *          fcblk.ucoeftab otherwise.
 *
 * @param[in] work
 *          Temporary memory buffer.
 *
 *******************************************************************************
 *
 * @return
 *          The number of static pivoting during factorization of the diagonal
 *          block.
 *
 *******************************************************************************/
static inline int
core_ztrsmsp_1d( int side, int uplo, int trans, int diag,
                       SolverCblk         *cblk,
                 const pastix_complex64_t *A,
                       pastix_complex64_t *C )
{
    SolverBlok *fblok;
    pastix_int_t M, N, lda;

    N     = cblk->lcolnum - cblk->fcolnum + 1;
    lda   = cblk->stride;
    fblok = cblk->fblokptr;  /* The diagonal block */

    /* vertical dimension */
    M = lda - N;

    /* if there is an extra-diagonal bloc in column block */
    assert( fblok + 1 < cblk[1].fblokptr );
    assert( blok_rownbr( fblok) == N );
    assert(!(cblk->cblktype & CBLK_SPLIT));

    /* first extra-diagonal bloc in column block address */
    C = C + fblok[1].coefind;

    cblas_ztrsm(CblasColMajor,
                side, uplo, trans, diag,
                M, N,
                CBLAS_SADDR(zone), A, lda,
                                   C, lda);

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_ztrsmsp_2d - Computes the updates associated to one off-diagonal block
 * between two cblk stored in 2D.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          If uplo == PastixLower, the contribution of:
 *          (block .. (cblk[1].fblokptr-1)) -by- block is computed and added to
 *          C, otherwise the contribution:
 *          (block+1 .. (cblk[1].fblokptr-1)) -by- block is computed and added
 *          to C.
 *          The pointer to the data structure that describes the panel from
 *          which we compute the contributions. Next column blok must be
 *          accessible through cblk[1].
 *
 * @param[in] trans
 *          Specify the transposition used for the B matrix. It has to be either
 *          PastixTrans or PastixConjTrans.
 *
 * @param[in] cblk
 *          The cblk structure to which block belongs to. The A and B pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[in] blok
 *          The block from which we compute the contributions.
 *
 * @param[in] fcblk
 *          The pointer to the data structure that describes the panel on which
 *          we compute the contributions. The C pointer must be one of the
 *          oceftab from this fcblk. Next column blok must be accessible through
 *          fcblk[1].
 *
 * @param[in] A
 *          The pointer to the coeftab of the cblk.lcoeftab matrix storing the
 *          coefficients of the panel when the Lower part is computed,
 *          cblk.ucoeftab otherwise. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in] B The pointer to the coeftab of the cblk.lcoeftab matrix storing
 *          the coefficients of the panel, if Symmetric/Hermitian cases or if
 *          upper part is computed; cblk.ucoeftab otherwise. Must be of size
 *          cblk.stride -by- cblk.width
 *
 * @param[in,out] C
 *          The pointer to the fcblk.lcoeftab if the lower part is computed,
 *          fcblk.ucoeftab otherwise.
 *
 * @param[in] work
 *          Temporary memory buffer.
 *
 * @param[in] tol
 *          Tolerance for low-rank compression kernels
 *
 *******************************************************************************
 *
 * @return
 *          The number of static pivoting during factorization of the diagonal
 *          block.
 *
 *******************************************************************************/
static inline int
core_ztrsmsp_2d( int side, int uplo, int trans, int diag,
                       SolverCblk         *cblk,
                 const pastix_complex64_t *A,
                       pastix_complex64_t *C )
{
    const SolverBlok *fblok, *lblok, *blok;
    pastix_int_t M, N, lda, ldc;
    pastix_complex64_t *blokC;

    N     = cblk->lcolnum - cblk->fcolnum + 1;
    fblok = cblk[0].fblokptr;  /* The diagonal block */
    lblok = cblk[1].fblokptr;  /* The diagonal block of the next cblk */
    lda   = blok_rownbr( fblok );

    assert( blok_rownbr(fblok) == N );
    assert( cblk->cblktype & CBLK_SPLIT );

    for (blok=fblok+1; blok<lblok; blok++) {

        blokC = C + blok->coefind;
        M   = blok_rownbr(blok);
        ldc = M;

        cblas_ztrsm(CblasColMajor,
                    side, uplo, trans, diag,
                    M, N,
                    CBLAS_SADDR(zone), A, lda,
                                       blokC, ldc);
    }

    return PASTIX_SUCCESS;
}

int
core_ztrsmsp_2dsub( int side, int uplo, int trans, int diag,
                          SolverCblk         *cblk,
                          pastix_int_t        fcblknum,
                    const pastix_complex64_t *A,
                          pastix_complex64_t *C )
{
    const SolverBlok *fblok, *lblok, *blok;
    pastix_int_t M, N, lda, ldc, offset;
    pastix_complex64_t *Cptr;

    N     = cblk->lcolnum - cblk->fcolnum + 1;
    fblok = cblk[0].fblokptr;  /* The diagonal block */
    lblok = cblk[1].fblokptr;  /* The diagonal block of the next cblk */
    lda   = blok_rownbr( fblok );

    assert( blok_rownbr(fblok) == N );
    assert( cblk->cblktype & CBLK_SPLIT );

    blok = fblok+1;
    while( (blok->fcblknm < fcblknum) &&
           (blok < lblok) )
    {
        blok++;
    }

    offset = blok->coefind;
    for (; (blok < lblok) && (blok->fcblknm == fcblknum); blok++) {

        Cptr = C + blok->coefind - offset;
        M   = blok_rownbr(blok);
        ldc = M;

        cblas_ztrsm(CblasColMajor,
                    side, uplo, trans, diag,
                    M, N,
                    CBLAS_SADDR(zone), A, lda,
                                       Cptr, ldc);
    }

    return PASTIX_SUCCESS;
}

static inline int
core_ztrsmsp_2dlr( int coef, int side, int uplo, int trans, int diag,
                   SolverCblk *cblk )
{
    const SolverBlok *fblok, *lblok, *blok;
    pastix_int_t M, N, lda;
    pastix_lrblock_t *lrA, *lrC;
    pastix_complex64_t *A;

    N     = cblk->lcolnum - cblk->fcolnum + 1;
    fblok = cblk[0].fblokptr;  /* The diagonal block */
    lblok = cblk[1].fblokptr;  /* The diagonal block of the next cblk */

    lrA   = fblok->LRblock + coef;
    A     = lrA->u;
    lda   = lrA->rkmax;

    assert( lrA->rk == -1 );
    assert( blok_rownbr(fblok) == N );
    assert( cblk->cblktype & CBLK_SPLIT );
    assert( !(cblk->cblktype & CBLK_DENSE) );

    for (blok=fblok+1; blok<lblok; blok++) {

        lrC = blok->LRblock + coef;

        if ( lrC->rk != 0 ) {
            if ( lrC->rk != -1 ) {
                cblas_ztrsm(CblasColMajor,
                            side, uplo, trans, diag,
                            lrC->rk, N,
                            CBLAS_SADDR(zone), A, lda,
                            lrC->v, lrC->rkmax);
            }
            else {
                M = blok_rownbr(blok);
                cblas_ztrsm(CblasColMajor,
                            side, uplo, trans, diag,
                            M, N,
                            CBLAS_SADDR(zone), A, lda,
                            lrC->u, lrC->rkmax);
            }
        }
    }

    return PASTIX_SUCCESS;
}

int
core_ztrsmsp_2dlrsub( int coef, int side, int uplo, int trans, int diag,
                          SolverCblk         *cblk,
                          pastix_int_t        fcblknum )
{
    const SolverBlok *fblok, *lblok, *blok;
    pastix_int_t M, N, lda;
    pastix_complex64_t *A;
    pastix_lrblock_t *lrA, *lrC;

    N     = cblk->lcolnum - cblk->fcolnum + 1;
    fblok = cblk[0].fblokptr;  /* The diagonal block */
    lblok = cblk[1].fblokptr;  /* The diagonal block of the next cblk */

    lrA   = fblok->LRblock + coef;
    A     = lrA->u;
    lda   = lrA->rkmax;

    assert( lrA->rk == -1 );
    assert( blok_rownbr(fblok) == N );
    assert( cblk->cblktype & CBLK_SPLIT );
    assert( !(cblk->cblktype & CBLK_DENSE) );

    blok = fblok+1;
    while( (blok->fcblknm < fcblknum) &&
           (blok < lblok) )
    {
        blok++;
    }

    for (; (blok < lblok) && (blok->fcblknm == fcblknum); blok++) {

        lrC = blok->LRblock + coef;

        if ( lrC->rk != 0 ) {
            if ( lrC->rk != -1 ) {
                cblas_ztrsm(CblasColMajor,
                            side, uplo, trans, diag,
                            lrC->rk, N,
                            CBLAS_SADDR(zone), A, lda,
                            lrC->v, lrC->rkmax);
            }
            else {
                M = blok_rownbr(blok);
                cblas_ztrsm(CblasColMajor,
                            side, uplo, trans, diag,
                            M, N,
                            CBLAS_SADDR(zone), A, lda,
                            lrC->u, lrC->rkmax);
            }
        }
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_ztrsmsp - Computes the updates associated to one off-diagonal block.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          If uplo == PastixLower, the contribution of:
 *          (block .. (cblk[1].fblokptr-1)) -by- block is computed and added to
 *          C, otherwise the contribution:
 *          (block+1 .. (cblk[1].fblokptr-1)) -by- block is computed and added
 *          to C.
 *          The pointer to the data structure that describes the panel from
 *          which we compute the contributions. Next column blok must be
 *          accessible through cblk[1].
 *
 * @param[in] trans
 *          Specify the transposition used for the B matrix. It has to be either
 *          PastixTrans or PastixConjTrans.
 *
 * @param[in] cblk
 *          The cblk structure to which block belongs to. The A and B pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[in] blok
 *          The block from which we compute the contributions.
 *
 * @param[in] fcblk
 *          The pointer to the data structure that describes the panel on which
 *          we compute the contributions. The C pointer must be one of the
 *          oceftab from this fcblk. Next column blok must be accessible through
 *          fcblk[1].
 *
 * @param[in] A
 *          The pointer to the coeftab of the cblk.lcoeftab matrix storing the
 *          coefficients of the panel when the Lower part is computed,
 *          cblk.ucoeftab otherwise. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in] B The pointer to the coeftab of the cblk.lcoeftab matrix storing
 *          the coefficients of the panel, if Symmetric/Hermitian cases or if
 *          upper part is computed; cblk.ucoeftab otherwise. Must be of size
 *          cblk.stride -by- cblk.width
 *
 * @param[in,out] C
 *          The pointer to the fcblk.lcoeftab if the lower part is computed,
 *          fcblk.ucoeftab otherwise.
 *
 * @param[in] work
 *          Temporary memory buffer.
 *
 * @param[in] tol
 *          Tolerance for low-rank compression kernels
 *
 *******************************************************************************
 *
 * @return
 *          The number of static pivoting during factorization of the diagonal
 *          block.
 *
 *******************************************************************************/
void core_ztrsmsp( int coef, int side, int uplo, int trans, int diag,
                         SolverCblk         *cblk,
                   const pastix_complex64_t *A,
                         pastix_complex64_t *C )
{
    if (  cblk[0].fblokptr + 1 < cblk[1].fblokptr ) {
        if ( !(cblk->cblktype & CBLK_DENSE) ) {
            core_ztrsmsp_2dlr( coef, side, uplo, trans, diag,
                               cblk );
        }
        else {
            if ( cblk->cblktype & CBLK_SPLIT ) {
                core_ztrsmsp_2d( side, uplo, trans, diag,
                                 cblk, A, C );
            }
            else {
                core_ztrsmsp_1d( side, uplo, trans, diag,
                                 cblk, A, C );
            }
        }
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_ztrsmsp - Computes the updates associated to one off-diagonal block.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          If uplo == PastixLower, the contribution of:
 *          (block .. (cblk[1].fblokptr-1)) -by- block is computed and added to
 *          C, otherwise the contribution:
 *          (block+1 .. (cblk[1].fblokptr-1)) -by- block is computed and added
 *          to C.
 *          The pointer to the data structure that describes the panel from
 *          which we compute the contributions. Next column blok must be
 *          accessible through cblk[1].
 *
 * @param[in] trans
 *          Specify the transposition used for the B matrix. It has to be either
 *          PastixTrans or PastixConjTrans.
 *
 * @param[in] cblk
 *          The cblk structure to which block belongs to. The A and B pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[in] blok
 *          The block from which we compute the contributions.
 *
 * @param[in] fcblk
 *          The pointer to the data structure that describes the panel on which
 *          we compute the contributions. The C pointer must be one of the
 *          oceftab from this fcblk. Next column blok must be accessible through
 *          fcblk[1].
 *
 * @param[in] A
 *          The pointer to the coeftab of the cblk.lcoeftab matrix storing the
 *          coefficients of the panel when the Lower part is computed,
 *          cblk.ucoeftab otherwise. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in] B The pointer to the coeftab of the cblk.lcoeftab matrix storing
 *          the coefficients of the panel, if Symmetric/Hermitian cases or if
 *          upper part is computed; cblk.ucoeftab otherwise. Must be of size
 *          cblk.stride -by- cblk.width
 *
 * @param[in,out] C
 *          The pointer to the fcblk.lcoeftab if the lower part is computed,
 *          fcblk.ucoeftab otherwise.
 *
 * @param[in] work
 *          Temporary memory buffer.
 *
 * @param[in] tol
 *          Tolerance for low-rank compression kernels
 *
 *******************************************************************************
 *
 * @return
 *          The number of static pivoting during factorization of the diagonal
 *          block.
 *
 *******************************************************************************/
void solve_ztrsmsp( int side, int uplo, int trans, int diag,
                    SolverMatrix *datacode, SolverCblk *cblk,
                    int nrhs, pastix_complex64_t *b, int ldb )
{
    SolverCblk *fcbk;
    SolverBlok *blok;
    pastix_complex64_t *A, *tmp;
    pastix_int_t j, tempm, tempn, lda;
    pastix_lrblock_t *lrA;

    tempn = cblk->lcolnum - cblk->fcolnum + 1;
    lda = (cblk->cblktype & CBLK_SPLIT) ? blok_rownbr( cblk->fblokptr ) : cblk->stride;

    /*
     *  Left / Upper / NoTrans
     */
    if (side == PastixLeft) {
        if (uplo == PastixUpper) {
            if ( cblk->cblktype & CBLK_DENSE ) {
                A = (pastix_complex64_t*)(cblk->ucoeftab);
            } else {
                A = (pastix_complex64_t*)(cblk->fblokptr->LRblock[1].u);
                assert( cblk->fblokptr->LRblock[1].rkmax == lda );
            }

            /*  We store U^t, so we swap uplo and trans */
            if (trans == PastixNoTrans) {

                /* Solve the diagonal block */
                cblas_ztrsm(
                    CblasColMajor, CblasLeft, CblasLower,
                    CblasTrans, (enum CBLAS_DIAG)diag,
                    tempn, nrhs,
                    CBLAS_SADDR(zone), A, lda,
                                       b + cblk->lcolidx, ldb );

                /* Apply the update */
                for (j = cblk[1].brownum-1; j>=cblk[0].brownum; j-- ) {
                    blok = datacode->bloktab + datacode->browtab[j];
                    fcbk = datacode->cblktab + blok->lcblknm;
                    tempm = fcbk->lcolnum - fcbk->fcolnum + 1;
                    tempn = blok->lrownum - blok->frownum + 1;
                    A   = (pastix_complex64_t*)(fcbk->ucoeftab);
                    lda = (fcbk->cblktype & CBLK_SPLIT) ? tempn : fcbk->stride;

                    if ( ! (fcbk->cblktype & CBLK_DENSE)) {
                        lrA = blok->LRblock + 1;

                        switch (lrA->rk){
                        case 0:
                            break;
                        case -1:
                            cblas_zgemm(
                                CblasColMajor, CblasTrans, CblasNoTrans,
                                tempm, nrhs, tempn,
                                CBLAS_SADDR(mzone), lrA->u, tempn,
                                                    b + cblk->lcolidx + blok->frownum - cblk->fcolnum, ldb,
                                CBLAS_SADDR(zone),  b + fcbk->lcolidx, ldb );
                                break;
                        default:
                            MALLOC_INTERN( tmp, lrA->rk * nrhs, pastix_complex64_t);
                            cblas_zgemm(
                                CblasColMajor, CblasTrans, CblasNoTrans,
                                lrA->rk, nrhs, tempn,
                                CBLAS_SADDR(zone),  lrA->u, tempn,
                                                    b + cblk->lcolidx + blok->frownum - cblk->fcolnum, ldb,
                                CBLAS_SADDR(zzero), tmp,  lrA->rk );

                            cblas_zgemm(
                                CblasColMajor, CblasTrans, CblasNoTrans,
                                tempm, nrhs, lrA->rk,
                                CBLAS_SADDR(mzone), lrA->v, lrA->rkmax,
                                                    tmp, lrA->rk,
                                CBLAS_SADDR(zone),  b + fcbk->lcolidx, ldb );
                            memFree_null(tmp);
                            break;
                        }
                    }
                    else{
                        cblas_zgemm(
                            CblasColMajor, CblasTrans, CblasNoTrans,
                            tempm, nrhs, tempn,
                            CBLAS_SADDR(mzone), A + blok->coefind, lda,
                                                b + cblk->lcolidx + blok->frownum - cblk->fcolnum, ldb,
                            CBLAS_SADDR(zone),  b + fcbk->lcolidx, ldb );
                    }
                }
            }
            /*
             *  Left / Upper / [Conj]Trans
             */
            else {
            }
        }
        else {
            A = (pastix_complex64_t*)(cblk->lcoeftab);

            /*
             *  Left / Lower / NoTrans
             */
            if (trans == PastixNoTrans) {

                /* In sequential */
                assert( cblk->fcolnum == cblk->lcolidx );

                if ( cblk->cblktype & CBLK_DENSE ) {

                    /* Solve the diagonal block */
                    cblas_ztrsm(
                        CblasColMajor, CblasLeft, CblasLower,
                        CblasNoTrans, (enum CBLAS_DIAG)diag,
                        tempn, nrhs,
                        CBLAS_SADDR(zone), A, lda,
                                           b + cblk->lcolidx, ldb );

                    /* Apply the update */
                    for (blok = cblk[0].fblokptr+1; blok < cblk[1].fblokptr; blok++ ) {
                        fcbk  = datacode->cblktab + blok->fcblknm;
                        tempm = blok->lrownum - blok->frownum + 1;

                        assert( blok->frownum >= fcbk->fcolnum );
                        assert( tempm <= (fcbk->lcolnum - fcbk->fcolnum + 1));

                        lda = (cblk->cblktype & CBLK_SPLIT) ? tempm : cblk->stride;

                        cblas_zgemm(
                            CblasColMajor, CblasNoTrans, CblasNoTrans,
                            tempm, nrhs, tempn,
                            CBLAS_SADDR(mzone), A + blok->coefind, lda,
                                                b + cblk->lcolidx, ldb,
                            CBLAS_SADDR(zone),  b + fcbk->lcolidx + blok->frownum - fcbk->fcolnum, ldb );
                    }

                } else {

                    /* Solve the diagonal block */
                    lrA = cblk->fblokptr->LRblock;
                    cblas_ztrsm(
                        CblasColMajor, CblasLeft, CblasLower,
                        CblasNoTrans, (enum CBLAS_DIAG)diag,
                        tempn, nrhs, CBLAS_SADDR(zone),
                        lrA->u, tempn,
                        b + cblk->lcolidx, ldb );

                    /* Apply the update */
                    for (blok = cblk[0].fblokptr+1; blok < cblk[1].fblokptr; blok++ ) {
                        fcbk  = datacode->cblktab + blok->fcblknm;
                        tempm = blok->lrownum - blok->frownum + 1;
                        lrA   = blok->LRblock;

                        assert( blok->frownum >= fcbk->fcolnum );
                        assert( tempm <= (fcbk->lcolnum - fcbk->fcolnum + 1));


                        switch (lrA->rk){
                        case 0:
                            break;
                        case -1:
                            assert( lrA->rkmax == tempm );
                            cblas_zgemm(
                                CblasColMajor, CblasNoTrans, CblasNoTrans,
                                tempm, nrhs, tempn,
                                CBLAS_SADDR(mzone), lrA->u, tempm,
                                                    b + cblk->lcolidx, ldb,
                                CBLAS_SADDR(zone),  b + fcbk->lcolidx + blok->frownum - fcbk->fcolnum, ldb );
                            break;
                        default:
                            MALLOC_INTERN( tmp, lrA->rk * nrhs, pastix_complex64_t);

                            cblas_zgemm(
                                CblasColMajor, CblasNoTrans, CblasNoTrans,
                                lrA->rk, nrhs, tempn,
                                CBLAS_SADDR(zone), lrA->v, lrA->rkmax,
                                                   b + cblk->lcolidx, ldb,
                                CBLAS_SADDR(zzero), tmp, lrA->rk);

                            cblas_zgemm(
                                CblasColMajor, CblasNoTrans, CblasNoTrans,
                                tempm, nrhs, lrA->rk,
                                CBLAS_SADDR(mzone), lrA->u, tempm,
                                                    tmp, lrA->rk,
                                CBLAS_SADDR(zone),  b + fcbk->lcolidx + blok->frownum - fcbk->fcolnum, ldb );

                            memFree_null(tmp);
                        }
                    }
                }
            }
            /*
             *  Left / Lower / [Conj]Trans
             */
            else {
                assert( cblk->cblktype & CBLK_DENSE );

                /* Solve the diagonal block */
                cblas_ztrsm(
                    CblasColMajor, CblasLeft, CblasLower,
                    (enum CBLAS_TRANSPOSE)trans, (enum CBLAS_DIAG)diag,
                    tempn, nrhs,
                    CBLAS_SADDR(zone), A, lda,
                                       b + cblk->lcolidx, ldb );

                /* Apply the update */
                for (j = cblk[1].brownum-1; j>=cblk[0].brownum; j-- ) {
                    blok = datacode->bloktab + datacode->browtab[j];
                    fcbk = datacode->cblktab + blok->lcblknm;
                    tempm = fcbk->lcolnum - fcbk->fcolnum + 1;
                    tempn = blok->lrownum - blok->frownum + 1;
                    A   = (pastix_complex64_t*)(fcbk->lcoeftab);
                    lda = (fcbk->cblktype & CBLK_SPLIT) ? tempn : fcbk->stride;
                    assert( fcbk->cblktype & CBLK_DENSE );

                    cblas_zgemm(
                        CblasColMajor, (enum CBLAS_TRANSPOSE)trans, CblasNoTrans,
                        tempm, nrhs, tempn,
                        CBLAS_SADDR(mzone), A + blok->coefind, lda,
                                            b + cblk->lcolidx + blok->frownum - cblk->fcolnum, ldb,
                        CBLAS_SADDR(zone),  b + fcbk->lcolidx, ldb );
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


