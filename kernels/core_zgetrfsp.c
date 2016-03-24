/**
 *
 * @file core_zgetrfsp.c
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
#include "pastix_zcores.h"
#include <cblas.h>
#include "../blend/solver.h"
#include <lapacke.h>

static pastix_complex64_t zone  =  1.;
static pastix_complex64_t mzone = -1.;
static pastix_complex64_t zzero =  0.;

#if defined(PASTIX_WITH_HODLR)
#include "cHODLR_Matrix.h"
#endif /* defined(PASTIX_WITH_HODLR) */

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zgetf2sp - Computes the sequential static pivoting LU factorization of
 * the matrix m-by-n A = L * U.
 *
 *******************************************************************************
 *
 * @param[in] m
 *          The number of rows and columns of the matrix A.
 *
 * @param[in] n
 *          The number of rows and columns of the matrix A.
 *
 * @param[in,out] A
 *          The matrix A to factorize with LU factorization. The matrix
 *          is of size lda -by- n.
 *
 * @param[in] lda
 *          The leading dimension of the matrix A.
 *
 * @param[in,out] nbpivot
 *          Pointer to the number of piovting operations made during
 *          factorization. It is updated during this call
 *
 * @param[in] criteria
 *          Threshold use for static pivoting. If diagonal value is under this
 *          threshold, its value is replaced by the threshold and the number of
 *          pivots is incremented.
 *
 *******************************************************************************
 *
 * @return
 *          This routine will fail if it discovers a 0. on the diagonal during
 *          factorization.
 *
 *******************************************************************************/
static void core_zgetf2sp(pastix_int_t        m,
                          pastix_int_t        n,
                          pastix_complex64_t *A,
                          pastix_int_t        lda,
                          pastix_int_t       *nbpivot,
                          double              criteria )
{
    pastix_int_t k, minMN;
    pastix_complex64_t *Akk, *Aik, alpha;

    minMN = pastix_imin( m, n );

    Akk = A;
    for (k=0; k<minMN; k++) {
        Aik = Akk + 1;

        if ( cabs(*Akk) < criteria ) {
            (*Akk) = (pastix_complex64_t)criteria;
            (*nbpivot)++;
        }

        /* A_ik = A_ik / A_kk, i = k+1 .. n */
        alpha = 1. / (*Akk);
        cblas_zscal(m-k-1, CBLAS_SADDR( alpha ), Aik, 1 );

        if ( k+1 < minMN ) {

            /* A_ij = A_ij - A_ik * A_kj, i,j = k+1..n */
            cblas_zgeru(CblasColMajor, m-k-1, n-k-1,
                        CBLAS_SADDR(mzone),
                        Aik,        1,
                        Akk+lda, lda,
                        Aik+lda, lda);
        }

        Akk += lda+1;
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zgetrfsp - Computes the block static pivoting LU factorization of
 * the matrix m-by-n A = L * U.
 *
 *******************************************************************************
 *
 * @param[in] m
 *          The number of rows and columns of the matrix A.
 *
 * @param[in] n
 *          The number of rows and columns of the matrix A.
 *
 * @param[in,out] A
 *          The matrix A to factorize with LU factorization. The matrix
 *          is of size lda -by- n.
 *
 * @param[in] lda
 *          The leading dimension of the matrix A.
 *
 * @param[in,out] nbpivot
 *          Pointer to the number of piovting operations made during
 *          factorization. It is updated during this call
 *
 * @param[in] criteria
 *          Threshold use for static pivoting. If diagonal value is under this
 *          threshold, its value is replaced by the threshold and the nu,ber of
 *          pivots is incremented.
 *
 *******************************************************************************
 *
 * @return
 *          This routine will fail if it discovers a 0. on the diagonal during
 *          factorization.
 *
 *******************************************************************************/
#define MAXSIZEOFBLOCKS 64

static void core_zgetrfsp(pastix_int_t        n,
                          pastix_complex64_t *A,
                          pastix_int_t        lda,
                          pastix_int_t       *nbpivot,
                          double              criteria)
{
    pastix_int_t k, blocknbr, blocksize, matrixsize, tempm;
    pastix_complex64_t *Akk, *Lik, *Ukj, *Aij;

    blocknbr = (pastix_int_t) ceil( (double)n/(double)MAXSIZEOFBLOCKS );

    Akk = A; /* Lk,k     */

    for (k=0; k<blocknbr; k++) {

        tempm = n - k * MAXSIZEOFBLOCKS;
        blocksize = pastix_imin(MAXSIZEOFBLOCKS, tempm);
        Lik = Akk + blocksize;
        Ukj = Akk + blocksize*lda;
        Aij = Ukj + blocksize;

        /* Factorize the diagonal block Akk*/
        core_zgetf2sp( tempm, blocksize, Akk, lda, nbpivot, criteria );

        matrixsize = tempm - blocksize;
        if ( matrixsize > 0 ) {

            /* Compute the column Ukk+1 */
            cblas_ztrsm(CblasColMajor,
                        CblasLeft, CblasLower,
                        CblasNoTrans, CblasUnit,
                        blocksize, matrixsize,
                        CBLAS_SADDR(zone), Akk, lda,
                        Ukj, lda);

            /* Update Ak+1,k+1 = Ak+1,k+1 - Lk+1,k*Uk,k+1 */
            cblas_zgemm(CblasColMajor,
                        CblasNoTrans, CblasNoTrans,
                        matrixsize, matrixsize, blocksize,
                        CBLAS_SADDR(mzone), Lik, lda,
                        Ukj, lda,
                        CBLAS_SADDR(zone),  Aij, lda);
        }

        Akk += blocksize * (lda+1);
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zgetrfsp1d_getrf - Computes the LU factorization of one panel.
 *
 *******************************************************************************
 *
 * @param[in] cblk
 *          Pointer to the structure representing the panel to factorize in the
 *          cblktab array.  Next column blok must be accessible through cblk[1].
 *
 * @param[in,out] L
 *          The pointer to the lower matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in,out] U
 *          The pointer to the upper matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in] criteria
 *          Threshold use for static pivoting. If diagonal value is under this
 *          threshold, its value is replaced by the threshold and the nu,ber of
 *          pivots is incremented.
 *
 *******************************************************************************
 *
 * @return
 *          This routine will fail if it discovers a 0. on the diagonal during
 *          factorization.
 *
 *******************************************************************************/
static pastix_int_t compute_cblklevel( pastix_int_t cblknum )
{
    pastix_int_t father = treetab[cblknum];
    if ( father == -1 ) {
        return 1;
    }
    else {
        return compute_cblklevel( father ) + 1;
    }
}

int core_zgetrfsp1d_getrf( SolverCblk         *cblk,
                           pastix_complex64_t *D,
                           pastix_complex64_t *nul,
                           double              criteria)
{
    (void) nul;
    pastix_int_t ncols, n;
    pastix_int_t nbpivot = 0;

    ncols = cblk->lcolnum - cblk->fcolnum + 1;

    /* n corresponds to the size of diagonal block (n-by-n supported) */
    n = cblk->lcolnum - cblk->fcolnum + 1;

    /* Level in the elimination tree of the supernode beeing factorized */
    pastix_int_t level = compute_cblklevel( current_cblk );
    /* printf("Current %ld size %ld\n", current_cblk, cblk->lcolnum - cblk->fcolnum); */
    (void) level;

    /* check if diagonal column block */
    assert( cblk->fcolnum == cblk->fblokptr->frownum );
    assert( cblk->lcolnum == cblk->fblokptr->lrownum );

    pastix_complex64_t tolerance;
    pastix_int_t       levels, split, threshold, hodlrtree;

    char *lvls = getenv("LEVELS");
    levels     = atoi(lvls);
    char *tol  = getenv("TOLERANCE");
    tolerance  = atof(tol);
    char *splt = getenv("SPLITSIZE");
    split      = atoi(splt);
    char *thr  = getenv("THRESHOLD");
    threshold  = atoi(thr);
    char *tree = getenv("HODLRTREE");
    hodlrtree  = atoi(tree);

    (void) levels;
    (void) tolerance;
    (void) split;
    (void) threshold;
    (void) hodlrtree;

    current_cblk++;

#if defined(PASTIX_WITH_HODLR)
    if (current_cblk == 0){
        printf(" Size | Level | Sparse Storage | HODLR_Storage | Dense_Storage |   Gain (Mo)   | HODLR compression+factorization | Dense factorization | TRSM on U | TRSM on L | Bloks size\n");
    }

    if (n > split && level <= levels){
        Clock timer;
        clockStart(timer);

        /* Print size, level */
        printf("%5ld  %5ld ", n, level);

        /* Get sparse matrix size */
        pastix_int_t sparse_storage = 0;
        pastix_int_t i, j;
        for (i=0; i<n; i++){
            for (j=0; j<n; j++){
                if (cabs(D[n*i+j]) != 0.0)
                    sparse_storage++;
            }
        }
        printf("%9ld   ", sparse_storage);


        /* Add sparse matrix and low-rank updates */
        /* pastix_complex64_t *D_HODLR = cblk->dcoeftab_HODLR; */
        /* core_zgeadd( CblasNoTrans, n, n, 1.0, */
        /*              D_HODLR, n, */
        /*              D, n ); */


        cHODLR schur_H;

        /* Number of partitions returned by SplitSymbol() */
        pastix_int_t split_size = cblk->split_size;

        pastix_int_t nb_parts;
        pastix_int_t *parts;

        if (split_size == 0){
            schur_H = newHODLR(n, n, D, threshold, "SVD", &nb_parts, &parts);
        }
        else{
            if (hodlrtree == 0){
                schur_H = newHODLR(n, n, D, threshold, "SVD", &nb_parts, &parts);
            }
            else{
                schur_H = newHODLR_Tree(n, n, D, threshold, "SVD",
                                        split_size, cblk->split,
                                        &nb_parts, &parts);
            }
        }

        cblk->nb_parts = nb_parts;
        cblk->parts    = parts;

        gain_D += cStoreLR(schur_H, tolerance, 1);

        cblk->is_HODLR = 1;
        cblk->cMatrix  = schur_H;

        clockStop(timer);

        /* Print compression+factorization time */
        printf("%32.3g ", (double)clockVal(timer));
    }
#endif /* defined(PASTIX_WITH_HODLR) */
    /* else */
    {
        Clock timer;
        clockStart(timer);

        core_zgetrfsp(ncols, D, n, &nbpivot, criteria);

        clockStop(timer);

#if defined(PASTIX_WITH_HODLR)
        if (cblk->is_HODLR == 1){
            /* Print dense factorisation time */
            printf("%20.3g ", (double)clockVal(timer));
        }
#endif /* defined(PASTIX_WITH_HODLR) */
    }

    /* Symmetric case */
    double local_memory = n*n*8/1000000.;
    local_memory += (cblk->stride - n)*n*8/1000000.;
    total_memory += local_memory;

    /* Unsymmetric case */
    local_memory   = n*n*8/1000000.;
    local_memory  += 2*(cblk->stride - n)*n*8/1000000.;
    total_memory2 += local_memory;
    return nbpivot;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zgetrfsp1d_trsm - Apply all the trsm updates on one panel.
 *
 *******************************************************************************
 *
 * @param[in] cblk
 *          Pointer to the structure representing the panel to factorize in the
 *          cblktab array.  Next column blok must be accessible through cblk[1].
 *
 * @param[in,out] L
 *          The pointer to the lower matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in,out] U
 *          The pointer to the upper matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width
 *
 *******************************************************************************
 *
 * @return
 *         \retval PASTIX_SUCCESS on successful exit.
 *
 *******************************************************************************/
int core_zgetrfsp1d_trsm( SolverCblk         *cblk,
                          pastix_complex64_t *L,
                          pastix_complex64_t *U)
{
    (void)L;
    SolverBlok *fblok, *lblok;
    pastix_complex64_t *fU;
    pastix_int_t dima, dimb, stride;

    dima    = cblk->lcolnum - cblk->fcolnum + 1;
    stride  = cblk->stride;
    fblok = cblk->fblokptr;   /* this diagonal block */
    lblok = cblk[1].fblokptr; /* the next diagonal block */

    /* vertical dimension */
    dimb = stride - dima;

    pastix_complex64_t *D = cblk->dcoeftab;
    pastix_int_t stride_D = cblk->lcolnum - cblk->fcolnum + 1;

    /* To compress by facing cblk */
#if defined(PASTIX_WITH_HODLR)
    char *face          = getenv("FACING");
    pastix_int_t facing = atoi(face);

    if ( fblok+1 < lblok && facing == 1)
    {
        SymbolBlok *blok;
        pastix_int_t totalsize = 0;
        pastix_int_t nb_bloks  = 0;

        pastix_int_t facing    = 0;
        pastix_int_t localsize = 0;
        pastix_int_t current   = 0;

        if ( fblok+1 < lblok ){
            while (totalsize != dimb){

                blok = &fblok[nb_bloks+1];
                pastix_int_t bloksize  = blok->lrownum - blok->frownum + 1;

                if (blok->fcblknm != facing){
                    if (facing != 0){
                        fU = U + fblok[1].coefind + current;
                        if (cblk->is_HODLR == 1){
                            int rank;
                            gain_U += compress_SVD(fU, dima, localsize, &rank);
                        }
                    }

                    current  += localsize;
                    facing    = blok->fcblknm;
                    localsize = bloksize;
                }
                else{
                    localsize += bloksize;
                }

                totalsize += bloksize;
                nb_bloks++;
            }

            if (facing != 0){
                fU = U + fblok[1].coefind + current;

                if (cblk->is_HODLR == 1){
                    int rank;
                    gain_U += compress_SVD(fU, dima, localsize, &rank);
                }
            }
        }
    }
#endif /* defined(PASTIX_WITH_HODLR) */

    /* if there is an extra-diagonal bloc in column block */
    if ( fblok+1 < lblok )
    {
        /* first extra-diagonal bloc in column block address */
        fU = U + fblok[1].coefind;

#if defined(PASTIX_WITH_HODLR)
        if (cblk->is_HODLR == 1 && facing == 0){
            int rank;
            gain_U += compress_SVD(fU, dima, dimb, &rank);
        }

        /* Solve on U */
        if (cblk->is_HODLR == 0)
#endif /* defined(PASTIX_WITH_HODLR) */
        {

            pastix_int_t total = dimb;
            SolverBlok *blok   = fblok;
            while (total != 0){
                blok++;
                pastix_int_t dimb = blok->lrownum - blok->frownum + 1;
                total -= dimb;

                /* Update the data blok */
                fU = U + blok->coefind;

                pastix_int_t rank = blok->rankU;
                if (rank != -1){
                    /* printf("TRSM LR version\n"); */
                    /*        blok->rankU, current_cblk); */

                    pastix_complex64_t *u = blok->coefU_u_LR;
                    pastix_complex64_t *v = blok->coefU_v_LR;

                    cblas_ztrsm(CblasColMajor,
                                CblasRight, CblasLower,
                                CblasTrans, CblasUnit,
                                rank, dima,
                                CBLAS_SADDR(zone), D, stride_D,
                                v, dima);
                    cblas_ztrsm(CblasColMajor,
                                CblasRight, CblasUpper,
                                CblasTrans, CblasNonUnit,
                                rank, dima,
                                CBLAS_SADDR(zone), D, stride_D,
                                v, dima);
                }
                else{
                    cblas_ztrsm(CblasColMajor,
                                CblasRight, CblasLower,
                                CblasTrans, CblasUnit,
                                dimb, dima,
                                CBLAS_SADDR(zone), D, stride_D,
                                fU, stride);
                    cblas_ztrsm(CblasColMajor,
                                CblasRight, CblasUpper,
                                CblasTrans, CblasNonUnit,
                                dimb, dima,
                                CBLAS_SADDR(zone), D, stride_D,
                                fU, stride);
                }
            }
        }
#if defined(PASTIX_WITH_HODLR)
        else{
            cLU_TRSM(cblk->cMatrix, fU, dimb, dima, stride);
        }
#endif /* defined(PASTIX_WITH_HODLR) */
    }

    return PASTIX_SUCCESS;
}

int core_zgetrfsp1d_uncompress( SolverCblk         *cblk,
                                pastix_complex64_t *L,
                                pastix_complex64_t *U)
{
    (void)L;
    SolverBlok *fblok, *lblok;
    pastix_complex64_t *fU, *fL;
    pastix_int_t dima, dimb, stride;

    dima    = cblk->lcolnum - cblk->fcolnum + 1;
    stride  = cblk->stride;
    fblok = cblk->fblokptr;   /* this diagonal block */
    lblok = cblk[1].fblokptr; /* the next diagonal block */

    /* vertical dimension */
    dimb = stride - dima;

    pastix_complex64_t *D = cblk->dcoeftab;
    pastix_int_t stride_D = cblk->lcolnum - cblk->fcolnum + 1;

    /* if there is an extra-diagonal bloc in column block */
    if ( fblok+1 < lblok )
    {
        pastix_int_t total = dimb;
        SolverBlok *blok   = fblok;
        while (total != 0){
            blok++;
            pastix_int_t dimb = blok->lrownum - blok->frownum + 1;
            total -= dimb;

            /* Update the data blok */
            fU = U + blok->coefind;

            double mem_dense = dima*dimb*8./1000000.;
            if (blok->rankU != -1){
                printf("Uncompress U block SIZE %ld %ld RANK %ld\n", dima, dimb, blok->rankU);
                pastix_complex64_t *u = blok->coefU_u_LR;
                pastix_complex64_t *v = blok->coefU_v_LR;

                core_z_uncompress_LR(fU, stride,
                                     dimb, dima,
                                     u, dimb,
                                     v, dima,
                                     blok->rankU);

                double mem_SVD = blok->rankU * (dima + dimb) * 8. / 1000000.;
                if (mem_SVD < mem_dense)
                    gain_U += mem_dense - mem_SVD;
            }


            fL = L + blok->coefind;

            if (blok->rankL != -1){
                printf("Uncompress L block SIZE %ld %ld RANK %ld\n", dima, dimb, blok->rankL);
                pastix_complex64_t *u = blok->coefL_u_LR;
                pastix_complex64_t *v = blok->coefL_v_LR;

                core_z_uncompress_LR(fL, stride,
                                     dimb, dima,
                                     u, dimb,
                                     v, dima,
                                     blok->rankL);

                double mem_SVD = blok->rankL * (dima + dimb) * 8. / 1000000.;
                if (mem_SVD < mem_dense)
                    gain_L += mem_dense - mem_SVD;
            }
        }
    }

    return PASTIX_SUCCESS;
}

int core_zgetrfsp1d_trsm2( SolverCblk         *cblk,
                           pastix_complex64_t *L,
                           pastix_complex64_t *U)
{
    (void)U;
    SolverBlok *fblok, *lblok;
    pastix_complex64_t *fL;
    pastix_int_t dima, dimb, stride;

    dima    = cblk->lcolnum - cblk->fcolnum + 1;
    stride  = cblk->stride;
    fblok = cblk->fblokptr;   /* this diagonal block */
    lblok = cblk[1].fblokptr; /* the next diagonal block */

    /* vertical dimension */
    dimb = stride - dima;
    (void) dimb;

    /* if there is an extra-diagonal bloc in column block */
    if ( fblok+1 < lblok )
    {
        /* first extra-diagonal bloc in column block address */
        fL = L + fblok[1].coefind;
        (void) fL;

        /* Solve on L */
#if defined(PASTIX_WITH_HODLR)
        if (cblk->is_HODLR == 1){
            int rank;
            gain_L += compress_SVD(fL, dima, dimb, &rank);
        }
#endif /* defined(PASTIX_WITH_HODLR) */

    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zgetrfsp1d - Computes the LU factorization of one panel and apply
 * all the trsm updates to this panel.
 *
 *******************************************************************************
 *
 * @param[in] cblk
 *          Pointer to the structure representing the panel to factorize in the
 *          cblktab array.  Next column blok must be accessible through cblk[1].
 *
 * @param[in,out] L
 *          The pointer to the lower matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in,out] U
 *          The pointer to the upper matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in] criteria
 *          Threshold use for static pivoting. If diagonal value is under this
 *          threshold, its value is replaced by the threshold and the nu,ber of
 *          pivots is incremented.
 *
 *******************************************************************************
 *
 * @return
 *          This routine will fail if it discovers a 0. on the diagonal during
 *          factorization.
 *
 *******************************************************************************/
int core_zgetrfsp1d_panel( SolverCblk         *cblk,
                           pastix_complex64_t *L,
                           pastix_complex64_t *U,
                           double              criteria)
{
    pastix_complex64_t *D = cblk->dcoeftab;
    pastix_int_t nbpivot  = core_zgetrfsp1d_getrf(cblk, D, NULL, criteria);

    /* Compress + TRSM U */
    core_zgetrfsp1d_trsm(cblk, L, U);
    core_zgetrfsp1d_trsm2(cblk, L, U);

    /* Compress L */
#if defined(PASTIX_WITH_HODLR)
    if (cblk->is_HODLR)
        core_zgetrfsp1d_trsm2(cblk, L, U);
#endif /* defined(PASTIX_WITH_HODLR) */

    return nbpivot;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zgetrfsp1d_gemm - Computes the Cholesky factorization of one panel and apply
 * all the trsm updates to this panel.
 *
 *******************************************************************************
 *
 * @param[in] cblk
 *          The pointer to the data structure that describes the panel from
 *          which we compute the contributions. Next column blok must be
 *          accessible through cblk[1].
 *
 * @param[in] blok
 *          The pointer to the data structure that describes the blok from which
 *          we compute the contributions.
 *
 * @param[in] fcblk
 *          The pointer to the data structure that describes the panel on
 *          which we compute the contributions. Next column blok must be
 *          accessible through fcblk[1].
 *
 * @param[in,out] L
 *          The pointer to the lower matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in,out] U
 *          The pointer to the upper matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in,out] Cl
 *          The pointer to the lower matrix storing the coefficients of the
 *          updated panel. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in,out] Cu
 *          The pointer to the upper matrix storing the coefficients of the
 *          updated panel. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in] work
 *          Temporary memory buffer.
 *
 *******************************************************************************
 *
 * @return
 *          The number of static pivoting during factorization of the diagonal block.
 *
 *******************************************************************************/
void core_zgetrfsp1d_gemm( SolverCblk         *cblk,
                           SolverBlok         *blok,
                           SolverCblk         *fcblk,
                           pastix_complex64_t *L,
                           pastix_complex64_t *U,
                           pastix_complex64_t *Cl,
                           pastix_complex64_t *Cu,
                           pastix_complex64_t *work )
{
    SolverBlok *iterblok;
    SolverBlok *fblok;
    SolverBlok *lblok;

    pastix_complex64_t *Aik, *Akj, *Aij, *C;
    pastix_complex64_t *wtmp;
    pastix_int_t stride, stridefc, indblok;
    pastix_int_t dimi, dimj, dima, dimb;

    stride  = cblk->stride;
    dima = cblk->lcolnum - cblk->fcolnum + 1;

    pastix_complex64_t *Cd = fcblk->dcoeftab;
    pastix_int_t stride_D  = fcblk->lcolnum - fcblk->fcolnum + 1;

    /* First blok */
    indblok = blok->coefind;

    dimj = blok->lrownum - blok->frownum + 1;
    dimi = stride - indblok;

    /* Matrix A = Aik */
    Aik = L + indblok;
    Akj = U + indblok;

    /*
     * Compute update on L
     */
    wtmp = work;
    cblas_zgemm( CblasColMajor, CblasNoTrans, CblasTrans,
                 dimi, dimj, dima,
                 CBLAS_SADDR(zone),  Aik,  stride,
                 Akj,  stride,
                 CBLAS_SADDR(zzero), wtmp, dimi  );

    /*
     * Add contribution to facing cblk
     */

    /* Get the first block of the distant panel */
    fblok = fcblk->fblokptr;

    /* Move the pointer to the top of the right column */
    stridefc = fcblk->stride;
    C = Cl + (blok->frownum - fcblk->fcolnum) * stridefc;

    lblok = cblk[1].fblokptr;

    /* for all following blocks in block column */
    for (iterblok=blok; iterblok<lblok; iterblok++) {

        /* Find facing blok */
        while (!is_block_inside_fblock( iterblok, fblok ))
        {
            fblok++;
            assert( fblok < fcblk[1].fblokptr );
        }

        Aij = C + fblok->coefind + iterblok->frownum - fblok->frownum;
        dimb = iterblok->lrownum - iterblok->frownum + 1;
        pastix_cblk_lock( fcblk );

        /* If the blok modifies a diagonal block */
        if (fblok->coefind + iterblok->frownum - fblok->frownum < stride_D){

            char *lvls             = getenv("LEVELS");
            pastix_int_t levels    = atoi(lvls);
            char *splt             = getenv("SPLITSIZE");
            pastix_int_t split     = atoi(splt);
            char *thr              = getenv("THRESHOLD");
            pastix_int_t threshold = atoi(thr);
            char *tree             = getenv("HODLRTREE");
            pastix_int_t hodlrtree = atoi(tree);
            char *comp             = getenv("COMPRESS");
            pastix_int_t compress  = atoi(comp);

            (void) levels;
            (void) split;
            (void) threshold;
            (void) hodlrtree;
            (void) compress;

#if defined(PASTIX_WITH_HODLR)
            /* COMPRESS / UNCOMPRESS */
            /* pastix_int_t local_level = compute_cblklevel(current_cblk); */

            pastix_int_t n = fcblk->lcolnum - fcblk->fcolnum + 1;
            if (compress == 1 && n > split){
                pastix_int_t i;
                pastix_int_t  nb_parts = fcblk->nb_parts;
                pastix_int_t *parts    = fcblk->parts;

                pastix_int_t current = 0;
                while (parts[current+1] < (iterblok->frownum - fblok->frownum)){
                    current += 2;
                }

                while (parts[current] <= (iterblok->lrownum - fblok->frownum)
                       && current < 2*nb_parts){
                    pastix_int_t start, end;
                    start = pastix_imax(parts[current]   - iterblok->frownum + fblok->frownum, 0);
                    end   = pastix_imin(parts[current+1] - iterblok->frownum + fblok->frownum,
                                        iterblok->lrownum - iterblok->frownum);

                    Aij = Cd + (blok->frownum - fcblk->fcolnum) * stride_D
                        + fblok->coefind + iterblok->frownum - fblok->frownum;

                    core_zgeadd( CblasNoTrans, end-start+1, dimj, -1.0,
                                 wtmp + start, dimi,
                                 Aij  + start, stride_D );
                    current+=2;

                    fcblk->nb_contributions++;
                    fcblk->surface += (end-start+1) * dimj;
                }
            }

            /* Contribution to a non-HODLR structure */
            else
#endif /* defined(PASTIX_WITH_HODLR) */
            {
                Aij = Cd + (blok->frownum - fcblk->fcolnum) * stride_D
                    + fblok->coefind + iterblok->frownum - fblok->frownum;
                core_zgeadd( CblasNoTrans, dimb, dimj, -1.0,
                             wtmp, dimi,
                             Aij,  stride_D );
            }
        }

        /* If the block modifies an off-diagonal blok */
        else{

            Aij = C + fblok->coefind + iterblok->frownum - fblok->frownum;
            dimb = iterblok->lrownum - iterblok->frownum + 1;
            dimj = blok->lrownum - blok->frownum + 1;

            core_zgeadd( CblasNoTrans, dimb, dimj, -1.0,
                         wtmp, dimi,
                         Aij,  stridefc );

            /* No more LR, we received a dense contribution */
            fblok->rankU = -1;
        }

        wtmp += dimb;

        pastix_cblk_unlock( fcblk );
    }

    /*
     * Compute update on U
     */
    Aik = U + indblok;
    Akj = L + indblok;

    /* Restore data */
    wtmp  = work;
    fblok = fcblk->fblokptr;

    cblas_zgemm( CblasColMajor, CblasNoTrans, CblasTrans,
                 dimi, dimj, dima,
                 CBLAS_SADDR(zone),  Aik,  stride,
                 Akj,  stride,
                 CBLAS_SADDR(zzero), wtmp, dimi  );

    wtmp += blok->lrownum - blok->frownum + 1;

    /*
     * Add contribution to facing cblk
     */
    C  = Cd + (blok->frownum - fcblk->fcolnum);
    Cu = Cu + (blok->frownum - fcblk->fcolnum) * stridefc;

    /*
     * Update on L part (blocks facing diagonal block)
     */
    for (iterblok=blok+1; iterblok<lblok; iterblok++) {

        /* Find facing bloknum */
        if (!is_block_inside_fblock( iterblok, fblok ))
            break;

        Aij = C + (iterblok->frownum - fblok->frownum)*stride_D;
        dimb = iterblok->lrownum - iterblok->frownum + 1;

        pastix_cblk_lock( fcblk );
        core_zgeadd( CblasTrans, dimj, dimb, -1.0,
                     wtmp, dimi,
                     Aij,  stride_D);
        pastix_cblk_unlock( fcblk );

        /* Displacement to next block */
        wtmp += dimb;
    }

    C = Cu;

    /*
     * For all following blocks in block column,
     * Update is done on U directly
     */
    /* Keep updating on U */
    for (; iterblok<lblok; iterblok++) {

        /* Find facing bloknum */
        while (!is_block_inside_fblock( iterblok, fblok ))
        {
            fblok++;
            assert( fblok < fcblk[1].fblokptr );
        }

        dimb = iterblok->lrownum - iterblok->frownum + 1;
        Aij = C + fblok->coefind + iterblok->frownum - fblok->frownum;

        pastix_cblk_lock( fcblk );
        core_zgeadd( CblasNoTrans, dimb, dimj, -1.0,
                     wtmp, dimi,
                     Aij,  stridefc );
        pastix_cblk_unlock( fcblk );

        /* Displacement to next block */
        wtmp += dimb;
    }
}









void core_zgetrfsp1d_gemm_LR( SolverCblk         *cblk,
                              SolverBlok         *blok,
                              SolverCblk         *fcblk,
                              pastix_complex64_t *L,
                              pastix_complex64_t *U,
                              pastix_complex64_t *Cl,
                              pastix_complex64_t *Cu,
                              pastix_complex64_t *work )
{
    SolverBlok *iterblok;
    SolverBlok *fblok;
    SolverBlok *lblok;

    pastix_complex64_t *Aik, *Akj, *Aij, *C;
    pastix_int_t stride, stridefc, indblok;
    pastix_int_t dimi, dimj, dima, dimb;

    stride  = cblk->stride;
    dima = cblk->lcolnum - cblk->fcolnum + 1;

    pastix_complex64_t *Cd = fcblk->dcoeftab;
    pastix_int_t stride_D  = fcblk->lcolnum - fcblk->fcolnum + 1;

    /* First blok */
    indblok = blok->coefind;

    dimj = blok->lrownum - blok->frownum + 1;
    dimi = stride - indblok;

    /* Matrix A = Aik */
    Aik = L + indblok;
    Akj = U + indblok;

    /* Get the first block of the distant panel */
    fblok = fcblk->fblokptr;

    /* Move the pointer to the top of the right column */
    stridefc = fcblk->stride;
    C = Cl + (blok->frownum - fcblk->fcolnum) * stridefc;

    lblok = cblk[1].fblokptr;

    /* for all following blocks in block column */
    for (iterblok=blok; iterblok<lblok; iterblok++) {

        /* Find facing blok */
        while (!is_block_inside_fblock( iterblok, fblok ))
        {
            fblok++;
            assert( fblok < fcblk[1].fblokptr );
        }

        Aij = C + fblok->coefind + iterblok->frownum - fblok->frownum;
        dimb = iterblok->lrownum - iterblok->frownum + 1;
        pastix_cblk_lock( fcblk );

        Aik = L + iterblok->coefind;

        /* If the blok modifies a diagonal block */
        if (fblok->coefind + iterblok->frownum - fblok->frownum < stride_D){

            core_zproduct_lr2dense(iterblok, Aik, stride,
                                   dima, L_side,
                                   blok, Akj, stride,
                                   dima, U_side,
                                   work, dimi);

            Aij = Cd + (blok->frownum - fcblk->fcolnum) * stride_D
                + fblok->coefind + iterblok->frownum - fblok->frownum;
            core_zgeadd( CblasNoTrans, dimb, dimj, -1.0,
                         work, dimi,
                         Aij,  stride_D );
        }

        /* If the block modifies an off-diagonal blok */
        else{
            pastix_int_t rankL, rankU, rankF;
            rankL = iterblok->rankL;
            rankU = blok->rankU;
            rankF = fblok->rankL;

            /* The blok receiving a contribution is LR */
            if (rankF != -1){
                core_zproduct_lr2lr(iterblok, Aik, stride,
                                    dima, L_side,
                                    blok, Akj, stride,
                                    dima, U_side,
                                    fblok, Cl + fblok->coefind,
                                    C + fblok->coefind + iterblok->frownum - fblok->frownum,
                                    stridefc,
                                    stride_D, L_side,
                                    fcblk->fcolnum,
                                    work);
            }
            /* The blok receiving a contribution is dense */
            else{
                core_zproduct_lr2dense(iterblok, Aik, stride,
                                       dima, L_side,
                                       blok, Akj, stride,
                                       dima, U_side,
                                       work, dimi);

                Aij = C + fblok->coefind + iterblok->frownum - fblok->frownum;

                core_zgeadd( CblasNoTrans, dimb, dimj, -1.0,
                             work, dimi,
                             Aij,  stridefc );
            }
        }

        pastix_cblk_unlock( fcblk );
    }

    /*
     * Compute update on U
     */
    Aik = U + indblok;
    Akj = L + indblok;

    /* Restore data */
    fblok = fcblk->fblokptr;

    /*
     * Add contribution to facing cblk
     */
    C  = Cd + (blok->frownum - fcblk->fcolnum);
    /* Cu = Cu + (blok->frownum - fcblk->fcolnum) * stridefc; */

    /*
     * Update on L part (blocks facing diagonal block)
     */
    for (iterblok=blok+1; iterblok<lblok; iterblok++) {

        /* Find facing bloknum */
        if (!is_block_inside_fblock( iterblok, fblok ))
            break;

        Aij = C + (iterblok->frownum - fblok->frownum)*stride_D;
        dimb = iterblok->lrownum - iterblok->frownum + 1;
        Aik = U + iterblok->coefind;

        core_zproduct_lr2dense(iterblok, Aik, stride,
                               dima, U_side,
                               blok, Akj, stride,
                               dima, L_side,
                               work, dimi);

        pastix_cblk_lock( fcblk );
        core_zgeadd( CblasTrans, dimj, dimb, -1.0,
                     work, dimi,
                     Aij,  stride_D);
        pastix_cblk_unlock( fcblk );
    }

    /* C = Cu; */
    C = Cu + (blok->frownum - fcblk->fcolnum) * stridefc;

    /*
     * For all following blocks in block column,
     * Update is done on U directly
     */
    /* Keep updating on U */
    for (; iterblok<lblok; iterblok++) {

        /* Find facing bloknum */
        while (!is_block_inside_fblock( iterblok, fblok ))
        {
            fblok++;
            assert( fblok < fcblk[1].fblokptr );
        }

        dimb = iterblok->lrownum - iterblok->frownum + 1;
        Aij = C + fblok->coefind + iterblok->frownum - fblok->frownum;
        Aik = U + iterblok->coefind;

        pastix_int_t rankL, rankU, rankF;
        rankL = iterblok->rankL;
        rankU = blok->rankU;
        rankF = fblok->rankU;

        /* The blok receiving a contribution is LR */
        if (rankF != -1){
            core_zproduct_lr2lr(iterblok, Aik, stride,
                                dima, L_side,
                                blok, Akj, stride,
                                dima, U_side,
                                fblok, Cu + fblok->coefind,
                                C + fblok->coefind + iterblok->frownum - fblok->frownum,
                                stridefc,
                                stride_D, U_side,
                                fcblk->fcolnum,
                                work);
        }
        /* The blok receiving a contribution is dense */
        else{
            core_zproduct_lr2dense(iterblok, Aik, stride,
                                   dima, U_side,
                                   blok, Akj, stride,
                                   dima, L_side,
                                   work, dimi);

            pastix_cblk_lock( fcblk );
            core_zgeadd( CblasNoTrans, dimb, dimj, -1.0,
                         work, dimi,
                         Aij,  stridefc );
            pastix_cblk_unlock( fcblk );
        }
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_zgetrfsp1d - Computes the Cholesky factorization of one panel, apply all
 * the trsm updates to this panel, and apply all updates to the right with no
 * lock.
 *
 *******************************************************************************
 *
 * @param[in] cblk
 *          Pointer to the structure representing the panel to factorize in the
 *          cblktab array.  Next column blok must be accessible through cblk[1].
 *
 * @param[in] criteria
 *          Threshold use for static pivoting. If diagonal value is under this
 *          threshold, its value is replaced by the threshold and the nu,ber of
 *          pivots is incremented.
 *
 *******************************************************************************
 *
 * @return
 *          The number of static pivoting during factorization of the diagonal block.
 *
 *******************************************************************************/

int
core_zgetrfsp1d( SolverMatrix       *solvmtx,
                 SolverCblk         *cblk,
                 double              criteria,
                 pastix_complex64_t *work)
{
    pastix_complex64_t *L = cblk->lcoeftab;
    pastix_complex64_t *U = cblk->ucoeftab;
    SolverCblk  *fcblk;
    SolverBlok  *blok, *lblk;
    pastix_int_t nbpivot;

    /* printf("CBLK %ld STRIDE %ld\n", current_cblk, cblk->stride); */
    nbpivot = core_zgetrfsp1d_panel(cblk, L, U, criteria);

    blok = cblk->fblokptr + 1; /* this diagonal block */
    lblk = cblk[1].fblokptr;   /* the next diagonal block */

    char *splt             = getenv("SPLITSIZE");
    pastix_int_t split     = atoi(splt);
    char *thr              = getenv("THRESHOLD");
    pastix_int_t threshold = atoi(thr);
    char *tree             = getenv("HODLRTREE");
    pastix_int_t hodlrtree = atoi(tree);

    (void) split;
    (void) threshold;
    (void) hodlrtree;

    blok = cblk->fblokptr+1;
    for( ; blok < lblk; blok++ )
    {
        fcblk = (solvmtx->cblktab + blok->fcblknm);

        core_zgetrfsp1d_gemm_LR( cblk, blok, fcblk,
                                 L, U, fcblk->lcoeftab, fcblk->ucoeftab, work );

    }

    /* Uncompress for sequential_ztrsm.c still implemented in dense fashion */
    core_zgetrfsp1d_uncompress(cblk, L, U);


#if defined(PASTIX_WITH_HODLR)
    pastix_int_t n = fcblk->lcolnum - fcblk->fcolnum + 1;
    if (n > split && fcblk->is_HODLR == 0){
        printf ("\033[34;01mFirst time we touch cblk of size %ld (contrib from %ld)\033[00m\n", n, current_cblk-1);

        pastix_int_t nb_parts;
        pastix_int_t *parts;
        cHODLR schur_H;

        if (fcblk->split_size == 0){
            /* printf ("\033[34;01mNo user tree used\033[00m\n"); */
            schur_H = newHODLR(n, n, fcblk->dcoeftab, threshold,
                               "SVD", &nb_parts, &parts);
        }
        else{
            if (hodlrtree == 0){
                /* printf ("\033[34;01mNo user tree used\033[00m\n"); */
                schur_H = newHODLR(n, n, fcblk->dcoeftab, threshold,
                                   "SVD", &nb_parts, &parts);
            }
            else{
                /* printf ("\033[32;01mUser tree with %ld parts used\033[00m\n", */
                /*         fcblk->split_size); */
                schur_H = newHODLR_Tree(n, n, fcblk->dcoeftab, threshold,
                                        "SVD", fcblk->split_size, fcblk->split,
                                        &nb_parts, &parts);
            }
        }

        fcblk->is_HODLR         = 1;
        fcblk->cMatrix          = schur_H;
        fcblk->nb_parts         = nb_parts;
        fcblk->parts            = parts;
        fcblk->nb_contributions = 0;
        fcblk->surface          = 0;
        /* delHODLR(schur_H); */
    }
#endif /* defined(PASTIX_WITH_HODLR) */

#if defined(PASTIX_WITH_HODLR)
    if (cblk->is_HODLR){
        printf(" %ld contributions, surface %ld \n",
               cblk->nb_contributions, cblk->surface);
    }
#endif /* defined(PASTIX_WITH_HODLR) */

    return nbpivot;
}

