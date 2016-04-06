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
#include <cblas.h>
#include "blend/solver.h"
#include "pastix_zcores.h"

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


                /* TODO: uncompress if some contributions remain */
                if (blok->Lupdates != 0 && blok->rankL != -1){
                    core_zupdate_lr(blok, L + blok->coefind,
                                    stride, dima, L_side,
                                    blok->Lxmin, blok->Lxmax, blok->Lymin, blok->Lymax);
                    blok->Lupdates = 0;
                }

                if (blok->Uupdates != 0 && blok->rankU != -1){
                    core_zupdate_lr(blok, U + blok->coefind,
                                    stride, dima, U_side,
                                    blok->Uxmin, blok->Uxmax, blok->Uymin, blok->Uymax);
                    blok->Uupdates = 0;
                }

                pastix_int_t rank = blok->rankU;
                if (rank != -1){
                    /* printf("TRSM LR version\n"); */
                    /*        blok->rankU, current_cblk); */

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
                printf("Uncompress U block SIZE %ld %ld RANK %ld (surface %.3g Mo)\n",
                       dima, dimb,
                       blok->rankU, blok->Usurface*8./1000000.);
                pastix_complex64_t *u = blok->coefU_u_LR;
                pastix_complex64_t *v = blok->coefU_v_LR;

                core_z_uncompress_LR(fU, stride,
                                     dimb, dima,
                                     u, dimb,
                                     v, dima,
                                     blok->rankU);
                blok->rankU = -1;

                double mem_SVD = blok->rankU * (dima + dimb) * 8. / 1000000.;
                if (mem_SVD < mem_dense)
                    gain_U += mem_dense - mem_SVD;

                tot_surface += blok->Usurface*8./1000000.;
            }


            fL = L + blok->coefind;

            if (blok->rankL != -1){
                printf("Uncompress L block SIZE %ld %ld RANK %ld (surface %.3g Mo)\n",
                       dima, dimb,
                       blok->rankL, blok->Lsurface*8./1000000.);
                pastix_complex64_t *u = blok->coefL_u_LR;
                pastix_complex64_t *v = blok->coefL_v_LR;

                core_z_uncompress_LR(fL, stride,
                                     dimb, dima,
                                     u, dimb,
                                     v, dima,
                                     blok->rankL);
                blok->rankL = -1;

                double mem_SVD = blok->rankL * (dima + dimb) * 8. / 1000000.;
                if (mem_SVD < mem_dense)
                    gain_L += mem_dense - mem_SVD;

                tot_surface += blok->Lsurface*8./1000000.;
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

    Clock timer;
    clockStart(timer);
    pastix_int_t nbpivot  = core_zgetrfsp1d_getrf(cblk, D, NULL, criteria);
    clockStop(timer);
    time_fact += clockVal(timer);

    /* Compress + TRSM U */
    clockStart(timer);
    core_zgetrfsp1d_trsm(cblk, L, U);
    core_zgetrfsp1d_trsm2(cblk, L, U);
    clockStop(timer);
    time_trsm += clockVal(timer);

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

    /* TODO: apply contributions by facing cblk !!! */

    /* for all following blocks in block column */
    SolverBlok *current, *iterblok2, *fblok2;
    current = NULL;
    pastix_int_t number = 0;
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

            /* The blok receiving a contribution is LR */
            if (fblok->rankL != -1){

                char *surf = getenv("SURFACE");
                pastix_int_t surface_max = atoi(surf);

                if (dimb < surface_max && dimj < surface_max){

                    pastix_int_t new_xmin, new_xmax, new_ymin, new_ymax;
                    pastix_int_t new_xsize, new_ysize;
                    new_xmin = pastix_imin(fblok->Lxmin, iterblok->frownum - fblok->frownum);
                    new_xmax = pastix_imax(fblok->Lxmax, iterblok->frownum - fblok->frownum + dimb);
                    new_ymin = pastix_imin(fblok->Lymin, blok->frownum - fcblk->fcolnum);
                    new_ymax = pastix_imax(fblok->Lymax, blok->frownum - fcblk->fcolnum + dimj);

                    new_xsize = new_xmax - new_xmin;
                    new_ysize = new_ymax - new_ymin;

                    /* Uncompress, blok would had a too large structure */
                    if ((new_xsize > surface_max || new_ysize > surface_max)
                        && fblok->Lupdates != 0){
                        /* printf("%ld updates\n", fblok->Lupdates); */
                        core_zupdate_lr(fblok, Cl + fblok->coefind,
                                        stridefc, stride_D, L_side,
                                        fblok->Lxmin, fblok->Lxmax, fblok->Lymin, fblok->Lymax);
                        fblok->Lupdates = 0;
                        fblok->Lxmin    = 1000000;
                        fblok->Lxmax    = 0;
                        fblok->Lymin    = 1000000;
                        fblok->Lymax    = 0;
                    }

                    fblok->Lupdates ++;

                    /* Columns impacted */
                    fblok->Lymin = pastix_imin(fblok->Lymin, blok->frownum - fcblk->fcolnum);
                    fblok->Lymax = pastix_imax(fblok->Lymax, blok->frownum - fcblk->fcolnum + dimj);

                    /* Rows impacted */
                    fblok->Lxmin = pastix_imin(fblok->Lxmin, iterblok->frownum - fblok->frownum);
                    fblok->Lxmax = pastix_imax(fblok->Lxmax, iterblok->frownum - fblok->frownum + dimb);

                    new_xsize = fblok->Lxmax - fblok->Lxmin;
                    new_ysize = fblok->Lymax - fblok->Lymin;
                    fblok->Lsurface = pastix_imax(fblok->Lsurface, new_xsize * new_ysize);

                    /* Add dense update */
                    core_zproduct_lr2dense(iterblok, Aik, stride,
                                           dima, L_side,
                                           blok, Akj, stride,
                                           dima, U_side,
                                           work, dimi);

                    Aij = C + fblok->coefind + iterblok->frownum - fblok->frownum;

                    /* Warning, the -1.0 will appear in adding with LR structure */
                    core_zgeadd( CblasNoTrans, dimb, dimj, 1.0,
                                 work, dimi,
                                 Aij,  stridefc );

                }

                /* LR blok receiving a contribution from a dense * dense product */
                if (dimb >= surface_max || dimj >= surface_max)
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

                /* We need to apply contributions before uncompressing */
                if (iterblok->rankL != -1 && iterblok->Lupdates != 0){
                    core_zupdate_lr(iterblok, Cl + iterblok->coefind,
                                    stridefc, stride_D, L_side,
                                    iterblok->Lxmin, iterblok->Lxmax,
                                    iterblok->Lymin, iterblok->Lymax);
                    iterblok->Lupdates = 0;
                    iterblok->Lxmin    = 1000000;
                    iterblok->Lxmax    = 0;
                    iterblok->Lymin    = 1000000;
                    iterblok->Lymax    = 0;
                }

                if (blok->rankU != -1 && blok->Uupdates != 0){
                    core_zupdate_lr(blok, Cu + blok->coefind,
                                    stridefc, stride_D, U_side,
                                    blok->Uxmin, blok->Uxmax,
                                    blok->Uymin, blok->Uymax);
                    blok->Uupdates = 0;
                    blok->Uxmin    = 1000000;
                    blok->Uxmax    = 0;
                    blok->Uymin    = 1000000;
                    blok->Uymax    = 0;
                }

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

        /* The blok receiving a contribution is LR */
        if (fblok->rankU != -1){

                char *surf = getenv("SURFACE");
                pastix_int_t surface_max = atoi(surf);

                if (dimb < surface_max && dimj < surface_max){

                    pastix_int_t new_xmin, new_xmax, new_ymin, new_ymax;
                    pastix_int_t new_xsize, new_ysize;
                    new_xmin = pastix_imin(fblok->Uxmin, iterblok->frownum - fblok->frownum);
                    new_xmax = pastix_imax(fblok->Uxmax, iterblok->frownum - fblok->frownum + dimb);
                    new_ymin = pastix_imin(fblok->Uymin, blok->frownum - fcblk->fcolnum);
                    new_ymax = pastix_imax(fblok->Uymax, blok->frownum - fcblk->fcolnum + dimj);

                    new_xsize = new_xmax - new_xmin;
                    new_ysize = new_ymax - new_ymin;

                    /* Uncompress, blok would had a too large structure */
                    if ((new_xsize > surface_max || new_ysize > surface_max)
                        && fblok->Uupdates != 0){
                        core_zupdate_lr(fblok, Cu + fblok->coefind,
                                        stridefc, stride_D, U_side,
                                        fblok->Uxmin, fblok->Uxmax, fblok->Uymin, fblok->Uymax);
                        fblok->Uupdates = 0;
                        fblok->Uxmin    = 1000000;
                        fblok->Uxmax    = 0;
                        fblok->Uymin    = 1000000;
                        fblok->Uymax    = 0;
                    }

                    fblok->Uupdates ++;

                    /* Columns impacted */
                    fblok->Uymin = pastix_imin(fblok->Uymin, blok->frownum - fcblk->fcolnum);
                    fblok->Uymax = pastix_imax(fblok->Uymax, blok->frownum - fcblk->fcolnum + dimj);

                    /* Rows impacted */
                    fblok->Uxmin = pastix_imin(fblok->Uxmin, iterblok->frownum - fblok->frownum);
                    fblok->Uxmax = pastix_imax(fblok->Uxmax, iterblok->frownum - fblok->frownum + dimb);

                    new_xsize = fblok->Uxmax - fblok->Uxmin;
                    new_ysize = fblok->Uymax - fblok->Uymin;
                    fblok->Usurface = pastix_imax(fblok->Usurface, new_xsize * new_ysize);

                    /* Add dense update */
                    core_zproduct_lr2dense(iterblok, Aik, stride,
                                           dima, U_side,
                                           blok, Akj, stride,
                                           dima, L_side,
                                           work, dimi);

                    /* Warning, the -1.0 will appear in adding with LR structure */
                    core_zgeadd( CblasNoTrans, dimb, dimj, 1.0,
                                 work, dimi,
                                 Aij,  stridefc );
                }

                /* LR blok receiving a contribution from a dense * dense product */
                if (dimb >= surface_max || dimj >= surface_max)
                    core_zproduct_lr2lr(iterblok, Aik, stride,
                                        dima, U_side,
                                        blok, Akj, stride,
                                        dima, L_side,
                                        fblok, Cu + fblok->coefind,
                                        C + fblok->coefind + iterblok->frownum - fblok->frownum,
                                        stridefc,
                                        stride_D, U_side,
                                        fcblk->fcolnum,
                                        work);
        }
        /* The blok receiving a contribution is dense */
        else{
            /* We need to apply contributions before uncompressing */
            if (blok->rankL != -1 && blok->Lupdates != 0){
                core_zupdate_lr(blok, Cl + blok->coefind,
                                stridefc, stride_D, L_side,
                                blok->Lxmin, blok->Lxmax,
                                blok->Lymin, blok->Lymax);
                blok->Lupdates = 0;
                blok->Lxmin    = 1000000;
                blok->Lxmax    = 0;
                blok->Lymin    = 1000000;
                blok->Lymax    = 0;
            }

            if (iterblok->rankU != -1 && iterblok->Uupdates != 0){
                core_zupdate_lr(iterblok, Cu + iterblok->coefind,
                                stridefc, stride_D, U_side,
                                iterblok->Uxmin, iterblok->Uxmax,
                                iterblok->Uymin, iterblok->Uymax);
                iterblok->Uupdates = 0;
                iterblok->Uxmin    = 1000000;
                iterblok->Uxmax    = 0;
                iterblok->Uymin    = 1000000;
                iterblok->Uymax    = 0;
            }

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

    nbpivot = core_zgetrfsp1d_panel(cblk, L, U, criteria);

    blok = cblk->fblokptr + 1; /* this diagonal block */
    lblk = cblk[1].fblokptr;   /* the next diagonal block */

    /* if there are off-diagonal supernodes in the column */
    for( ; blok < lblk; blok++ )
    {
        fcblk = (solvmtx->cblktab + blok->fcblknm);

        /* Update on L */
        core_zgemmsp( PastixLower, PastixTrans, cblk, blok, fcblk,
                      L, U, fcblk->lcoeftab, work );

        /* Update on U */
        if ( blok+1 < lblk ) {
            core_zgemmsp( PastixUpper, PastixTrans, cblk, blok, fcblk,
                          U, L, fcblk->ucoeftab, work );
        }
        pastix_atomic_dec_32b( &(fcblk->ctrbcnt) );
    }

    return nbpivot;
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
core_zgetrfsp1d_LR( SolverMatrix       *solvmtx,
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

    Clock timer;
    clockStart(timer);

    /* if there are off-diagonal supernodes in the column */
    for( ; blok < lblk; blok++ )
    {
        fcblk = (solvmtx->cblktab + blok->fcblknm);

        core_zgetrfsp1d_gemm_LR( cblk, blok, fcblk,
                                 L, U, fcblk->lcoeftab, fcblk->ucoeftab, work );

        pastix_atomic_dec_32b( &(fcblk->ctrbcnt) );
    }

    clockStop(timer);
    time_update += clockVal(timer);

    /* Uncompress for sequential_ztrsm.c still implemented in dense fashion */
    clockStart(timer);
    core_zgetrfsp1d_uncompress(cblk, L, U);
    clockStop(timer);
    time_uncomp += clockVal(timer);

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

