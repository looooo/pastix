/**
 *
 * @file compute_z.h
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/

/***************************************************************************//**
 *  Macro for matrix conversion / Lapack interface
 **/
#define plasma_zdesc_alloc( descA, mb, nb, lm, ln, i, j, m, n, free)   \
    descA = plasma_desc_init(                                          \
        PlasmaComplexDouble, (mb), (nb), ((mb)*(nb)),                  \
        (m), (n), (i), (j), (m), (n));                                 \
    if ( plasma_desc_mat_alloc( &(descA) ) ) {                         \
        plasma_error( __func__, "plasma_shared_alloc() failed");       \
        {free;};                                                       \
        return PLASMA_ERR_OUT_OF_RESOURCES;                            \
    }

#define plasma_zooplap2tile( descA, A, mb, nb, lm, ln, i, j, m, n, seq, req, free) \
    descA = plasma_desc_init(                                           \
        PlasmaComplexDouble, (mb), (nb), ((mb)*(nb)),                   \
        (lm), (ln), (i), (j), (m), (n));                                \
    if ( plasma_desc_mat_alloc( &(descA) ) ) {                          \
        plasma_error( __func__, "plasma_shared_alloc() failed");        \
        {free;};                                                        \
        return PLASMA_ERR_OUT_OF_RESOURCES;                             \
    }                                                                   \
    plasma_parallel_call_5(                                             \
        plasma_pzlapack_to_tile,                                        \
        PLASMA_Complex64_t*, (A),                                       \
        int,                 (lm),                                      \
        PLASMA_desc,         (descA),                                   \
        PLASMA_sequence*,    (seq),                                     \
        PLASMA_request*,     (req));

#define plasma_ziplap2tile( descA, A, mb, nb, lm, ln, i, j, m, n, seq, req) \
    descA = plasma_desc_init(                                         \
        PlasmaComplexDouble, (mb), (nb), ((mb)*(nb)),                 \
        (lm), (ln), (i), (j), (m), (n));                              \
    descA.mat = A;                                                    \
    PLASMA_zgecfi_Async((lm), (ln), (A), PlasmaCM, (mb), (nb),        \
                        PlasmaCCRB, (mb), (nb), (seq), (req));


#define plasma_zooptile2lap( descA, A, mb, nb, lm, ln, seq, req)    \
    plasma_parallel_call_5(plasma_pztile_to_lapack,                 \
                           PLASMA_desc,         (descA),            \
                           PLASMA_Complex64_t*, (A),                \
                           int,                 (lm),               \
                           PLASMA_sequence*,    (seq),              \
                           PLASMA_request*,     (req));

#define plasma_ziptile2lap( descA, A, mb, nb, lm, ln, seq, req)         \
    PLASMA_zgecfi_Async((lm), (ln), (A), PlasmaCCRB, (mb), (nb),        \
                        PlasmaCM, (mb), (nb), (seq), (req));

/***************************************************************************//**
 *  Declarations of parallel functions (static scheduling) - alphabetical order
 **/
void plasma_pzgeadd  (plasma_context_t *plasma);
void plasma_pzgelqf (plasma_context_t *plasma);
void plasma_pzgemm  (plasma_context_t *plasma);
void plasma_pzgeqrf (plasma_context_t *plasma);
void plasma_pzgerbb (plasma_context_t *plasma);
void plasma_pzgetmi2(plasma_context_t *plasma);
void plasma_pzgetrf_incpiv(plasma_context_t *plasma);
#ifdef COMPLEX
void plasma_pzhemm  (plasma_context_t *plasma);
void plasma_pzherk  (plasma_context_t *plasma);
void plasma_pzher2k (plasma_context_t *plasma);
#endif
void plasma_pzlacpy (plasma_context_t *plasma);
void plasma_pzlag2c (plasma_context_t *plasma);
void plasma_pzlange (plasma_context_t *plasma);
#ifdef COMPLEX
void plasma_pzlanhe (plasma_context_t *plasma);
#endif
void plasma_pzlansy (plasma_context_t *plasma);
void plasma_pzpack  (plasma_context_t *plasma);
void plasma_pzplghe (plasma_context_t *plasma);
void plasma_pzplgsy (plasma_context_t *plasma);
void plasma_pzplrnt (plasma_context_t *plasma);
void plasma_pzpotrf (plasma_context_t *plasma);
void plasma_pzshift (plasma_context_t *plasma);
void plasma_pzsymm  (plasma_context_t *plasma);
void plasma_pzsyrk  (plasma_context_t *plasma);
void plasma_pzsyr2k (plasma_context_t *plasma);
void plasma_pztrmm  (plasma_context_t *plasma);
void plasma_pztrsm  (plasma_context_t *plasma);
void plasma_pztrsmpl(plasma_context_t *plasma);
void plasma_pztrsmrv(plasma_context_t *plasma);
void plasma_pzunglq (plasma_context_t *plasma);
void plasma_pzungqr (plasma_context_t *plasma);
void plasma_pzungqrrh(plasma_context_t *plasma);
void plasma_pzunmlq (plasma_context_t *plasma);
void plasma_pzunmqr (plasma_context_t *plasma);
void plasma_pzunpack(plasma_context_t *plasma);
void plasma_pzhetrd_hb2st(plasma_context_t *plasma);
void plasma_pzunmqr_blgtrd(plasma_context_t *plasma);
void plasma_pzlarft_blgtrd(plasma_context_t *plasma);

/***************************************************************************//**
 *  Declarations of internal sequential functions 
 **/
int plasma_zshift(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A,
                  int nprob, int me, int ne, int L,
                  PLASMA_sequence *sequence, PLASMA_request *request);

/***************************************************************************//**
 *  Declarations of parallel functions (dynamic scheduling) - alphabetical order
 **/
void plasma_pzgeadd_quark(PLASMA_Complex64_t alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzbarrier_tl2pnl_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzbarrier_pnl2tl_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzbarrier_tl2row_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzbarrier_row2tl_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzgebrd_tb2bd_quark(PLASMA_enum uplo, PLASMA_desc A, double *D, double *E, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzgebrd_ge2tb_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzgelqf_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzgelqfrh_quark(PLASMA_desc A, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzgemm_quark(PLASMA_enum transA, PLASMA_enum transB, PLASMA_Complex64_t alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_Complex64_t beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);

void plasma_pzgeqp3_quark( PLASMA_desc A, int *jpvt, PLASMA_Complex64_t *tau,
                           PLASMA_Complex64_t *work, double *rwork,
                           PLASMA_sequence *sequence, PLASMA_request *request );

void plasma_pzgeqrf_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzgeqrfrh_quark(PLASMA_desc A, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzgerbh_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzgerbbrh_quark(PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzgetmi2_quark(PLASMA_enum idep, PLASMA_enum odep, PLASMA_enum storev, int m, int n, int mb, int nb, PLASMA_Complex64_t *A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzgetrf_incpiv_quark(PLASMA_desc A, PLASMA_desc L, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzgetrf_reclap_quark(PLASMA_desc A, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzgetrf_rectil_quark(PLASMA_desc A, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzhbcpy_t2bl_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_Complex64_t *AB, int LDAB, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzhegst_quark(PLASMA_enum itype, PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
#ifdef COMPLEX
void plasma_pzhemm_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_Complex64_t alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_Complex64_t beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzherk_quark(PLASMA_enum uplo, PLASMA_enum trans, double alpha, PLASMA_desc A, double beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzher2k_quark(PLASMA_enum uplo, PLASMA_enum trans, PLASMA_Complex64_t alpha, PLASMA_desc A, PLASMA_desc B, double beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
void plasma_pzhetrd_he2hb_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzhetrd_hb2st_quark(PLASMA_enum uplo, int N, int NB, int Vblksiz, PLASMA_Complex64_t *A, int LDA, PLASMA_Complex64_t *V, PLASMA_Complex64_t *TAU, double *D, double *E, int WANTZ, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzlacpy_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzlag2c_quark(PLASMA_desc A, PLASMA_desc SB, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzlange_quark(PLASMA_enum norm, PLASMA_desc A, double *work, double *result, PLASMA_sequence *sequence, PLASMA_request *request);
#ifdef COMPLEX
void plasma_pzlanhe_quark(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc A, double *work, double *result, PLASMA_sequence *sequence, PLASMA_request *request);
#endif
void plasma_pzlansy_quark(PLASMA_enum norm, PLASMA_enum uplo, PLASMA_desc A, double *work, double *result, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzlaset_quark( PLASMA_enum uplo, PLASMA_Complex64_t alpha, PLASMA_Complex64_t beta, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzlaset2_quark(PLASMA_enum uplo, PLASMA_Complex64_t alpha,                          PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzlaswp_quark(PLASMA_desc B, int *IPIV, int inc, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzlaswpc_quark(PLASMA_desc B, int *IPIV, int inc, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzlauum_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzplghe_quark(double bump, PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzplgsy_quark(PLASMA_Complex64_t bump, PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pzplrnt_quark(PLASMA_desc A, unsigned long long int seed, PLASMA_sequence *sequence, PLASMA_request *request );
void plasma_pzpotrf_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzshift_quark(int, int, int, PLASMA_Complex64_t *, int *, int, int, int, PLASMA_sequence*, PLASMA_request*);
void plasma_pzsymm_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_Complex64_t alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_Complex64_t beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzsyrk_quark(PLASMA_enum uplo, PLASMA_enum trans, PLASMA_Complex64_t alpha, PLASMA_desc A, PLASMA_Complex64_t beta,  PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzsyr2k_quark(PLASMA_enum uplo, PLASMA_enum trans, PLASMA_Complex64_t alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_Complex64_t beta, PLASMA_desc C, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pztrmm_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, PLASMA_Complex64_t alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pztrsm_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, PLASMA_Complex64_t alpha, PLASMA_desc A, PLASMA_desc B, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pztrsmpl_quark(PLASMA_desc A, PLASMA_desc B, PLASMA_desc L, int *IPIV, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pztrsmrv_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag, PLASMA_Complex64_t alpha, PLASMA_desc A, PLASMA_desc W, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pztrtri_quark(PLASMA_enum uplo, PLASMA_enum diag, PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzungbr_quark(PLASMA_enum side, PLASMA_desc A, PLASMA_desc O, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzungbrrh_quark(PLASMA_enum side, PLASMA_desc A, PLASMA_desc O, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzungqr_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzungqrrh_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzunglq_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzunglqrh_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzungtr_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzunmqr_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzunmqrrh_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzunmlq_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzunmlqrh_quark(PLASMA_enum side, PLASMA_enum trans, PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, int BS, PLASMA_sequence *sequence, PLASMA_request *request);

void plasma_pzgetrf_nopiv_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request);
void plasma_pzgetrf_tntpiv_quark(PLASMA_desc A, int *IPIV, PLASMA_desc W, int *Wi, PLASMA_sequence *sequence, PLASMA_request *request);
