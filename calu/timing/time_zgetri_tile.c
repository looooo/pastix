/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_zgetri_Tile"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_GETRF(M, N) + FMULS_GETRI( N ))
#define _FADDS (FADDS_GETRF(M, N) + FADDS_GETRI( N ))

//#define GETRI_SYNC

#include "./timing.c"

/*------------------------------------------------------------------------
 *  Check the factorization of the matrix A2
 */
#if 0
static int check_getri_factorization(PLASMA_desc *descA1, PLASMA_desc *descA2, int *IPIV)
{
    int info_factorization;
    double Rnorm, Anorm, Xnorm, Bnorm, result;
    double eps = LAPACKE_dlamch_work('e');
    PLASMA_desc        *descB, *descX;
    PLASMA_Complex64_t *b = (PLASMA_Complex64_t *)malloc((descA1->m)*sizeof(PLASMA_Complex64_t));
    PLASMA_Complex64_t *x = (PLASMA_Complex64_t *)malloc((descA1->m)*sizeof(PLASMA_Complex64_t));

    PLASMA_Desc_Create(&descB, b, PlasmaComplexDouble, descA1->mb, descA1->nb, descA1->bsiz,
                       descA1->m, 1, 0, 0, descA1->m, 1);
    PLASMA_Desc_Create(&descX, x, PlasmaComplexDouble, descA1->mb, descA1->nb, descA1->bsiz,
                       descA1->m, 1, 0, 0, descA1->m, 1);

    PLASMA_zplrnt_Tile( descX, 537 );
    PLASMA_zlacpy_Tile( PlasmaUpperLower, descX, descB);

    PLASMA_zgetrs_Tile( PlasmaNoTrans, descA2, IPIV, descX );

    Xnorm = PLASMA_zlange_Tile(PlasmaInfNorm, descX);
    Anorm = PLASMA_zlange_Tile(PlasmaInfNorm, descA1);
    Bnorm = PLASMA_zlange_Tile(PlasmaInfNorm, descB);

    PLASMA_zgemm_Tile( PlasmaNoTrans, PlasmaNoTrans,
                       (PLASMA_Complex64_t)1.,  descA1, descX,
                       (PLASMA_Complex64_t)-1., descB);

    Rnorm = PLASMA_zlange_Tile(PlasmaInfNorm, descB);

    if (getenv("PLASMA_TESTING_VERBOSE"))
      printf( "||A||_oo=%f\n||X||_oo=%f\n||B||_oo=%f\n||A X - B||_oo=%e\n", Anorm, Xnorm, Bnorm, Rnorm );

    result = Rnorm / ( (Anorm*Xnorm+Bnorm)*(descA1->m)*eps ) ;
    printf("============\n");
    printf("Checking the Residual of the solution \n");
    printf("-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps) = %e \n", result);

    if (  isnan(Xnorm) || isinf(Xnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        printf("-- The factorization is suspicious ! \n");
        info_factorization = 1;
     }
    else{
        printf("-- The factorization is CORRECT ! \n");
        info_factorization = 0;
    }
    free(x); free(b);
    PLASMA_Desc_Destroy(&descB);
    PLASMA_Desc_Destroy(&descX);

    return info_factorization;
}
#endif

/*------------------------------------------------------------------------
 *  Check the accuracy of the computed inverse
 */

static int check_getri_inverse(PLASMA_desc *descA1, PLASMA_desc *descA2, int *IPIV, double *dparam )
{
    double Rnorm, Anorm, Ainvnorm, result;
    PLASMA_Complex64_t *work = (PLASMA_Complex64_t *)malloc(descA1->n*descA1->n*sizeof(PLASMA_Complex64_t));
    double eps = LAPACKE_dlamch_work('e');
    PLASMA_desc        *descW;

    PLASMA_Desc_Create(&descW, work, PlasmaComplexDouble,  descA1->mb, descA1->nb, descA1->bsiz,
                       descA1->m, descA1->n, 0, 0, descA1->m, descA1->n);

    PLASMA_zlaset_Tile( PlasmaUpperLower, (PLASMA_Complex64_t)0., (PLASMA_Complex64_t)1., descW);
    PLASMA_zgemm_Tile( PlasmaNoTrans, PlasmaNoTrans,
                       (PLASMA_Complex64_t)-1., descA2, descA1,
                       (PLASMA_Complex64_t)1.,  descW);

    Anorm    = PLASMA_zlange_Tile(PlasmaInfNorm, descA1);
    Ainvnorm = PLASMA_zlange_Tile(PlasmaInfNorm, descA2);
    Rnorm    = PLASMA_zlange_Tile(PlasmaInfNorm, descW);

    dparam[IPARAM_ANORM] = Anorm;
    dparam[IPARAM_BNORM] = Ainvnorm;

    result = Rnorm / ( (Anorm*Ainvnorm)*descA1->m*eps ) ;
    dparam[IPARAM_RES] = Rnorm;

    if (  isnan(Ainvnorm) || isinf(Ainvnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        dparam[IPARAM_XNORM] = -1.;
    }
    else{
        dparam[IPARAM_XNORM] = 0.;
    }

    PLASMA_Desc_Destroy(&descW);
    free(work);

    return PLASMA_SUCCESS;
}

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_)
{
    PLASMA_desc descW;
    int ret = 0;
    PASTE_CODE_IPARAM_LOCALS( iparam );

    if ( M != N ) {
        fprintf(stderr, "This timing works only with M == N\n");
        return -1;
    }

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA,      1, PLASMA_Complex64_t, PlasmaComplexDouble, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA2, check, PLASMA_Complex64_t, PlasmaComplexDouble, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX( piv, 1, int, N, 1 );

    PLASMA_Alloc_Workspace_zgetri_Tile_Async(descA, &descW);
    PLASMA_zplrnt_Tile( descA, 3453 );

    if ( check ) {
        PLASMA_zlacpy_Tile( PlasmaUpperLower, descA, descA2 );
    }

    /* PLASMA ZGETRF / ZTRTRI / ZTRSMRV  */
    {
#if defined(TRACE_BY_SEQUENCE)
        PLASMA_sequence *sequence[4];
        PLASMA_request request[4] = { PLASMA_REQUEST_INITIALIZER,
                                      PLASMA_REQUEST_INITIALIZER,
                                      PLASMA_REQUEST_INITIALIZER,
                                      PLASMA_REQUEST_INITIALIZER };

        PLASMA_Sequence_Create(&sequence[0]);
        PLASMA_Sequence_Create(&sequence[1]);
        PLASMA_Sequence_Create(&sequence[2]);
        PLASMA_Sequence_Create(&sequence[3]);

        if ( ! iparam[IPARAM_ASYNC] ) {

            START_TIMING();
            PLASMA_zgetrf_Tile_Async(descA, piv, sequence[0], &request[0]);
            PLASMA_Sequence_Wait(sequence[0]);

            PLASMA_ztrtri_Tile_Async(PlasmaUpper, PlasmaNonUnit, descA, sequence[1], &request[1]);
            PLASMA_Sequence_Wait(sequence[1]);

            PLASMA_ztrsmrv_Tile_Async(PlasmaRight, PlasmaLower, PlasmaNoTrans, PlasmaUnit,
                                      (PLASMA_Complex64_t) 1.0, descA, &descW, sequence[2], &request[2]);
            PLASMA_Sequence_Wait(sequence[2]);

            PLASMA_zlaswpc_Tile_Async(descA, 1, descA->m, piv, -1, sequence[3], &request[3]);
            PLASMA_Sequence_Wait(sequence[3]);
            STOP_TIMING();

        } else {

            START_TIMING();
            PLASMA_zgetrf_Tile_Async( descA, piv, sequence[0], &request[0]);
            PLASMA_ztrtri_Tile_Async( PlasmaUpper, PlasmaNonUnit,
                                      descA, sequence[1], &request[1]);
            PLASMA_ztrsmrv_Tile_Async(PlasmaRight, PlasmaLower, PlasmaNoTrans, PlasmaUnit,
                                      (PLASMA_Complex64_t) 1.0,
                                      descA, &descW, sequence[2], &request[2]);
            PLASMA_zlaswpc_Tile_Async(descA, 1, descA->m, piv, -1,
                                      sequence[3], &request[3]);

            /* Wait for everything */
            PLASMA_Sequence_Wait(sequence[0]);
            PLASMA_Sequence_Wait(sequence[1]);
            PLASMA_Sequence_Wait(sequence[2]);
            PLASMA_Sequence_Wait(sequence[3]);
            STOP_TIMING();

        }

        PLASMA_Sequence_Destroy(sequence[0]);
        PLASMA_Sequence_Destroy(sequence[1]);
        PLASMA_Sequence_Destroy(sequence[2]);
        PLASMA_Sequence_Destroy(sequence[3]);

#else
        if ( ! iparam[IPARAM_ASYNC] ) {

            START_TIMING();
            PLASMA_zgetrf_Tile(descA, piv);
            PLASMA_ztrtri_Tile(PlasmaUpper, PlasmaNonUnit, descA);
            PLASMA_ztrsmrv_Tile(PlasmaRight, PlasmaLower, PlasmaNoTrans, PlasmaUnit,
                                (PLASMA_Complex64_t) 1.0, descA, &descW);
            PLASMA_zlaswpc_Tile(descA, 1, descA->m, piv, -1);
            STOP_TIMING();

        } else {

            PLASMA_sequence *sequence;
            PLASMA_request request[2] = { PLASMA_REQUEST_INITIALIZER,
                                          PLASMA_REQUEST_INITIALIZER };

            PLASMA_Sequence_Create(&sequence);

            START_TIMING();
            PLASMA_zgetrf_Tile_Async(descA, piv, sequence, &request[0]);
            PLASMA_zgetri_Tile_Async(descA, piv, &descW, sequence, &request[1]);
            PLASMA_Sequence_Wait(sequence);
            STOP_TIMING();

            PLASMA_Sequence_Destroy(sequence);
        }
#endif
    }

    /* Check the solution */
    if ( check )
    {
        ret = check_getri_inverse(descA2, descA, piv, dparam);

        PASTE_CODE_FREE_MATRIX( descA2 );
    }

    PASTE_CODE_FREE_MATRIX( descA );
    free(descW.mat);
    free( piv );

    return ret;
}
