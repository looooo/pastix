/**
 *
 * @file pzhetrd_hb2st.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Azzam Haidar
 * @date 2011-05-15
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#undef REAL
#define COMPLEX

/***************************************************************************/
/**
 *  Parallel bulge chasing column-wise - static scheduling
 **/
/***************************************************************************/
void plasma_pzhetrd_hb2st(plasma_context_t *plasma)
{
    int my_core_id = PLASMA_RANK;
    int cores_num  = plasma->world_size;
    /*===========================*/
    int grsiz = 1;
    int N, NB, LDA, Vblksiz, WANTZ;
    double *D;
    double *E;
    PLASMA_enum uplo;
    PLASMA_Complex64_t *A;
    PLASMA_Complex64_t *V;
    PLASMA_Complex64_t *TAU;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    /*===========================
     *   local variables
     *===========================*/
    int nbtiles;
    PLASMA_Complex64_t *WORK;
    int coreid, sweepid, myid, shift, stt, st, ed, stind, edind;
    int blklastind, colpt, stepercol;
    int i,j,m,k;
    int thgrsiz, thgrnb, thgrid, thed;
    int colblktile,maxrequiredcores,colpercore,allcoresnb;

    /* variables used when this function is standalone */
#ifdef COMPLEX
    int standalonework = 0 ;
    static PLASMA_Complex64_t zone  = (PLASMA_Complex64_t) 1.0;
#endif

    plasma_unpack_args_13(uplo, N, NB, Vblksiz, A, LDA, V, TAU, D, E, WANTZ, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;

    if ( uplo != PlasmaLower ) {
        plasma_request_fail(sequence, request, PLASMA_ERR_NOT_SUPPORTED);
        return;
    }

    /* Quick return */
    if (N == 0) {
        return;
    }

    /*
     * The matrix is diagonal. We just copy it (abs for complex) and store it in D
     * sequential code here is better so only core 0 will work
     */
    if ( NB == 0 ) {
        if( my_core_id == 0 ){
            memset(E, 0, (N-1)*sizeof(double));
#ifdef COMPLEX
            for (i=0; i<N; i++){
                D[i] = cabs(A[LDA*i]); /* A is band storage*/
            }
#else
            for (i=0; i<N; i++){
                D[i] = A[LDA*i];
            }
#endif
        } /* end if my_core_id=0 */
        return;
    }

    /*
     * The matrix is already Bidiagonal. We do not need to apply the bulge chasing.
     * We make diagonal and superdiagonal elements real, and store them in
     * D and E.
     * If uplo == PlasmaLower, we first transform the lower bidiagonal form
     * to an upper bidiagonal form by applying plane rotations/ Householder
     * from the left, overwriting superdiagonal elements, then we make
     * elements real on the resulting upper Bidiagonal.
     * If uplo == PlasmaUpper, we make its elements real.
     * For Q, PT: ZSCAL should be done if WANTQ == 1.
     */
    /* =======================================================
     * NOTE in CASE USED STANDALONE FUNCTION
     * =======================================================*/
    /* in case where the first stage has been done we are sure
     * that all element of E are real so no need to go through
     * making the off diagonal element real just copy them to D and E.
     * However in case where this function is used for only
     * band matrices meaning only stage2 want to be called then it requires
     * to make all off-diagonal elements real and so remove this portion below
     * and uncomment the NB=1 case below.*/

    /* sequential code here so only core 0 will work */
    if ( NB == 1 ) {
        if( my_core_id == 0 ) {

            /* Case used when this function is standalone and so the off diagonal
             * element require to be real. (Only for complex) */
#ifdef COMPLEX
            if( standalonework != 0 ) {
                /*
                 * PlasmaLower
                 */
                if( uplo == PlasmaLower ) {
                    PLASMA_Complex64_t ztmp, absztmp;

                    TAU[0] = zone;
                    for (i=0; i<N-1; i++) {
                        D[i]       = creal(A[LDA*i]);
                        ztmp       = A[LDA*i+1];
                        absztmp    = cabs(ztmp);
                        A[LDA*i+1] = absztmp;
                        E[i]       = absztmp;
                        if( absztmp != 0. )
                            ztmp = (PLASMA_Complex64_t) (ztmp / (PLASMA_Complex64_t) absztmp);
                        else
                            ztmp = (PLASMA_Complex64_t)1.;
                        if( i < N-2 ) {
                            A[LDA*(i+1)+1] *= ztmp;
                        }
                        // for Q: ZSCAL should be done in case of WANTQ
                        TAU[i+1] = ztmp;
                    }
                    D[N-1] = creal(A[LDA*(N-1)]);
                }
                /*
                 * PlasmaUpper
                 */
                else {
                    plasma_error("pzhetrd_hb2st", "PlasmaUpper is not supported\n");
                    plasma_request_fail(sequence, request, PLASMA_ERR_NOT_SUPPORTED);
                    return;
                }
            }
            else
#endif
            {
                /* Case used when first stage has been applied or all real cases */
                if( uplo == PlasmaLower ){
                    for (i=0; i < N-1; i++) {
                        D[i] = creal(A[i*LDA]);
                        E[i] = creal(A[i*LDA+1]);
                    }
                    D[N-1] = creal(A[(N-1)*LDA]);
                }
                else {
                    for (i=0; i < N-1; i++) {
                        D[i] = creal(A[i*LDA+NB]);
                        E[i] = creal(A[i*LDA+NB-1]);
                    }
                    D[N-1] = creal(A[(N-1)*LDA+NB]);
                    E[N-1] = 0.;
                }
            }
        }
        return;
    }

    /*
     * General case:
     *
     * As I store V in the V vector there are overlap between
     * tasks so shift is now 4 where group need to be always
     * multiple of 2 (or shift=5 if not multiple of 2),
     * because as example if grs=1 task 2 from
     * sweep 2 can run with task 6 sweep 1., but task 2 sweep 2
     * will overwrite the V of tasks 5 sweep 1 which are used by
     * task 6, so keep in mind that group need to be multiple of 2,
     * and thus tasks 2 sweep 2 will never run with task 6 sweep 1.
     * OR if we allocate V as V(N,2) and we switch between the storing of
     * sweep's like odd in V(N,1) and even in V(N,2) then no overlap and so
     * shift is 3.
     * when storing V in matrix style, shift could be back to 3.
     * */

    /* Some tunning for the bulge chasing code
     * see technical report for details */
    nbtiles    = plasma_ceildiv(N,NB);
    allcoresnb = cores_num;
    /* grsiz   = 2;*/
    if( WANTZ == 0 ) {
        shift = 4;
    } else {
        shift = 3; /* it was 5 before see above for explanation*/
    }

    if( grsiz == 1 )
        colblktile = 1;
    else
        colblktile = grsiz/2;

    maxrequiredcores = max( nbtiles/colblktile, 1 );
    colpercore = colblktile*NB;
    allcoresnb = min( allcoresnb, maxrequiredcores );
    thgrsiz = N;

    /* allocate small workspace for the kernel to avoid an alloc/dealloc at each call */
    WORK    = (PLASMA_Complex64_t*) plasma_private_alloc(plasma, NB, PlasmaComplexDouble);
    /* Initialize static scheduler progress table */
    ss_init(2*nbtiles+shift+cores_num+10, 1, 0);


    /* main bulge chasing code */
    i = shift/grsiz;
    stepercol =  i*grsiz == shift ? i:i+1;
    i       = (N-1)/thgrsiz;
    thgrnb  = i*thgrsiz == (N-1) ? i:i+1;
    for (thgrid = 1; thgrid<=thgrnb; thgrid++){
        stt  = (thgrid-1)*thgrsiz+1;
        thed = min( (stt + thgrsiz -1), (N-1));
        for (i = stt; i <= N-1; i++){
            ed = min(i,thed);
            if(stt>ed) break;
            for (m = 1; m <=stepercol; m++){
                st=stt;
                for (sweepid = st; sweepid <=ed; sweepid++){

                    for (k = 1; k <=grsiz; k++){
                        myid = (i-sweepid)*(stepercol*grsiz) +(m-1)*grsiz + k;
                        if(myid%2 ==0){
                            colpt      = (myid/2)*NB+1+sweepid-1;
                            stind      = colpt-NB+1;
                            edind      = min(colpt,N);
                            blklastind = colpt;
                        } else {
                            colpt      = ((myid+1)/2)*NB + 1 +sweepid -1 ;
                            stind      = colpt-NB+1;
                            edind      = min(colpt,N);
                            if( (stind>=edind-1) && (edind==N) )
                                blklastind=N;
                            else
                                blklastind=0;
                        }
                        coreid = (stind/colpercore)%allcoresnb;

                        if(my_core_id==coreid) {
                            if(myid==1) {

                                ss_cond_wait(myid+shift-1, 0, sweepid-1);
                                CORE_zhbtype1cb(N, NB, A, LDA, V, TAU, stind-1, edind-1, sweepid-1, Vblksiz, WANTZ, WORK);
                                ss_cond_set(myid, 0, sweepid);

                                if(blklastind >= (N-1)) {
                                    for (j = 1; j <= shift; j++)
                                        ss_cond_set(myid+j, 0, sweepid);
                                }
                            } else {
                                ss_cond_wait(myid-1,       0, sweepid);
                                ss_cond_wait(myid+shift-1, 0, sweepid-1);
                                if(myid%2 == 0)
                                    CORE_zhbtype2cb(N, NB, A, LDA, V, TAU, stind-1, edind-1, sweepid-1, Vblksiz, WANTZ, WORK);
                                else
                                    CORE_zhbtype3cb(N, NB, A, LDA, V, TAU, stind-1, edind-1, sweepid-1, Vblksiz, WANTZ, WORK);

                                ss_cond_set(myid, 0, sweepid);
                                if(blklastind >= (N-1)) {
                                    for (j = 1; j <= shift+allcoresnb; j++)
                                        ss_cond_set(myid+j, 0, sweepid);
                                }
                            } /* END if myid==1 */
                        } /* END if my_core_id==coreid */

                        if(blklastind >= (N-1)) {
                            stt++;
                            break;
                        }
                    } /* END for k=1:grsiz */
                } /* END for sweepid=st:ed */
            } /* END for m=1:stepercol */
        } /* END for i=1:N-1 */
    } /* END for thgrid=1:thgrnb */

    /* finalize static sched */
    ss_finalize();

    /*================================================
     *  store resulting diag and lower diag D and E
     *  note that D and E are always real
     *================================================*/
    /*
     * STORE THE RESULTING diagonal/off-diagonal in D AND E
     */
    /* Make diagonal and superdiagonal elements real,
     * storing them in D and E
     */
    /* In complex case, the off diagonal element are
     * not necessary real. we have to make off-diagonal
     * elements real and copy them to E.
     * When using HouseHolder elimination,
     * the ZLARFG give us a real as output so, all the
     * diagonal/off-diagonal element except the last one are already
     * real and thus we need only to take the abs of the last
     * one.
     *  */
    /* sequential code here so only core 0 will work */
    if( my_core_id == 0 ) {
        if( uplo == PlasmaLower ) {
            for (i=0; i < N-1 ; i++) {
                D[i] = creal(A[i*LDA]);
                E[i] = creal(A[i*LDA+1]);
            }
            D[N-1] = creal(A[(N-1)*LDA]);
        } else { /* PlasmaUpper not tested yet */
            for (i=0; i<N-1; i++) {
                D[i] = creal(A[i*LDA+NB]);
                E[i] = creal(A[i*LDA+NB-1]);
            }
            D[N-1] = creal(A[(N-1)*LDA+NB]);
        } /* end PlasmaUpper */
    }

    plasma_private_free(plasma, WORK);
    return;
}
/***************************************************************************/

/***************************************************************************/
/**
 *  Parallel bulge chasing column-wise - dynamic scheduling
 **/
/***************************************************************************/
void plasma_pzhetrd_hb2st_quark(PLASMA_enum uplo, int N, int NB, int Vblksiz,
                                PLASMA_Complex64_t *A, int LDA,
                                PLASMA_Complex64_t *V, PLASMA_Complex64_t *TAU,
                                double *D, double *E, int WANTZ,
                                PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    /*===========================
     *   local variables
     *===========================*/
    int grsiz, nbtiles, blksweep, lcsweep, lcNB, blkid, lcgrsiz;
    int i;

    /* variables used when this function is standalone */
#ifdef COMPLEX
    int standalonework = 0 ;
    static PLASMA_Complex64_t zone  = (PLASMA_Complex64_t) 1.0;
#endif

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;

    if ( uplo != PlasmaLower ) {
        plasma_request_fail(sequence, request, PLASMA_ERR_NOT_SUPPORTED);
        return;
    }

    /* Quick return */
    if (N == 0) {
        return;
    }
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    /*
     * The matrix is diagonal. We just copy it (abs for complex) and store it in D
     * sequential code here is better so only core 0 will work
     */
    if ( NB == 0 ) {
        memset(E, 0, (N-1)*sizeof(double));
#ifdef COMPLEX
        for (i=0; i<N; i++){
            D[i] = cabs(A[LDA*i]); /* A is band storage*/
        }
#else
        for (i=0; i<N; i++){
            D[i] = A[LDA*i];
        }
#endif
        return;
    }

    /*
     * The matrix is already Bidiagonal. We do not need to apply the bulge chasing.
     * We make diagonal and superdiagonal elements real, and store them in
     * D and E.
     * If uplo == PlasmaLower, we first transform the lower bidiagonal form
     * to an upper bidiagonal form by applying plane rotations/ Householder
     * from the left, overwriting superdiagonal elements, then we make
     * elements real on the resulting upper Bidiagonal.
     * If uplo == PlasmaUpper, we make its elements real.
     * For Q, PT: ZSCAL should be done if WANTQ == 1.
     */
    /* =======================================================
     * NOTE in CASE USED STANDALONE FUNCTION
     * =======================================================*/
    /* in case where the first stage has been done we are sure
     * that all element of E are real so no need to go through
     * making the off diagonal element real just copy them to D and E.
     * However in case where this function is used for only
     * band matrices meaning only stage2 want to be called then it requires
     * to make all off-diagonal elements real and so remove this portion below
     * and uncomment the NB=1 case below.*/

    /* sequential code here so only core 0 will work */
    if ( NB == 1 ) {
        /* Case used when this function is standalone and so the off diagonal
         * element require to be real. (Only for complex) */
#ifdef COMPLEX
        if( standalonework != 0 ) {
            /*
             * PlasmaLower
             */
            if( uplo == PlasmaLower ) {
                PLASMA_Complex64_t ztmp, absztmp;

                TAU[0] = zone;
                for (i=0; i<N-1; i++) {
                    D[i]       = creal(A[LDA*i]);
                    ztmp       = A[LDA*i+1];
                    absztmp    = cabs(ztmp);
                    A[LDA*i+1] = absztmp;
                    E[i]       = absztmp;
                    if( absztmp != 0. )
                        ztmp = (PLASMA_Complex64_t) (ztmp / (PLASMA_Complex64_t) absztmp);
                    else
                        ztmp = (PLASMA_Complex64_t)1.;
                    if( i < N-2 ) {
                        A[LDA*(i+1)+1] *= ztmp;
                    }
                    // for Q: ZSCAL should be done in case of WANTQ
                    TAU[i+1] = ztmp;
                }
                D[N-1] = creal(A[LDA*(N-1)]);
            }
            /*
             * PlasmaUpper
             */
            else {
                plasma_error("pzhetrd_hb2st", "PlasmaUpper is not supported\n");
                plasma_request_fail(sequence, request, PLASMA_ERR_NOT_SUPPORTED);
                return;
            }
        }
        else
#endif
        {
            /* Case used when first stage has been applied or all real cases */
            if( uplo == PlasmaLower ){
                for (i=0; i < N-1; i++) {
                    D[i] = creal(A[i*LDA]);
                    E[i] = creal(A[i*LDA+1]);
                }
                D[N-1] = creal(A[(N-1)*LDA]);
            }
            else {
                for (i=0; i < N-1; i++) {
                    D[i] = creal(A[i*LDA+NB]);
                    E[i] = creal(A[i*LDA+NB-1]);
                }
                D[N-1] = creal(A[(N-1)*LDA+NB]);
                E[N-1] = 0.;
            }
        }
        return;
    }

    /*
     * General case:
     *
     * As I store V in the V vector there are overlap between
     * tasks so shift is now 4 where group need to be always
     * multiple of 2 (or shift=5 if not multiple of 2),
     * because as example if grs=1 task 2 from
     * sweep 2 can run with task 6 sweep 1., but task 2 sweep 2
     * will overwrite the V of tasks 5 sweep 1 which are used by
     * task 6, so keep in mind that group need to be multiple of 2,
     * and thus tasks 2 sweep 2 will never run with task 6 sweep 1.
     * OR if we allocate V as V(N,2) and we switch between the storing of
     * sweep's like odd in V(N,1) and even in V(N,2) then no overlap and so
     * shift is 3.
     * when storing V in matrix style, shift could be back to 3.
     * */

    /* Some tunning for the bulge chasing code
     * see technical report for details */
    nbtiles    = plasma_ceildiv(N,NB);
    /*
     * Initialisation of local parameter. those parameter should be
     * input or tuned parameter.
     */
    grsiz = 1;
    if( NB > 160 ) {
        grsiz = 1;
    }
    else if( NB > 100 ) {
        grsiz = 1;
        /*
        if( N < 5000 )
            grsiz = 1;
        else
            grsiz = 2;
        */
    } else {
        grsiz = 2;
    }
    grsiz = max(1, grsiz);

    for (blksweep = 0; blksweep<nbtiles; blksweep++){
        lcNB =  blksweep == nbtiles-1 ? N-blksweep*NB : NB;
        for (lcsweep = 0; lcsweep<lcNB; lcsweep++){
            for (blkid = blksweep; blkid<nbtiles; blkid=blkid+grsiz){
                lcgrsiz = (blkid+1) < nbtiles ? grsiz : nbtiles-blkid;
                QUARK_CORE_ztrdalg(
                    plasma->quark, &task_flags,
                    uplo, N, NB, A, LDA, V, TAU, Vblksiz, WANTZ,
                    lcgrsiz, lcsweep, blkid, blksweep);
            }
        }
    }
    /*
     * Barrier used only for now, to be sure that everything
     * is done before copying the D and E and free workspace.
     * this will be removed later when D and E are directly filled
     * during the bulge process.
     */
    QUARK_Barrier(plasma->quark);
    /*================================================
     *  store resulting diag and lower diag D and E
     *  note that D and E are always real
     *================================================*/
    /*
     * STORE THE RESULTING diagonal/off-diagonal in D AND E
     */
    /* Make diagonal and superdiagonal elements real,
     * storing them in D and E
     */
    /* In complex case, the off diagonal element are
     * not necessary real. we have to make off-diagonal
     * elements real and copy them to E.
     * When using HouseHolder elimination,
     * the ZLARFG give us a real as output so, all the
     * diagonal/off-diagonal element except the last one are already
     * real and thus we need only to take the abs of the last
     * one.
     *  */
    /* sequential code here so only core 0 will work */
    if( uplo == PlasmaLower ) {
        for (i=0; i < N-1 ; i++) {
            D[i] = creal(A[i*LDA]);
            E[i] = creal(A[i*LDA+1]);
        }
        D[N-1] = creal(A[(N-1)*LDA]);
    } else { /* PlasmaUpper not tested yet */
        for (i=0; i<N-1; i++) {
            D[i] = creal(A[i*LDA+NB]);
            E[i] = creal(A[i*LDA+NB-1]);
        }
        D[N-1] = creal(A[(N-1)*LDA+NB]);
    } /* end PlasmaUpper */

    return;
}
/***************************************************************************/
