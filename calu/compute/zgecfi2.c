/**
 *
 * @file zgecfi2.c
 *
 *  PLASMA InPlaceTransformation module
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 *  This work is the implementation of an inplace transformation 
 *  based on the GKK algorithm by Gustavson, Karlsson, Kagstrom 
 *  and its fortran implementation.
 *
 * @version 2.4.6
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 * @precisions normal z -> c d s
 *
 **/

#include <sys/types.h>
#include "common.h"
#include "zgecfi2.h"

#define PLASMA_pzgetmi2(idep, odep, storev, m, n, mb, nb, A)     \
    plasma_parallel_call_10(                                     \
        plasma_pzgetmi2,                                         \
        PLASMA_enum,  (idep),                                    \
        PLASMA_enum,  (odep),                                    \
        PLASMA_enum,  (storev),                                  \
        int,  (m),                                               \
        int,  (n),                                               \
        int,  (mb),                                              \
        int,  (nb),                                              \
        PLASMA_Complex64_t*, (A),                                \
        PLASMA_sequence*, sequence,                              \
        PLASMA_request*, request);

#define PLASMA_zshift(m, n, mb, nb, A)                          \
    plasma_zshift(plasma, (m), (n), (A),                        \
                  ( (n) / (nb) ), ( (m) / (mb) ), (nb), (mb),   \
                  sequence, request);

#define PLASMA_zshiftr(m, n, mb, nb, A)                         \
    plasma_zshift(plasma, (m), (n), (A),                        \
                  ( (n) / (nb) ), (nb), ( (m) / (mb) ), (mb),   \
                  sequence, request);

/** ****************************************************************************
 *
 * @ingroup InPlaceTransformation
 *
 *  ipt_zcm2ccrb converts a matrix from CM format to CCRB format
 *
 *******************************************************************************
 *
 * @param[in] plasma
 *          Plasma context to which this call belong to.
 *
 * @param[in] m
 *         Number of rows of matrix A
 *
 * @param[in] n
 *         Number of columns of matrix A
 *
 * @param[in,out] A
 *         Matrix of size m*n.
 *
 * @param[in] mb
 *         Number of rows of each block
 *
 * @param[in] nb
 *         Number of columns of each block
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 ******************************************************************************/


/** ****************************************************************************
 *
 *  Self-contained functions
 *
 *******************************************************************************/
/*
 * Shift inside panels
 */
int ipt_zcm2ccrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request) 
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    
    PLASMA_zshift(m, n, mb, nb, A);
    ipt_zpanel2tile(plasma, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_zccrb2cm(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request) 
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;

    ipt_ztile2panel(plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_zshiftr(m, n, mb, nb, A);

    return PLASMA_SUCCESS;
}

/*
 * Transpose each tile
 */
int ipt_zccrb2crrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
                   PLASMA_sequence *sequence, PLASMA_request *request) 
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;

    PLASMA_pzgetmi2(idep, odep, PlasmaColumnwise, m, n, mb, nb, A);

    return PLASMA_SUCCESS;
}

int ipt_zcrrb2ccrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
                   PLASMA_sequence *sequence, PLASMA_request *request) 
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;

    PLASMA_pzgetmi2(idep, odep, PlasmaRowwise, n, m, nb, mb, A);

    return PLASMA_SUCCESS;
}

int ipt_zrcrb2rrrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
                   PLASMA_sequence *sequence, PLASMA_request *request) 
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;

    PLASMA_pzgetmi2(idep, odep, PlasmaRowwise, m, n, mb, nb, A);

    return PLASMA_SUCCESS;
}

int ipt_zrrrb2rcrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
                   PLASMA_sequence *sequence, PLASMA_request *request) 
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;

    PLASMA_pzgetmi2(idep, odep, PlasmaColumnwise, n, m, nb, mb, A);

    return PLASMA_SUCCESS;
}

/*
 * Transpose all tiles
 */
int ipt_zccrb2rcrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
                   PLASMA_sequence *sequence, PLASMA_request *request) 
{
    int M_, N_;

    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;

    M_ = m / mb;
    N_ = n / nb;

    /* quick return */
    if( (M_ < 2) || (N_ < 2) ) {
        return PLASMA_SUCCESS;
    }

    plasma_zshift(plasma, m, n, A, 1, ( m / mb ), ( n / nb ), (mb*nb),
                  sequence, request);

    return PLASMA_SUCCESS;
}

/** ****************************************************************************
 *
 *  Composition of 2 sub-routines
 *
 *******************************************************************************/
/*
 * Shift inside panels + Transpose all tiles
 */
int ipt_zcm2rcrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request) 
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    //ipt_zcm2ccrb(  plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_zshift(m, n, mb, nb, A);
    ipt_zpanel2all(plasma, m, n, A, mb, nb, sequence, request);
    ipt_zccrb2rcrb(plasma, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_zrcrb2cm(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request) 
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_zrcrb2ccrb(plasma, m, n, A, mb, nb, sequence, request);
    ipt_zall2panel(plasma, m, n, A, mb, nb, sequence, request);
    //ipt_zccrb2cm(  plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_zshiftr(m, n, mb, nb, A);

    return PLASMA_SUCCESS;
}

/*
 * Transpose each tile  + Transpose all tiles
 */
int ipt_zccrb2rrrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
                   PLASMA_sequence *sequence, PLASMA_request *request)
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_zccrb2rcrb(plasma, m, n, A, mb, nb, sequence, request);
    ipt_zrcrb2rrrb(plasma, PlasmaIPT_All, PlasmaIPT_NoDep, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_zrrrb2ccrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request) 
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_zrrrb2rcrb(plasma, PlasmaIPT_NoDep, PlasmaIPT_All, m, n, A, mb, nb, sequence, request);
    ipt_zrcrb2ccrb(plasma, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_zrcrb2crrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request) 
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_zrcrb2ccrb(plasma, m, n, A, mb, nb, sequence, request);
    ipt_zccrb2crrb(plasma, PlasmaIPT_All, PlasmaIPT_NoDep, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_zcrrb2rcrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
                   PLASMA_sequence *sequence, PLASMA_request *request) 
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_zcrrb2ccrb(plasma, PlasmaIPT_NoDep, PlasmaIPT_All, m, n, A, mb, nb, sequence, request);
    ipt_zccrb2rcrb(plasma, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

/*
 * Transpose each tile  + Shift inside panels
 */
int ipt_zcm2crrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request) 
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    //ipt_zcm2ccrb(  plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_zshift(m, n, mb, nb, A);
    ipt_zccrb2crrb(plasma, PlasmaIPT_Panel, PlasmaIPT_NoDep, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_zcrrb2cm(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request) 
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_zcrrb2ccrb(plasma, PlasmaIPT_NoDep, PlasmaIPT_Panel, m, n, A, mb, nb, sequence, request);
    //ipt_zccrb2cm(  plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_zshiftr(m, n, mb, nb, A);

    return PLASMA_SUCCESS;
}

int ipt_zrm2rcrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request) 
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_zrm2rrrb(  plasma, m, n, A, mb, nb, sequence, request);
    ipt_zrrrb2rcrb(plasma, PlasmaIPT_Panel, PlasmaIPT_NoDep, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_zrcrb2rm(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request) 
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_zrcrb2rrrb(plasma, PlasmaIPT_NoDep, PlasmaIPT_Panel, m, n, A, mb, nb, sequence, request);
    ipt_zrrrb2rm(  plasma, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

/** ****************************************************************************
 *
 *  Composition of 3 sub-routines
 *
 *******************************************************************************/
/*
 * Shift inside panels + Transpose all tiles + Transpose inside each tile
 */
int ipt_zcm2rrrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request) 
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    //ipt_zcm2ccrb(  plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_zshift(m, n, mb, nb, A);
    ipt_zccrb2crrb(plasma, PlasmaIPT_Panel, PlasmaIPT_All, m, n, A, mb, nb, sequence, request);
    ipt_zcrrb2rrrb(plasma, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_zrrrb2cm(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request) 
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_zrrrb2crrb(plasma, m, n, A, mb, nb, sequence, request);
    ipt_zcrrb2ccrb(plasma, PlasmaIPT_All, PlasmaIPT_Panel, m, n, A, mb, nb, sequence, request);
    //ipt_zccrb2cm(  plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_zshiftr(m, n, mb, nb, A);

    return PLASMA_SUCCESS;
}

int ipt_zccrb2rm(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request) 
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_zccrb2rcrb(plasma, m, n, A, mb, nb, sequence, request);
    ipt_zrcrb2rrrb(plasma, PlasmaIPT_All, PlasmaIPT_Panel, m, n, A, mb, nb, sequence, request);
    ipt_zrrrb2rm(  plasma, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_zrm2ccrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request) 
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_zrm2rrrb(  plasma, m, n, A, mb, nb, sequence, request);
    ipt_zrrrb2rcrb(plasma, PlasmaIPT_Panel, PlasmaIPT_All, m, n, A, mb, nb, sequence, request);
    ipt_zrcrb2ccrb(plasma, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

/** ****************************************************************************
 *
 *  Composition of 4 sub-routines
 *
 *******************************************************************************/
/*
 * Shift inside panels + Transpose all tiles 
 *    + Transpose inside each tile + Shift inside panels
 */
int ipt_zcm2rm(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
               PLASMA_sequence *sequence, PLASMA_request *request) 
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    //ipt_zcm2ccrb(  plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_zshift(m, n, mb, nb, A);
    ipt_zpanel2all(plasma, m, n, A, mb, nb, sequence, request);
    ipt_zccrb2rcrb(plasma, m, n, A, mb, nb, sequence, request);
    ipt_zrcrb2rrrb(plasma, PlasmaIPT_All, PlasmaIPT_Panel, m, n, A, mb, nb, sequence, request);
    ipt_zrrrb2rm(  plasma, m, n, A, mb, nb, sequence, request);

    return PLASMA_SUCCESS;
}

int ipt_zrm2cm(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
                 PLASMA_sequence *sequence, PLASMA_request *request) 
{
    if( (m == 0) || (n == 0) )
        return PLASMA_SUCCESS;
    ipt_zrm2rrrb(  plasma, m, n, A, mb, nb, sequence, request);
    ipt_zrrrb2rcrb(plasma, PlasmaIPT_Panel, PlasmaIPT_All, m, n, A, mb, nb, sequence, request);
    ipt_zrcrb2ccrb(plasma, m, n, A, mb, nb, sequence, request);
    ipt_zall2panel(plasma, m, n, A, mb, nb, sequence, request);
    //ipt_zccrb2cm(  plasma, m, n, A, mb, nb, sequence, request);
    PLASMA_zshiftr(m, n, mb, nb, A);

    return PLASMA_SUCCESS;
}


/** ****************************************************************************
 *
 *  Barriers
 *
 *******************************************************************************/

int ipt_ztile2panel( plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
                     PLASMA_sequence *sequence, PLASMA_request *request) 
{
    if( PLASMA_SCHEDULING != PLASMA_DYNAMIC_SCHEDULING )
        return PLASMA_SUCCESS;

    PLASMA_Complex64_t *Al;
    int i,j;
    int M_ = m / mb;
    int N_ = n / nb;
    int bsiz = mb*nb;
    int psiz = m*nb;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    plasma_dynamic_spawn();
    for(j=0; j<N_; j++) {
        Al = &(A[psiz*j]);

        for(i=1; i<M_; i++) {

#ifdef TRACE_IPT
            char str[30];
            sprintf(str, "Foo2 C2RI %d", i*m*nb);
#endif
            QUARK_Insert_Task(plasma->quark, CORE_foo2_quark, &task_flags,
                              sizeof(PLASMA_Complex64_t)*psiz,  Al,           INOUT | GATHERV,
                              sizeof(PLASMA_Complex64_t)*bsiz, &(Al[i*bsiz]), INOUT,
#ifdef TRACE_IPT
                              30, str,   VALUE | TASKLABEL,
                              4, "red",  VALUE | TASKCOLOR,
#endif
                              0);
        }
    }

    return PLASMA_SUCCESS;
}

int ipt_zpanel2tile( plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
                     PLASMA_sequence *sequence, PLASMA_request *request )
{
    if( PLASMA_SCHEDULING != PLASMA_DYNAMIC_SCHEDULING )
        return PLASMA_SUCCESS;

    PLASMA_Complex64_t *Al;
    int i,j;
    int M_ = m / mb;
    int N_ = n / nb;
    int bsiz = mb*nb;
    int psiz = m*nb;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    plasma_dynamic_spawn();
    for(j=0; j<N_; j++) {
        Al = &(A[psiz*j]);

        for(i=1; i<M_; i++) {

#ifdef TRACE_IPT
            char str[30];
            sprintf(str, "Foo2 C2RI %d", i*m*nb);
#endif
            QUARK_Insert_Task(plasma->quark, CORE_foo2_quark, &task_flags,
                              sizeof(PLASMA_Complex64_t)*psiz,  Al,           INPUT,
                              sizeof(PLASMA_Complex64_t)*bsiz, &(Al[i*bsiz]), INOUT,
#ifdef TRACE_IPT
                              30, str,   VALUE | TASKLABEL,
                              4, "red",  VALUE | TASKCOLOR,
#endif
                              0);
        }
    }

    return PLASMA_SUCCESS;
}

int ipt_zpanel2all(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
                   PLASMA_sequence *sequence, PLASMA_request *request) 
{
    if (PLASMA_SCHEDULING != PLASMA_DYNAMIC_SCHEDULING)
        return PLASMA_SUCCESS;

    int i;
    int N_ = n / nb;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    if ( N_ > 1 ) {
        plasma_dynamic_spawn();
        for(i=1; i<N_; i++) {
#ifdef TRACE_IPT
            char str[30];
            sprintf(str, "Foo2 C2RI %d", i*m*nb);
#endif
            QUARK_Insert_Task(plasma->quark, CORE_foo2_quark, &task_flags,
                              sizeof(PLASMA_Complex64_t)*m*n,  A,            INOUT | GATHERV,
                              sizeof(PLASMA_Complex64_t)*m*nb, &(A[i*m*nb]), INPUT,
#ifdef TRACE_IPT
                              30, str,   VALUE | TASKLABEL,
                              4, "red",  VALUE | TASKCOLOR,
#endif
                              0);
        }
    }

    return PLASMA_SUCCESS;
}

int ipt_zall2panel(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb,
                   PLASMA_sequence *sequence, PLASMA_request *request) 
{
    if (PLASMA_SCHEDULING != PLASMA_DYNAMIC_SCHEDULING)
        return PLASMA_SUCCESS;

    int i;
    int N_ = n / nb;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    if ( N_ > 1 ) {
        plasma_dynamic_spawn();
        for(i=1; i<N_; i++) {
#ifdef TRACE_IPT
            char str[30];
            sprintf(str, "Foo2 C2RI %d", i*m*nb);
#endif
            QUARK_Insert_Task(plasma->quark, CORE_foo2_quark, &task_flags,
                              sizeof(PLASMA_Complex64_t)*m*n,  A,            INPUT,
                              sizeof(PLASMA_Complex64_t)*m*nb, &(A[i*m*nb]), INOUT,
#ifdef TRACE_IPT
                              30, str,   VALUE | TASKLABEL,
                              4, "red",  VALUE | TASKCOLOR,
#endif
                              0);
        }
    }
    return PLASMA_SUCCESS;
}
