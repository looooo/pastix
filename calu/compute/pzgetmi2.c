/**
 *
 * @file pzgetmi2.c
 *
 *  PLASMA InPlaceTransformation module
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 *  This work is the implementation of an inplace transformation
 *  based on the GKK algorithm by Gustavson, Karlsson, Kagstrom
 *  and its fortran implementation.
 *
 * @version 2.5.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

/** ****************************************************************************
 *
 * @ingroup InPlaceTransformation
 *
 * plasma_pzgetmi2 - realises nprob independant transpositions. Each
 * subproblem is a tile of mb-by-nb elements.
 * This function use an extra space of PLASMA_SIZE*(mb*nb).
 *
 *******************************************************************************
 *
 * @param[in] plasma
 *          Plasma context to which this call belong to.
 *
 *******************************************************************************
 *
 * @sa plasma_pzgetmi2_quark
 *
 ******************************************************************************/
void plasma_pzgetmi2(plasma_context_t *plasma) {
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    PLASMA_Complex64_t *A, *Al, *work;
    PLASMA_enum storev, idep, odep;
    int         i, m, n, mb, nb, nprob;
    int         size, bsiz;

    plasma_unpack_args_10(idep, odep, storev, m, n, mb, nb, A, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;

    /* quick return */
    if( (mb < 2) || (nb < 2) ) {
        return ;
    }

    size  = PLASMA_SIZE;
    bsiz  = mb*nb;
    nprob = ( m / mb ) * ( n / nb );

    work = (PLASMA_Complex64_t*)plasma_private_alloc(plasma, mb*nb, PlasmaComplexDouble);

    for (i=PLASMA_RANK; i<nprob; i+=size) {
        Al = &(A[ i * bsiz]);
        CORE_zgetrip(mb, nb, Al, work);
    }

    plasma_private_free(plasma, work);
}

/** ****************************************************************************
 *
 * @ingroup InPlaceTransformation
 *
 * plasma_pzgetmi2_quark - realises nprob independant transpositions. Each
 * subproblem is a tile of mb-by-nb elements.
 * This function use an extra space of PLASMA_SIZE*(mb*nb). This is a
 * maximum in case of dynamic scheduling.
 * 
 *
 *******************************************************************************
 *
 * @param[in] idep
 *       PlasmaIPT_Nodep: No fake dependencies are added.
 *       PlasmaIPT_Panel: A gatherv is added on each panel and panel size is m*nb.
 *       PlasmaIPT_All:   A gatherv is added on the whole matrix.
 *
 * @param[in] odep
 *       PlasmaIPT_Nodep: No fake dependencies are added.
 *       PlasmaIPT_Panel: A gatherv is added on each panel and panel size is m*nb.
 *       PlasmaIPT_All:   A gatherv is added on the whole matrix.
 *
 * @param[in] storev
 *       PlasmaColumnWise: Data stored in column major.
 *       PlasmaRowWise: Data stored in row major.
 *
 * @param[in] m
 *       Number of row of A if tiles are sorted in column major format,
 *       number of columns otherwise.
 *
 * @param[in] n
 *       Number of columns of A if tiles are sorted in column major format,
 *       number of rows otherwise.
 *
 * @param[in] mb
 *       Number of rows in each individual subproblem if storev == PlasmaColumnWise, 
 *       number of columns otherwise. m%mb must be 0.
 *
 * @param[in] nb
 *       Number of columns in each individual subproblem if storev == PlasmaColumnWise, 
 *       number of rows otherwise. n%nb must be 0.
 *
 * @param[in,out] A
 *       Matrix of size m*n.
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa plasma_pzgetmi2
 *
 ******************************************************************************/
void plasma_pzgetmi2_quark(PLASMA_enum idep, PLASMA_enum odep, PLASMA_enum storev, 
                           int m, int n, int mb, int nb, PLASMA_Complex64_t *A,
                           PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    PLASMA_Complex64_t *Al, *Ap;
    int i, j, nprob, mt, nt;
    int bsiz, psiz, size;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;

    /* quick return */
    if( (mb < 2) || (nb < 2) ) {
        return ;
    }

    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    bsiz  = mb*nb;
    if ( storev == PlasmaColumnwise ) {
        psiz = m*nb;
        mt = ( m / mb );
        nt = ( n / nb );
    } else {
        psiz = n*mb;
        mt = ( n / nb );
        nt = ( m / mb );    
    }

    size = m*n;

    switch ( idep ) {
/*
 * Dependencies on each panel as input
 */
    case PlasmaIPT_Panel:
        switch ( odep ) {
        case PlasmaIPT_Panel:
            for (j=0; j<nt; j++) {
                Ap = A + (psiz*j);
                for (i=0; i<mt; i++) {
                    Al = Ap + i*bsiz;
                    QUARK_CORE_zgetrip_f1(
                        plasma->quark, &task_flags,
                        mb, nb, Al, bsiz, 
                        Ap, psiz, INOUT|GATHERV);
                }
            }
            break;
              
        case PlasmaIPT_All:
            for (j=0; j<nt; j++) {
                Ap = A + (psiz*j);
                for (i=0; i<mt; i++) {
                    Al = Ap + i*bsiz;
                    
                    QUARK_CORE_zgetrip_f2(plasma->quark, &task_flags,
                                          mb, nb, Al, bsiz, 
                                          Ap, size, INPUT,
                                          A,  size, INOUT|GATHERV);
                }
            }
            break;

        case PlasmaIPT_NoDep:
        default:
            for (j=0; j<nt; j++) {
                Ap = A + (psiz*j);
                for (i=0; i<mt; i++) {
                    Al = Ap + i*bsiz;
                    QUARK_CORE_zgetrip_f1(
                        plasma->quark, &task_flags,
                        mb, nb, Al, bsiz, 
                        Ap, psiz, INPUT);
                }
            }
        }
        break;
        
/*
 * Dependency on all the matrix as input
 */
    case PlasmaIPT_All:
        switch ( odep ) {
        case PlasmaIPT_Panel:
            for (j=0; j<nt; j++) {
                Ap = A + (psiz*j);
                for (i=0; i<mt; i++) {
                    Al = Ap + i*bsiz;
                    QUARK_CORE_zgetrip_f2(
                        plasma->quark, &task_flags,
                        mb, nb, Al, bsiz, 
                        A,  size, INPUT,
                        Ap, psiz, INOUT|GATHERV);
                }
            }
            break;
              
        case PlasmaIPT_All:
            nprob = mt*nt;
            for (i=0; i<nprob; i++) {
                QUARK_CORE_zgetrip_f1(plasma->quark, &task_flags,
                                      mb, nb, &(A[ i*bsiz ]), bsiz, 
                                      A, size, INOUT|GATHERV);
            }
            break;
            
        case PlasmaIPT_NoDep:
        default:
            nprob = mt*nt;
            for (i=0; i<nprob; i++) {
                QUARK_CORE_zgetrip_f1(plasma->quark, &task_flags,
                                      mb, nb, &(A[ i*bsiz ]), bsiz, 
                                      A, size, INPUT);
            }
        }
        break;
        
/*
 * No Dependencies as input
 */
    case PlasmaIPT_NoDep:
    default:
        switch ( odep ) {
        case PlasmaIPT_Panel:
            for (j=0; j<nt; j++) {
                Ap = A + (psiz*j);
                for (i=0; i<mt; i++) {
                    Al = Ap + i*bsiz;
                    QUARK_CORE_zgetrip_f1(
                        plasma->quark, &task_flags,
                        mb, nb, Al, bsiz, 
                        Ap, psiz, INOUT|GATHERV);
                }
            }
            break;
              
        case PlasmaIPT_All:
            nprob = mt*nt;
            for (i=0; i<nprob; i++) {
                QUARK_CORE_zgetrip_f1(plasma->quark, &task_flags,
                                      mb, nb, &(A[ i*bsiz ]), bsiz, 
                                      A, size, INOUT|GATHERV);
            }
            break;

        case PlasmaIPT_NoDep:
        default:
            nprob = mt*nt;
            for (i=0; i<nprob; i++) {
                QUARK_CORE_zgetrip(plasma->quark, &task_flags,
                                   mb, nb,  &(A[ i*bsiz ]), bsiz);
            }
        }
    }
}
