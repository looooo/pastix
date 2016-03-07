/**
 *
 * @file: z_raff_functions.c
 *
 *  Functions computing operations for refinement methods
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Theophile Terraz
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "z_spm.h"
#include "bcsc.h"
#include "z_bcsc.h"
#include "sopalin_thread.h"
#include "sopalin_data.h"
#include "solver.h"
#include "z_raff_functions.h"

/* Alloue un vecteur de taille size octets */
void *z_Pastix_Malloc(size_t size)
{
    void *x = NULL;
    MALLOC_INTERN(x, size, char);
    memset(x, 0, size);
    return x;
}

/* Libere un vecteur */
void z_Pastix_Free( void *x)
{
    memFree_null(x);
}


/*** GESTION DE L'INTERFACE ***/

/* Affichage à chaque itération et communication de certaines informations à la structure */
void z_Pastix_Verbose(double t0, double t3, double tmp, pastix_int_t nb_iter)
{
    double stt;

    stt = t3 - t0;
    fprintf(stdout, OUT_ITERRAFF_ITER, (int)nb_iter);
    fprintf(stdout, OUT_ITERRAFF_TTT, stt);
    fprintf(stdout, OUT_ITERRAFF_ERR, tmp);
}

/* Affichage final */
void z_Pastix_End(pastix_data_t *pastix_data, pastix_complex64_t tmp, pastix_int_t nb_iter, double t, void *x, pastix_complex64_t *gmresx)
{
    (void) tmp;
    (void) nb_iter;
    (void) t;
    pastix_complex64_t *xptr = (pastix_complex64_t *)x;
    pastix_int_t        n = pastix_data->bcsc->gN;
    pastix_int_t i;

    for (i=0; i<n; i++)
        xptr[i] = gmresx[i];


    //   if (sopar->iparm[IPARM_PRODUCE_STATS] == API_YES) {
    //     PASTIX_FLOAT *r, *s;
    //
    //     MONOTHREAD_BEGIN;
    //     MALLOC_INTERN(r, UPDOWN_SM2XSZE, PASTIX_FLOAT);
    //     MALLOC_INTERN(s, UPDOWN_SM2XSZE, PASTIX_FLOAT);
    //     sopalin_data->ptr_raff[0] = (void *)r;
    //     sopalin_data->ptr_raff[1] = (void *)s;
    //     MONOTHREAD_END;
    //     SYNCHRO_THREAD;
    //
    //     r = (PASTIX_FLOAT *)sopalin_data->ptr_raff[0];
    //     s = (PASTIX_FLOAT *)sopalin_data->ptr_raff[1];
    //     MULTITHREAD_BEGIN;
    //     /* compute r = b - Ax */
    //     CscbMAx(sopalin_data, me, r, sopar->b, sopar->cscmtx,
    //             &(datacode->updovct), datacode, PASTIX_COMM,
    //             sopar->iparm[IPARM_TRANSPOSE_SOLVE]);
    //     /* |A||x| + |b| */
    //     CscAxPb( sopalin_data, me, s, sopar->b, sopar->cscmtx,
    //              &(datacode->updovct), datacode, PASTIX_COMM,
    //              sopar->iparm[IPARM_TRANSPOSE_SOLVE]);
    //     CscBerr(sopalin_data, me, r, s, UPDOWN_SM2XSZE,
    //             1, &(sopalin_data->sopar->dparm[DPARM_SCALED_RESIDUAL]),
    //             PASTIX_COMM);
    //     MULTITHREAD_END(1);
    //   }
}

/* Vecteur solution X */
void z_Pastix_X(pastix_data_t *pastix_data, void *x, pastix_complex64_t *gmresx)
{
    pastix_int_t        i;
    pastix_int_t        n = pastix_data->bcsc->gN;
    pastix_complex64_t *xptr = (pastix_complex64_t *)x;

    if (pastix_data->iparm[IPARM_ONLY_RAFF] == API_NO)
    {
        for (i=0; i<n; i++, xptr++)
            gmresx[i]= *xptr;
    }
    else
    {
        for (i=0; i<n; i++, xptr++)
            gmresx[i]=0.0;
    }
}

/* Taille d'un vecteur */
pastix_int_t z_Pastix_n(pastix_data_t *pastix_data)
{
    return pastix_data->bcsc->gN;
}

/* Second membre */
void z_Pastix_B(void *b, pastix_complex64_t *raffb, pastix_int_t n)
{
    pastix_complex64_t *bptr = (pastix_complex64_t *)b;
    pastix_int_t i;

    for (i=0; i<n; i++, bptr++)
    {
        raffb[i]= *bptr;
    }
}

/* Epsilon */
pastix_complex64_t z_Pastix_Eps(pastix_data_t *pastix_data)
{
    return pastix_data->dparm[DPARM_EPSILON_REFINEMENT];
}

/* Itermax */
pastix_int_t z_Pastix_Itermax(pastix_data_t *pastix_data)
{
    return pastix_data->iparm[IPARM_ITERMAX];
}


/* Itermax */
pastix_int_t z_Pastix_Krylov_Space(pastix_data_t *pastix_data)
{
    return pastix_data->iparm[IPARM_GMRES_IM];
}

/*** OPERATIONS DE BASE ***/
/* Multiplication pour plusieurs second membres */
// void z_Pastix_Mult(void *arg, pastix_complex64_t *alpha, pastix_complex64_t *beta, pastix_complex64_t *zeta, int flag)
// {
//   sopthread_data_t *argument     = (sopthread_data_t *)arg;
//   sopalin_data_t   *sopalin_data = (sopalin_data_t *)(argument->data);
//   SolverMatrix     *datacode     = sopalin_data->datacode;
//   pastix_int_t        me           = argument->me;
// //   MONOTHREAD_BEGIN;
// #ifdef MULT_SMX_RAFF
//   {
//     pastix_int_t itersmx;
//     for(itersmx=0; itersmx<UPDOWN_SM2XNBR;itersmx++)
//       {
//         zeta[itersmx]=alpha[itersmx]*beta[itersmx];
//       }
//   }
// #else
//   zeta[0]=alpha[0]*beta[0];
// #endif
// //   MONOTHREAD_END;
// //   if (flag)
// //     SYNCHRO_THREAD;
// }

/* Division pour plusieurs second membres */
// void z_Pastix_Div(void *arg, pastix_complex64_t *alpha, pastix_complex64_t *beta, pastix_complex64_t *zeta, int flag)
// {
//   sopthread_data_t *argument     = (sopthread_data_t *)arg;
//   sopalin_data_t   *sopalin_data = (sopalin_data_t *)(argument->data);
//   SolverMatrix     *datacode     = sopalin_data->datacode;
//   pastix_int_t        me           = argument->me;
// //   MONOTHREAD_BEGIN;
// #ifdef MULT_SMX_RAFF
//   {
//     pastix_int_t itersmx;
//     for(itersmx=0; itersmx<UPDOWN_SM2XNBR;itersmx++)
//       {
//         zeta[itersmx]=alpha[itersmx]/beta[itersmx];
//       }
//   }
// #else
//   zeta[0]=alpha[0]/beta[0];
// #endif
// //   MONOTHREAD_END;
// //   if (flag)
// //     SYNCHRO_THREAD;
// }

/* Calcul de la norme de frobenius */
pastix_complex64_t z_Pastix_Norm2(pastix_complex64_t *x, pastix_int_t n)
{
    double normx;
    void *xptr = (void*)x;
    normx = z_vectFrobeniusNorm(xptr, n);
    return normx;
}

/* Copie d'un vecteur */
// void z_Pastix_Copy(void *arg, pastix_complex64_t *s, pastix_complex64_t *d, int flag)
// {
//   sopthread_data_t *argument     = (sopthread_data_t *)arg;
//   sopalin_data_t   *sopalin_data = (sopalin_data_t *)(argument->data);
//   SolverMatrix     *datacode     = sopalin_data->datacode;
//   MPI_Comm          pastix_comm  = PASTIX_COMM;
//   pastix_int_t        me           = argument->me;
// //   MULTITHREAD_BEGIN;
//   z_CscCopy(sopalin_data, me, s, d,
//           UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, pastix_comm);
// //   MULTITHREAD_END(0);
//
// //   if (flag)
// //     SYNCHRO_THREAD;
// }

/* Application du préconditionneur */
void z_Pastix_Precond(pastix_data_t *pastix_data, pastix_complex64_t *s, pastix_complex64_t *d)
{
    pastix_int_t n = pastix_data->bcsc->gN;
    pastix_int_t nrhs = 1;
    void* bptr = (void*)d;

    memcpy(d, s, n * sizeof( pastix_complex64_t ));
    if (pastix_data->iparm[IPARM_ONLY_RAFF] == API_NO)
    {
        sopalin_data_t sopalin_data;
        sopalin_data.solvmtx = pastix_data->solvmatr;

        switch ( pastix_data->iparm[IPARM_FACTORIZATION] ){
        case PastixFactLLT:
            sequential_ztrsm( pastix_data, PastixLeft, PastixLower, PastixNoTrans,   PastixNonUnit, &sopalin_data, nrhs, bptr, n );
            sequential_ztrsm( pastix_data, PastixLeft, PastixLower, PastixConjTrans, PastixNonUnit, &sopalin_data, nrhs, bptr, n );
            break;

        case PastixFactLDLT:
            sequential_ztrsm( pastix_data, PastixLeft, PastixLower, PastixNoTrans, PastixUnit, &sopalin_data, nrhs, bptr, n );
            sequential_zdiag( pastix_data, &sopalin_data, nrhs, bptr, n );
            sequential_ztrsm( pastix_data, PastixLeft, PastixLower, PastixTrans,   PastixUnit, &sopalin_data, nrhs, bptr, n );
            break;

        case PastixFactLDLH:
            sequential_ztrsm( pastix_data, PastixLeft, PastixLower, PastixNoTrans,   PastixUnit, &sopalin_data, nrhs, bptr, n );
            sequential_zdiag( pastix_data, &sopalin_data, nrhs, bptr, n );
            sequential_ztrsm( pastix_data, PastixLeft, PastixLower, PastixConjTrans, PastixUnit, &sopalin_data, nrhs, bptr, n );
            break;

        case PastixFactLU:
        default:
            sequential_ztrsm( pastix_data, PastixLeft, PastixLower, PastixNoTrans, PastixUnit,    &sopalin_data, nrhs, bptr, n );
            sequential_ztrsm( pastix_data, PastixLeft, PastixUpper, PastixNoTrans, PastixNonUnit, &sopalin_data, nrhs, bptr, n );
            break;
        }
    }
}

/* Calcul de alpha * x */
void z_Pastix_Scal(pastix_int_t n, pastix_complex64_t alpha, pastix_complex64_t *x)
{
    z_bcscScal( x, alpha, n, 1);
}

/* Calcul du produit scalaire */
#if defined(PRECISION_z) || defined(PRECISION_c)
void z_Pastix_Dotc(pastix_int_t n, pastix_complex64_t *x, pastix_complex64_t *y, pastix_complex64_t *r)
{
    *r = z_bcscDotc(n, x, y);
}
#endif

void z_Pastix_Dotu(pastix_int_t n, pastix_complex64_t *x, pastix_complex64_t *y, pastix_complex64_t *r)
{
    *r = z_bcscDotu(n, x, y);
}

/* Produit matrice vecteur */
void z_Pastix_Ax(pastix_bcsc_t *bcsc, pastix_complex64_t *x, pastix_complex64_t *r)
{
    pastix_int_t alpha = 1.0;
    pastix_int_t beta = 0.0;
    void* xptr = (void*)x;
    void* yptr = (void*)r;

    z_bcscGemv(PastixNoTrans, alpha, bcsc, xptr, beta, yptr );
}

void z_Pastix_bMAx(pastix_bcsc_t *bcsc, pastix_complex64_t *b, pastix_complex64_t *x, pastix_complex64_t *r)
{
    pastix_int_t alpha = -1.0;
    pastix_int_t beta = 1.0;
    void* xptr = (void*)x;
    void* yptr = (void*)r;

    memcpy(r, b, bcsc->gN * sizeof( pastix_complex64_t ));
    z_bcscGemv(PastixNoTrans, alpha, bcsc, xptr, beta, yptr );
}

/* x = y + beta * x */
void z_Pastix_BYPX(pastix_int_t n, pastix_complex64_t *beta, pastix_complex64_t *y, pastix_complex64_t *x)
{
    void *yptr = (void*)y;
    void *xptr = (void*)x;

    z_bcscScal( xptr, *beta, n, 1);
    z_bcscAxpy( n, 1, 1., yptr, xptr );
}


void z_Pastix_AXPY(pastix_int_t n, double coeff, pastix_complex64_t *alpha, pastix_complex64_t *x, pastix_complex64_t *y)
{
    void *yptr = (void*)y;
    void *xptr = (void*)x;
    z_bcscAxpy( n, 1, coeff*(*alpha), yptr, xptr );
}


pastix_int_t z_Pastix_me(void *arg)
{
    sopthread_data_t *argument = (sopthread_data_t *)arg;
    pastix_int_t        me       = argument->me;
    return me;
}

void z_Pastix_Solveur(struct z_solver *solveur)
{
    /*** ALLOCATIONS ET SYNCHRONISATIONS ***/
    solveur->Malloc      = &z_Pastix_Malloc;
    solveur->Free        = &z_Pastix_Free;

    /*** GESTION DE L'INTERFACE ***/
    solveur->Verbose = &z_Pastix_Verbose;
    solveur->End     = &z_Pastix_End;
    solveur->X       = &z_Pastix_X;
    solveur->N       = &z_Pastix_n;
    solveur->B       = &z_Pastix_B;
    solveur->Eps     = &z_Pastix_Eps;
    solveur->Itermax = &z_Pastix_Itermax;
    solveur->me      = &z_Pastix_me;
    solveur->Krylov_Space = &z_Pastix_Krylov_Space;

    /*** OPERATIONS DE BASE ***/

    solveur->Norm    = &z_Pastix_Norm2;
    solveur->Precond = &z_Pastix_Precond;

    solveur->Scal    = &z_Pastix_Scal;
    solveur->Dotc    = &z_Pastix_Dotc;
    solveur->Ax      = &z_Pastix_Ax;

    solveur->AXPY    = &z_Pastix_AXPY;
    solveur->bMAx    = &z_Pastix_bMAx;
    solveur->BYPX    = &z_Pastix_BYPX;
}
