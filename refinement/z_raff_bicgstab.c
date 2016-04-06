/**
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Theophile Terraz
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "bcsc.h"
#include "z_raff_functions.h"
#include "solver.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_raff
 *
 * z_bicgstab_smp - Function computing bicgstab iterative refinement.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The PaStiX data structure that describes the solver instance.
 *
 * @param[in] x
 *          The solution vector.
 *
 * @param[in] b
 *          The right hand side member (only one).
 *
 *******************************************************************************/
void z_bicgstab_smp (pastix_data_t *pastix_data, void *x, void *b)
{
    /* Choix du solveur */
    struct z_solver solveur;
    z_Pastix_Solveur(&solveur);

    pastix_bcsc_t      * bcsc    = pastix_data->bcsc;
    pastix_int_t         n       = bcsc->gN;
    Clock                raff_clk;
    double               t0      = 0;
    double               t3      = 0;
    pastix_complex64_t   normb   = 0.0;
    pastix_complex64_t   normr   = 0.0;
    int                  nb_iter = 0;
    pastix_complex64_t   epsilon = solveur.Eps(pastix_data);
    pastix_int_t         itermax = solveur.Itermax(pastix_data);
    pastix_complex64_t   tmp     = 0.0;

    pastix_complex64_t * gradb  = NULL; /* Second membre b */
    pastix_complex64_t * gradr  = NULL; /* Solution actuelle */
    pastix_complex64_t * gradr2 = NULL; /* Condition initiale bis r^ */
    pastix_complex64_t * gradp  = NULL;
    pastix_complex64_t * grady  = NULL;
    pastix_complex64_t * gradv  = NULL;
    pastix_complex64_t * grads  = NULL;
    pastix_complex64_t * gradz  = NULL;
    pastix_complex64_t * gradt  = NULL;
    pastix_complex64_t * grad2  = NULL; /* Vecteurs de transition */
    pastix_complex64_t * grad3  = NULL;

    /* Alpha et Beta ne sont utilisés que par le thread 0 */
    pastix_complex64_t * alpha = NULL;
    pastix_complex64_t * beta  = NULL;
    pastix_complex64_t * v1    = NULL;
    pastix_complex64_t * v2    = NULL;
    pastix_complex64_t * w     = NULL;
    pastix_complex64_t * gradx = NULL;

    gradb  = (pastix_complex64_t *)solveur.Malloc(n * sizeof(pastix_complex64_t));
    gradr  = (pastix_complex64_t *)solveur.Malloc(n * sizeof(pastix_complex64_t));
    gradr2 = (pastix_complex64_t *)solveur.Malloc(n * sizeof(pastix_complex64_t));
    gradp  = (pastix_complex64_t *)solveur.Malloc(n * sizeof(pastix_complex64_t));
    grady  = (pastix_complex64_t *)solveur.Malloc(n * sizeof(pastix_complex64_t));
    gradv  = (pastix_complex64_t *)solveur.Malloc(n * sizeof(pastix_complex64_t));
    grads  = (pastix_complex64_t *)solveur.Malloc(n * sizeof(pastix_complex64_t));
    gradz  = (pastix_complex64_t *)solveur.Malloc(n * sizeof(pastix_complex64_t));
    gradt  = (pastix_complex64_t *)solveur.Malloc(n * sizeof(pastix_complex64_t));
    grad2  = (pastix_complex64_t *)solveur.Malloc(n * sizeof(pastix_complex64_t));
    grad3  = (pastix_complex64_t *)solveur.Malloc(n * sizeof(pastix_complex64_t));
    alpha  = (pastix_complex64_t *)solveur.Malloc(    sizeof(pastix_complex64_t));
    beta   = (pastix_complex64_t *)solveur.Malloc(    sizeof(pastix_complex64_t));
    v1     = (pastix_complex64_t *)solveur.Malloc(    sizeof(pastix_complex64_t));
    v2     = (pastix_complex64_t *)solveur.Malloc(    sizeof(pastix_complex64_t));
    w      = (pastix_complex64_t *)solveur.Malloc(    sizeof(pastix_complex64_t));
    gradx  = (pastix_complex64_t *)solveur.Malloc(n * sizeof(pastix_complex64_t));

    clockInit(raff_clk);clockStart(raff_clk);

    solveur.B(b, gradb, n);
    solveur.X(pastix_data, x, gradx);

    /* r = b - Ax */
    solveur.bMAx(bcsc, gradb, gradx, gradr);
    normb = solveur.Norm(gradb, n);
    normr = solveur.Norm(gradr, n);

    /* r2 = r */
    memcpy(gradr2, gradr, n * sizeof( pastix_complex64_t ));
    /* p = r */
    memcpy(gradp, gradr, n * sizeof( pastix_complex64_t ));

    /* tmp = ||r|| / ||b|| */
    tmp = normr / normb;

    while (((double)tmp > (double)epsilon) && (nb_iter < itermax))
    {
        clockStop((raff_clk));
        t0 = clockGet();
        nb_iter++;

        /* y = M-1 * p */
        solveur.Precond(pastix_data, gradp, grady);

        /* v = Ay */
        solveur.Ax(bcsc, grady, gradv);

        /* alpha = (r, r2) / (v, r2) */
        /* alpha = (v, r2) */
        solveur.Dotc(n, gradv, gradr2, alpha);
        /* beta = (r, r2) */
        solveur.Dotc(n, gradr, gradr2, beta);

        /* alpha = beta / alpha : alpha = (r, r2) / (v, r2) */
        //       solveur.Div(arg, beta, alpha, alpha, 0);
        alpha[0] = beta[0] / alpha[0];

        /* s = r - alpha * v */
        memcpy(grads, gradr, n * sizeof( pastix_complex64_t ));
        solveur.AXPY(n, -1, alpha, grads, gradv);

        /* z = M-1s */
        solveur.Precond(pastix_data, grads, gradz);

        /* t = Az */
        solveur.Ax(bcsc, gradz, gradt);

        /* w = (M-1t, M-1s) / (M-1t, M-1t) */
        /* grad2 = M-1t */
        solveur.Precond(pastix_data, gradt, grad2);

        /* v1 = (M-1t, M-1s) */
        /* v2 = (M-1t, M-1t) */
        solveur.Dotc(n, gradz, grad2, v1);
        solveur.Dotc(n, grad2, grad2, v2);

        //       solveur.Div(arg, v1, v2, w, 1);
        w[0] = v1[0] / v2[0];

        /* x = x + alpha * y + w * z */
        /* x = x + alpha * y */
        solveur.AXPY(n, 1, alpha, gradx, grady);

        /* x = x + w * z */
        solveur.AXPY(n, 1, w, gradx, gradz);

        /* r = s - w * t*/
        memcpy(gradr, grads, n * sizeof( pastix_complex64_t ));
        solveur.AXPY(n, -1, w, gradr, gradt);

        /* beta = (r', r2) / (r, r2) * (alpha / w) */
        /* v1 = (r', r2) */
        solveur.Dotc(n, gradr, gradr2, v1);

        /* v2 = alpha / w */
        //       solveur.Div(arg, alpha, w, v2, 0);
        v2[0] = alpha[0] / w[0];

        /* beta = v1 / beta */
        //       solveur.Div(arg, v1, beta, beta, 0);
        beta[0] = v1[0] / beta[0];

        /* beta = beta * v2 */
        //       solveur.Mult(arg, beta, v2, beta, 1);
        beta[0] = beta[0] * v2[0];

        /* p = r + beta * (p - w * v) */
        /* p = p - w * v */
        solveur.AXPY(n, -1, w, gradp, gradv);

        /* p = r + beta * p */
        solveur.BYPX(n, beta, gradr, gradp);

        normr = solveur.Norm(gradr, n);

        clockStop((raff_clk));
        t3 = clockGet();

        tmp = normr / normb;
        solveur.Verbose(t0, t3, tmp, nb_iter);
    }

    solveur.End(pastix_data, tmp, nb_iter, t3, x, gradx);

    solveur.Free((void*) gradb);
    solveur.Free((void*) gradr);
    solveur.Free((void*) gradr2);
    solveur.Free((void*) gradp);
    solveur.Free((void*) grady);
    solveur.Free((void*) gradv);
    solveur.Free((void*) grads);
    solveur.Free((void*) gradz);
    solveur.Free((void*) gradt);
    solveur.Free((void*) grad2);
    solveur.Free((void*) grad3);
    solveur.Free((void*) alpha);
    solveur.Free((void*) beta);
    solveur.Free((void*) v1);
    solveur.Free((void*) v2);
    solveur.Free((void*) w);
    solveur.Free((void*) gradx);
}
