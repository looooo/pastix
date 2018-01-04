/**
 *
 * @file order_regular_graphs.c
 *
 *  PaStiX order routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * Contains the routine to compute ordering on regular grids/cubes
 *
 * @version 5.1.0
 * @author Gregoire Pichon
 * @date 2016-05-19
 *
 **/
#include "common.h"
#include "spm.h"
#include "graph.h"
#include "pastix/order.h"

void
order_grid2D_classic( pastix_int_t *peritab,
                      pastix_int_t x0,
                      pastix_int_t xn,
                      pastix_int_t y0,
                      pastix_int_t yn,
                      pastix_int_t *max_number,
                      pastix_int_t ldax,
                      pastix_int_t lday )
{
    pastix_int_t dir, i;
    pastix_int_t nx = xn-x0;
    pastix_int_t ny = yn-y0;

    /* The subgraph is small enough */
    if (nx*ny < 50){
        pastix_int_t j;
        pastix_int_t current = 0;
        for (i=0; i<nx; i++){
            for (j=0; j<ny; j++){
                pastix_int_t index = (x0 + i) * ldax + (y0 + j) * lday;
                peritab[index] = max_number[0] - current;
                current++;
            }
        }
        max_number[0] -= current;
        return;
    }

    /* In which direction do we cut? 0 for x, 1 for y */
    dir = 0;
    if (ny > nx)
        dir = 1;

    /* If we cut in direction x */
    if (dir == 0){
        pastix_int_t current = 0;
        for (i=0; i<ny; i++){
            pastix_int_t index = (x0 + nx/2) * ldax + (y0+i) * lday;
            peritab[index] = max_number[0] - current;
            current++;
        }
        max_number[0] -= current;

        order_grid2D_classic(peritab,
                             x0, x0 + nx/2, y0, yn, max_number,
                             ldax, lday);

        order_grid2D_classic(peritab,
                             x0+nx/2+1, xn, y0, yn, max_number,
                             ldax, lday);
    }

    /* If we cut in direction y */
    else if (dir == 1){
        pastix_int_t current = 0;
        for (i=0; i<nx; i++){
            pastix_int_t index = (x0 + i)*ldax + lday * (y0+ny/2);
            peritab[index] = max_number[0] - current;
            current++;
        }
        max_number[0] -= current;

        order_grid2D_classic(peritab,
                             x0, xn, y0, y0+ny/2, max_number,
                             ldax, lday);

        order_grid2D_classic(peritab,
                             x0, xn, y0+ny/2+1, yn, max_number,
                             ldax, lday);
    }
}

void
order_grid3D_classic( pastix_int_t *rangtab,
                      pastix_int_t *peritab,
                      pastix_int_t *cblknbr,
                      pastix_int_t x0,
                      pastix_int_t xn,
                      pastix_int_t y0,
                      pastix_int_t yn,
                      pastix_int_t z0,
                      pastix_int_t zn,
                      pastix_int_t *max_number,
                      pastix_int_t *current_rangtab,
                      pastix_int_t *treetab,
                      pastix_int_t current_treetab,
                      pastix_int_t ldax,
                      pastix_int_t lday,
                      pastix_int_t ldaz )
{
    pastix_int_t dir, i, j;
    pastix_int_t nx = xn-x0;
    pastix_int_t ny = yn-y0;
    pastix_int_t nz = zn-z0;

    /* The subgraph is small enough */
    if (nx*ny*nz < 50){
        pastix_int_t k;
        pastix_int_t current = 0;
        cblknbr[0] ++;
        for (i=x0; i<xn; i++){
            for (j=y0; j<yn; j++){
                for (k=z0; k<zn; k++){
                    pastix_int_t index = i + lday * j + ldaz*ldaz * k;
                    peritab[index] = max_number[0] - current;
                    current++;
                }
            }
        }

        treetab[current_rangtab[0]] = current_treetab;
        rangtab[current_rangtab[0]] = max_number[0];
        max_number[0] -= current;
        current_rangtab[0]++;
        return;
    }

    cblknbr[0] ++;

    /* In which direction do we cut? 0 for x, 1 for y */
    dir = 0;
    if (ny > nx)
        dir = 1;
    if (nz > nx && nz > ny)
        dir = 2;

    /* If we cut in direction x */
    if (dir == 0){
        treetab[current_rangtab[0]] = current_treetab;
        rangtab[current_rangtab[0]] = max_number[0];
        current_rangtab[0]++;

        /* Order separator */
        pastix_int_t *peritab_separator = peritab + x0 + nx/2;
        order_grid2D_classic(peritab_separator,
                             y0, yn, z0, zn,
                             max_number,
                             lday, ldaz*ldaz);

        /* Order nested dissection subparts */
        order_grid3D_classic(rangtab, peritab, cblknbr,
                             x0, x0 + nx/2, y0, yn, z0, zn, max_number,
                             current_rangtab,
                             treetab, current_treetab+1,
                             ldax, lday, ldaz);

        order_grid3D_classic(rangtab, peritab, cblknbr,
                             x0+nx/2+1, xn, y0, yn, z0, zn, max_number,
                             current_rangtab,
                             treetab, current_treetab+1,
                             ldax, lday, ldaz);

    }

    /* If we cut in direction y */
    else if (dir == 1){
        treetab[current_rangtab[0]] = current_treetab;
        rangtab[current_rangtab[0]] = max_number[0];
        current_rangtab[0]++;

        /* Order separator */
        pastix_int_t *peritab_separator = peritab + lday * (y0 + ny / 2);
        order_grid2D_classic(peritab_separator,
                             x0, xn, z0, zn,
                             max_number,
                             1, ldaz*ldaz);

        /* Order nested dissection subparts */
        order_grid3D_classic(rangtab, peritab, cblknbr,
                             x0, xn, y0, y0+ny/2, z0, zn, max_number,
                             current_rangtab,
                             treetab, current_treetab+1,
                             ldax, lday, ldaz);

        order_grid3D_classic(rangtab, peritab, cblknbr,
                             x0, xn, y0+ny/2+1, yn, z0, zn, max_number,
                             current_rangtab,
                             treetab, current_treetab+1,
                             ldax, lday, ldaz);
    }

    /* If we cut in direction z */
    else{

        treetab[current_rangtab[0]] = current_treetab;
        rangtab[current_rangtab[0]] = max_number[0];
        current_rangtab[0]++;

        /* Order separator */
        pastix_int_t *peritab_separator = peritab + ldaz * ldaz * (z0 + nz/2);
        order_grid2D_classic(peritab_separator,
                             x0, xn, y0, yn,
                             max_number,
                             1, lday);

        /* Order nested dissection subparts */
        order_grid3D_classic(rangtab, peritab, cblknbr,
                             x0, xn, y0, yn, z0, z0+nz/2, max_number,
                             current_rangtab,
                             treetab, current_treetab+1,
                             ldax, lday, ldaz);

        order_grid3D_classic(rangtab, peritab, cblknbr,
                             x0, xn, y0, yn, z0+nz/2+1, zn, max_number,
                             current_rangtab,
                             treetab, current_treetab+1,
                             ldax, lday, ldaz);
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_ordering
 *
 * orderComputeOptimal - Compute the ordering of a regular 2D or 3D laplacian,
 * with an optimal strategy.
 *
 *******************************************************************************
 *
 * @param[out] myorder
 *          Personal ordering for regular laplacians
 *
 * @param[in] nx
 *          The number of vertices for the first dimension of the graph
 *
 * @param[in] ny
 *          The number of vertices for the second dimension of the graph
 *
 * @param[in] nz
 *          The number of vertices for the third dimension of the graph
 *
 *******************************************************************************
 *
 * @return
 *          \retval PASTIX_SUCCESS on successful exit
 *
 *******************************************************************************/
int
orderComputeOptimal( pastix_order_t **myorder,
                     pastix_int_t     nx,
                     pastix_int_t     ny,
                     pastix_int_t     nz )
{
    pastix_order_t *ordemesh = *myorder;
    pastix_int_t n = nx * ny * nz;

    pastixOrderAlloc(ordemesh, n, n);

    pastix_int_t *rangtab = ordemesh->rangtab;
    pastix_int_t *permtab = ordemesh->permtab;
    pastix_int_t *peritab = ordemesh->peritab;
    pastix_int_t *treetab = ordemesh->treetab;

    pastix_int_t *saved_rangtab, *saved_treetab;
    pastix_int_t i;

    pastix_int_t current_rangtab = 0;

    /* Graphs for using classical separators */
    if (nx == ny && ny == nz){
        pastix_int_t i = 2;
        while (i != nx && i < nx+1){
            i = 2*i+1;
        }
        if (i != nx){
            printf("\n\nThe given graph size is not correct for optimal manual ordering on 2D regular grid or 3D regular cube. Closer valid sizes are %ld %ld\n\n", i, 2*i+1);
        }
    }

    ordemesh->cblknbr = 0;

    pastix_int_t current_number = n-1;
    order_grid3D_classic(rangtab, permtab, &ordemesh->cblknbr,
                         0, nx, 0, ny, 0, nz, &current_number, &current_rangtab,
                         treetab, 1,
                         nx, ny, nz);

    for (i=0; i<n; i++){
        peritab[permtab[i]] = i;
    }

    saved_rangtab = malloc(n*sizeof(pastix_int_t));
    memcpy(saved_rangtab, rangtab, n*sizeof(pastix_int_t));
    saved_treetab = malloc(n*sizeof(pastix_int_t));
    memcpy(saved_treetab, treetab, n*sizeof(pastix_int_t));

    rangtab[0] = 0;
    for (i=0; i<ordemesh->cblknbr; i++){
        rangtab[i+1] = saved_rangtab[ordemesh->cblknbr - i - 1]+1;
        treetab[i]   = saved_treetab[ordemesh->cblknbr - i - 1];
    }
    for (i=0; i<ordemesh->cblknbr-1; i++){
        pastix_int_t j;
        for (j=i+1; j<ordemesh->cblknbr; j++){
            if (treetab[j] < treetab[i]){
                treetab[i] = j;
                break;
            }
        }
    }
    treetab[ordemesh->cblknbr-1] = -1;

    rangtab =
        (pastix_int_t *) memRealloc (rangtab,
                                     (ordemesh->cblknbr + 1)*sizeof (pastix_int_t));
    treetab =
        (pastix_int_t *) memRealloc (treetab,
                                     (ordemesh->cblknbr)*sizeof (pastix_int_t));

    return PASTIX_SUCCESS;
}

