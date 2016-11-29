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
#include "order.h"

void
order_grid2D_wide( pastix_int_t *rangtab,
                   pastix_int_t *peritab,
                   pastix_int_t *cblknbr,
                   pastix_int_t x0,
                   pastix_int_t xn,
                   pastix_int_t y0,
                   pastix_int_t yn,
                   pastix_int_t max_number,
                   pastix_int_t lda,
                   pastix_int_t *current_rangtab )
{
    pastix_int_t dir, i;
    pastix_int_t nx = xn-x0;
    pastix_int_t ny = yn-y0;

    /* The subgraph is small enough */
    if (nx <= 4 && ny <= 4){
        pastix_int_t j;
        pastix_int_t current = 0;
        cblknbr[0] ++;
        for (i=x0; i<xn; i++){
            for (j=y0; j<yn; j++){
                pastix_int_t index = i + lda * j;
                peritab[index] = max_number - current;
                current++;
            }
        }

        rangtab[current_rangtab[0]] = max_number;
        current_rangtab[0]++;
        return;
    }

    cblknbr[0] ++;

    /* In which direction do we cut? 0 for x, 1 for y */
    dir = 0;
    if (ny > nx)
        dir = 1;

    /* If we cut in direction x */
    if (dir == 0){
        rangtab[current_rangtab[0]] = max_number;
        current_rangtab[0]++;

        for (i=0; i<ny; i++){
            pastix_int_t index = x0 + nx/2 + lda*(y0+i) - 1;
            peritab[index] = max_number - i;
        }
        for (i=0; i<ny; i++){
            pastix_int_t index = x0 + nx/2 + lda*(y0+i);
            peritab[index] = max_number - ny - i;
        }

        order_grid2D_wide(rangtab, peritab, cblknbr,
                          x0, xn-nx/2-1, y0, yn, max_number - 2*ny,
                          lda, current_rangtab);

        order_grid2D_wide(rangtab, peritab, cblknbr,
                          x0+nx/2+1, xn, y0, yn, max_number - 2*ny - ny*(nx/2-1),
                     lda, current_rangtab);
    }

    /* If we cut in direction y */
    else{
        rangtab[current_rangtab[0]] = max_number;
        current_rangtab[0]++;

        for (i=0; i<nx; i++){
            pastix_int_t index = lda*(y0+ny/2-1) + x0 + i;
            peritab[index] = max_number - i;
        }
        for (i=0; i<nx; i++){
            pastix_int_t index = lda*(y0+ny/2) + x0 + i;
            peritab[index] = max_number - nx - i;
        }

        order_grid2D_wide(rangtab, peritab, cblknbr,
                          x0, xn, y0, yn-ny/2-1, max_number - 2*nx,
                          lda, current_rangtab);

        order_grid2D_wide(rangtab, peritab, cblknbr,
                          x0, xn, y0+ny/2+1, yn, max_number - 2*nx - nx*(ny/2-1),
                          lda, current_rangtab);
    }
}

void
order_grid3D_wide( pastix_int_t *rangtab,
                   pastix_int_t *peritab,
                   pastix_int_t *cblknbr,
                   pastix_int_t x0,
                   pastix_int_t xn,
                   pastix_int_t y0,
                   pastix_int_t yn,
                   pastix_int_t z0,
                   pastix_int_t zn,
                   pastix_int_t *max_number,
                   pastix_int_t lda,
                   pastix_int_t *current_rangtab,
                   pastix_int_t *treetab,
                   pastix_int_t current_treetab )
{
    pastix_int_t dir, i, j;
    pastix_int_t nx = xn-x0;
    pastix_int_t ny = yn-y0;
    pastix_int_t nz = zn-z0;

    /* The subgraph is small enough */
    if (nx*ny*nz < 64){
        cblknbr[0] ++;
        pastix_int_t k;
        pastix_int_t current = 0;
        for (i=x0; i<xn; i++){
            for (j=y0; j<yn; j++){
                for (k=z0; k<zn; k++){
                    pastix_int_t index = i + lda * j + lda*lda * k;
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
        pastix_int_t current = 0;
        treetab[current_rangtab[0]] = current_treetab;
        rangtab[current_rangtab[0]] = max_number[0];
        current_rangtab[0]++;

        for (i=0; i<ny; i++){
            for (j=0; j<nz; j++){
                pastix_int_t index = x0 + nx/2 + lda*(y0+i) - 1 + lda*lda*(z0+j);
                peritab[index] = max_number[0] - current;
                current++;
            }
        }
        for (i=0; i<ny; i++){
            for (j=0; j<nz; j++){
                pastix_int_t index = x0 + nx/2 + lda*(y0+i) + lda*lda*(z0+j);
                peritab[index] = max_number[0] - current;
                current++;
            }
        }
        max_number[0] -= current;

        order_grid3D_wide(rangtab, peritab, cblknbr,
                          x0, x0 + nx/2 - 1, y0, yn, z0, zn, max_number,
                          lda, current_rangtab,
                          treetab, current_treetab+1);

        order_grid3D_wide(rangtab, peritab, cblknbr,
                          x0+nx/2+1, xn, y0, yn, z0, zn, max_number,
                          lda, current_rangtab,
                          treetab, current_treetab+1);

    }

    /* If we cut in direction y */
    else if (dir == 1){
        pastix_int_t current = 0;
        treetab[current_rangtab[0]] = current_treetab;
        rangtab[current_rangtab[0]] = max_number[0];
        current_rangtab[0]++;

        for (i=0; i<nx; i++){
            for (j=0; j<nz; j++){
                pastix_int_t index = lda*(y0+ny/2-1) + x0 + i  + lda*lda*(z0+j);
                peritab[index] = max_number[0] - current;
                current++;
            }
        }
        for (i=0; i<nx; i++){
            for (j=0; j<nz; j++){
                pastix_int_t index = lda*(y0+ny/2) + x0 + i + lda*lda*(z0+j);
                peritab[index] = max_number[0] - current;
                current++;
            }
        }
        max_number[0] -= current;

        order_grid3D_wide(rangtab, peritab, cblknbr,
                          x0, xn, y0, y0+ny/2-1, z0, zn, max_number,
                          lda, current_rangtab,
                          treetab, current_treetab+1);

        order_grid3D_wide(rangtab, peritab, cblknbr,
                          x0, xn, y0+ny/2+1, yn, z0, zn, max_number,
                          lda, current_rangtab,
                          treetab, current_treetab+1);
    }

    /* If we cut in direction z */
    else{
        pastix_int_t current = 0;
        treetab[current_rangtab[0]] = current_treetab;
        rangtab[current_rangtab[0]] = max_number[0];
        current_rangtab[0]++;

        for (i=0; i<nx; i++){
            for (j=0; j<ny; j++){
                pastix_int_t index = lda*lda*(z0+nz/2-1) + x0 + i + lda * (y0+j);
                peritab[index] = max_number[0] - current;
                current++;
            }
        }
        for (i=0; i<nx; i++){
            for (j=0; j<ny; j++){
                pastix_int_t index = lda*lda*(z0+nz/2) + x0 + i + lda * (y0+j);
                peritab[index] = max_number[0] - current;
                current++;
            }
        }
        max_number[0] -= current;

        order_grid3D_wide(rangtab, peritab, cblknbr,
                          x0, xn, y0, yn, z0, z0+nz/2-1, max_number,
                          lda, current_rangtab,
                          treetab, current_treetab+1);

        order_grid3D_wide(rangtab, peritab, cblknbr,
                          x0, xn, y0, yn, z0+nz/2+1, zn, max_number,
                          lda, current_rangtab,
                          treetab, current_treetab+1);
    }
}

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
                      pastix_int_t lda,
                      pastix_int_t *current_rangtab,
                      pastix_int_t *treetab,
                      pastix_int_t current_treetab )
{
    pastix_int_t dir, i, j;
    pastix_int_t nx = xn-x0;
    pastix_int_t ny = yn-y0;
    pastix_int_t nz = zn-z0;

    /* The subgraph is small enough */
    if (nx*ny*nz < 10){
        pastix_int_t k;
        pastix_int_t current = 0;
        cblknbr[0] ++;
        for (i=x0; i<xn; i++){
            for (j=y0; j<yn; j++){
                for (k=z0; k<zn; k++){
                    pastix_int_t index = i + lda * j + lda*lda * k;
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
                             lda, lda*lda);

        /* Order nested dissection subparts */
        order_grid3D_classic(rangtab, peritab, cblknbr,
                             x0, x0 + nx/2, y0, yn, z0, zn, max_number,
                             lda, current_rangtab,
                             treetab, current_treetab+1);

        order_grid3D_classic(rangtab, peritab, cblknbr,
                             x0+nx/2+1, xn, y0, yn, z0, zn, max_number,
                             lda, current_rangtab,
                             treetab, current_treetab+1);

    }

    /* If we cut in direction y */
    else if (dir == 1){
        treetab[current_rangtab[0]] = current_treetab;
        rangtab[current_rangtab[0]] = max_number[0];
        current_rangtab[0]++;

        /* Order separator */
        pastix_int_t *peritab_separator = peritab + lda * (y0 + ny / 2);
        order_grid2D_classic(peritab_separator,
                             x0, xn, z0, zn,
                             max_number,
                             1, lda*lda);

        /* Order nested dissection subparts */
        order_grid3D_classic(rangtab, peritab, cblknbr,
                             x0, xn, y0, y0+ny/2, z0, zn, max_number,
                             lda, current_rangtab,
                             treetab, current_treetab+1);

        order_grid3D_classic(rangtab, peritab, cblknbr,
                             x0, xn, y0+ny/2+1, yn, z0, zn, max_number,
                             lda, current_rangtab,
                             treetab, current_treetab+1);
    }

    /* If we cut in direction z */
    else{

        treetab[current_rangtab[0]] = current_treetab;
        rangtab[current_rangtab[0]] = max_number[0];
        current_rangtab[0]++;

        /* Order separator */
        pastix_int_t *peritab_separator = peritab + lda * lda * (z0 + nz/2);
        order_grid2D_classic(peritab_separator,
                             x0, xn, y0, yn,
                             max_number,
                             1, lda);

        /* Order nested dissection subparts */
        order_grid3D_classic(rangtab, peritab, cblknbr,
                             x0, xn, y0, yn, z0, z0+nz/2, max_number,
                             lda, current_rangtab,
                             treetab, current_treetab+1);

        order_grid3D_classic(rangtab, peritab, cblknbr,
                             x0, xn, y0, yn, z0+nz/2+1, zn, max_number,
                             lda, current_rangtab,
                             treetab, current_treetab+1);
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
 * @param[in,out] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *          On exit, the field oerdemesh is initialize with the result of the
 *          ordering realized by Scotch.
 *
 * @param[in] n
 *          The number of vertices of the graph
 *
 *******************************************************************************
 *
 * @return
 *          \retval PASTIX_SUCCESS on successful exit
 *
 *******************************************************************************/
int
orderComputeOptimal( pastix_data_t *pastix_data,
                     pastix_int_t n )
{
    pastix_int_t *iparm    = pastix_data->iparm;
    Order        *ordemesh = pastix_data->ordemesh;

    orderInit(ordemesh, n, n);

    pastix_int_t *rangtab = ordemesh->rangtab;
    pastix_int_t *permtab = ordemesh->permtab;
    pastix_int_t *peritab = ordemesh->peritab;
    pastix_int_t *treetab = ordemesh->treetab;

    pastix_int_t *saved_rangtab, *saved_treetab;

    pastix_int_t sep = iparm[IPARM_OPTIMAL_ORDERING];
    pastix_int_t i   = 2;

    pastix_int_t current_rangtab = 0;

    /* Graphs for using classical separators */
    while (i != sep && i < sep+1){
        i = 2*i+1;
    }
    if (i != sep){
        printf("\n\nThe given graph size is not correct for optimal manual ordering on 2D regular grid or 3D regular cube. Closer valid sizes are %ld %ld\n\n", i, 2*i+1);
    }

    ordemesh->cblknbr = 0;

    if (sep * sep == n){
        order_grid2D_wide(rangtab, permtab, &ordemesh->cblknbr,
                          0, sep, 0, sep, n-1, sep, &current_rangtab);
    }
    else{
        pastix_int_t current_number = n-1;
        order_grid3D_wide(rangtab, permtab, &ordemesh->cblknbr,
                          0, sep, 0, sep, 0, sep, &current_number, sep, &current_rangtab,
                          treetab, 1);
    }

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

