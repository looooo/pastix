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


void order_grid2D_wide(pastix_int_t *rangtab,
                       pastix_int_t *peritab,
                       pastix_int_t *cblknbr,
                       pastix_int_t x0,
                       pastix_int_t xn,
                       pastix_int_t y0,
                       pastix_int_t yn,
                       pastix_int_t max_number,
                       pastix_int_t lda,
                       pastix_int_t *current_rangtab){

    pastix_int_t nx = xn-x0;
    pastix_int_t ny = yn-y0;
    /* printf("Treating a subgraph of size %ld %ld\n", nx, ny); */

    /* The subgraph is small enough */
    if (nx <= 4 && ny <= 4){
        cblknbr[0] ++;
        /* printf("Treating cblk %ld\n", cblknbr[0]); */
        pastix_int_t i, j;
        pastix_int_t current = 0;
        for (i=x0; i<xn; i++){
            for (j=y0; j<yn; j++){
                pastix_int_t index = i + lda * j;
                peritab[index] = max_number - current;
                /* printf("Fitting 1a %ld with %ld\n", index, max_number - current); */
                current++;
            }
        }

        rangtab[current_rangtab[0]] = max_number;
        current_rangtab[0]++;
        return;
    }

    cblknbr[0] ++;
    /* printf("Treating cblk %ld\n", cblknbr[0]); */

    /* In which direction do we cut? 0 for x, 1 for y */
    pastix_int_t dir = 0;
    if (ny > nx)
        dir = 1;

    /* If we cut in direction x */
    if (dir == 0){

        rangtab[current_rangtab[0]] = max_number;
        current_rangtab[0]++;

        pastix_int_t i;
        for (i=0; i<ny; i++){
            pastix_int_t index = x0 + nx/2 + lda*(y0+i) - 1;
            peritab[index] = max_number - i;
            /* printf("Fitting 1a %ld with %ld\n", index, max_number - i); */
        }
        for (i=0; i<ny; i++){
            pastix_int_t index = x0 + nx/2 + lda*(y0+i);
            peritab[index] = max_number - ny - i;
            /* printf("Fitting 1b %ld with %ld\n", index, max_number - ny - i); */
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

        pastix_int_t i;
        for (i=0; i<nx; i++){
            pastix_int_t index = lda*(y0+ny/2-1) + x0 + i;
            peritab[index] = max_number - i;
            /* printf("Fitting 2a %ld with %ld\n", index, max_number - i); */
        }
        for (i=0; i<nx; i++){
            pastix_int_t index = lda*(y0+ny/2) + x0 + i;
            peritab[index] = max_number - nx - i;
            /* printf("Fitting 2b %ld with %ld\n", index, max_number - nx - i); */
        }

        /* printf("Parameters %ld %ld %ld %ld\n", x0, xn, y0, yn/2-1); */
        order_grid2D_wide(rangtab, peritab, cblknbr,
                          x0, xn, y0, yn-ny/2-1, max_number - 2*nx,
                          lda, current_rangtab);

        order_grid2D_wide(rangtab, peritab, cblknbr,
                          x0, xn, y0+ny/2+1, yn, max_number - 2*nx - nx*(ny/2-1),
                          lda, current_rangtab);
    }
}

void order_grid3D_wide(pastix_int_t *rangtab,
                       pastix_int_t *peritab,
                       pastix_int_t *cblknbr,
                       pastix_int_t x0,
                       pastix_int_t xn,
                       pastix_int_t y0,
                       pastix_int_t yn,
                       pastix_int_t z0,
                       pastix_int_t zn,
                       pastix_int_t max_number,
                       pastix_int_t lda,
                       pastix_int_t *current_rangtab,
                       pastix_int_t *treetab,
                       pastix_int_t current_treetab){

    pastix_int_t nx = xn-x0;
    pastix_int_t ny = yn-y0;
    pastix_int_t nz = zn-z0;
    /* printf("Treating a subgraph of size %ld %ld %ld\n", nx, ny, nz); */

    /* The subgraph is small enough */
    if (nx <= 4 && ny <= 4 && nz <= 4){
        cblknbr[0] ++;
        /* printf("Treating cblk %ld\n", cblknbr[0]); */
        pastix_int_t i, j, k;
        pastix_int_t current = 0;
        for (i=x0; i<xn; i++){
            for (j=y0; j<yn; j++){
                for (k=z0; k<zn; k++){
                    pastix_int_t index = i + lda * j + lda*lda * k;
                    peritab[index] = max_number - current;
                    current++;
                }
            }
        }

        treetab[current_rangtab[0]] = current_treetab;
        rangtab[current_rangtab[0]] = max_number;
        current_rangtab[0]++;
        return;
    }

    cblknbr[0] ++;
    /* printf("Treating cblk %ld\n", cblknbr[0]); */

    /* In which direction do we cut? 0 for x, 1 for y */
    pastix_int_t dir = 0;
    if (ny > nx)
        dir = 1;
    if (nz > nx && nz > ny)
        dir = 2;

    /* If we cut in direction x */
    if (dir == 0){

        treetab[current_rangtab[0]] = current_treetab;
        rangtab[current_rangtab[0]] = max_number;
        current_rangtab[0]++;

        pastix_int_t i, j;
        pastix_int_t current = 0;
        for (i=0; i<ny; i++){
            for (j=0; j<nz; j++){
                pastix_int_t index = x0 + nx/2 + lda*(y0+i) - 1 + lda*lda*(z0+j);
                peritab[index] = max_number - current;
                current++;
            }
        }
        for (i=0; i<ny; i++){
            for (j=0; j<nz; j++){
                pastix_int_t index = x0 + nx/2 + lda*(y0+i) + lda*lda*(z0+j);
                peritab[index] = max_number - current;
                current++;
            }
        }

        order_grid3D_wide(rangtab, peritab, cblknbr,
                          x0, xn-nx/2-1, y0, yn, z0, zn, max_number - 2*ny*nz,
                          lda, current_rangtab,
                          treetab, current_treetab+1);

        order_grid3D_wide(rangtab, peritab, cblknbr,
                          x0+nx/2+1, xn, y0, yn, z0, zn, max_number - nx*ny*nz /2 - ny*nz,
                          lda, current_rangtab,
                          treetab, current_treetab+1);

    }

    /* If we cut in direction y */
    else if (dir == 1){

        treetab[current_rangtab[0]] = current_treetab;
        rangtab[current_rangtab[0]] = max_number;
        current_rangtab[0]++;

        pastix_int_t i, j;
        pastix_int_t current = 0;
        for (i=0; i<nx; i++){
            for (j=0; j<nz; j++){
                pastix_int_t index = lda*(y0+ny/2-1) + x0 + i  + lda*lda*(z0+j);
                peritab[index] = max_number - current;
                current++;
            }
        }
        for (i=0; i<nx; i++){
            for (j=0; j<nz; j++){
                pastix_int_t index = lda*(y0+ny/2) + x0 + i + lda*lda*(z0+j);
                peritab[index] = max_number - current;
                current++;
            }
        }

        order_grid3D_wide(rangtab, peritab, cblknbr,
                          x0, xn, y0, yn-ny/2-1, z0, zn, max_number - 2*nx*nz,
                          lda, current_rangtab,
                          treetab, current_treetab+1);

        order_grid3D_wide(rangtab, peritab, cblknbr,
                          x0, xn, y0+ny/2+1, yn, z0, zn, max_number - nx*ny*nz /2 - nx*nz,
                          lda, current_rangtab,
                          treetab, current_treetab+1);
    }

    /* If we cut in direction y */
    else{

        treetab[current_rangtab[0]] = current_treetab;
        rangtab[current_rangtab[0]] = max_number;
        current_rangtab[0]++;

        pastix_int_t i, j;
        pastix_int_t current = 0;
        for (i=0; i<nx; i++){
            for (j=0; j<ny; j++){
                pastix_int_t index = lda*lda*(z0+nz/2-1) + x0 + i + lda * (y0+j);
                peritab[index] = max_number - current;
                current++;
            }
        }
        for (i=0; i<nx; i++){
            for (j=0; j<ny; j++){
                pastix_int_t index = lda*lda*(z0+nz/2) + x0 + i + lda * (y0+j);
                peritab[index] = max_number - current;
                current++;
            }
        }

        order_grid3D_wide(rangtab, peritab, cblknbr,
                          x0, xn, y0, yn, z0, zn-nz/2-1, max_number - 2*nx*ny,
                          lda, current_rangtab,
                          treetab, current_treetab+1);

        order_grid3D_wide(rangtab, peritab, cblknbr,
                          x0, xn, y0, yn, z0+nz/2+1, zn, max_number - nx*ny*nz /2 - nx*ny,
                          lda, current_rangtab,
                          treetab, current_treetab+1);
    }
}

void order_grid2D_classic(pastix_int_t *peritab,
                          pastix_int_t x0,
                          pastix_int_t xn,
                          pastix_int_t y0,
                          pastix_int_t yn,
                          pastix_int_t *max_number,
                          pastix_int_t ldax,
                          pastix_int_t lday){

    pastix_int_t nx = xn-x0;
    pastix_int_t ny = yn-y0;

    /* The subgraph is small enough */
    /* if (nx <= 4 && ny <= 4 && nz <= 4){ */
    if (nx*ny < 10){
        pastix_int_t i, j;
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
    pastix_int_t dir = 0;
    if (ny > nx)
        dir = 1;

    /* If we cut in direction x */
    if (dir == 0){

        pastix_int_t i;
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

        pastix_int_t i;
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

void order_grid3D_classic(pastix_int_t *rangtab,
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
                          pastix_int_t current_treetab){

    pastix_int_t nx = xn-x0;
    pastix_int_t ny = yn-y0;
    pastix_int_t nz = zn-z0;

    /* printf("Treating a subgraph of size %ld %ld %ld\n", nx, ny, nz); */

    /* The subgraph is small enough */
    /* if (nx <= 4 && ny <= 4 && nz <= 4){ */
    if (nx*ny*nz < 10){
        cblknbr[0] ++;
        /* printf("Treating small clique %ld\n", nx*ny*nz); */
        pastix_int_t i, j, k;
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
    pastix_int_t dir = 0;
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
