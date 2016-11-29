/**
 *
 * @file order_compute_scotch.c
 *
 *  PaStiX order routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * Contains functions to perform ordering with Scotch library.
 *
 * @version 5.1.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 */
#include "common.h"
#include <string.h>
#include "graph.h"
#include "order.h"
#if defined(PASTIX_ORDERING_PTSCOTCH)
#include <ptscotch.h>
#elif defined(PASTIX_ORDERING_SCOTCH)
#include <scotch.h>
#endif /* defined(PASTIX_ORDERING_PTSCOTCH) */

#define COLOR 3

static pastix_int_t dichotomic_search( pastix_int_t node, const Order *order )
{
    pastix_int_t cblknum;
    pastix_int_t first, last, middle;

    first   = 0;
    last    = order->cblknbr;
    middle  = (last+first) / 2;
    cblknum = -1;

    while( last - first > 0 ) {
        if ( node >= order->rangtab[middle] ) {
            if ( node < order->rangtab[middle+1] ) {
                cblknum = middle;
                break;
            }
            first = middle;
        }
        else {
            last = middle;
        }
        middle = (last+first) / 2;
    }
    assert( (order->rangtab[cblknum]   <= node) &&
            (order->rangtab[cblknum+1] >  node) );
    return cblknum;
}

static inline pastix_int_t
compute_cblklevel( const pastix_int_t *treetab,
                   pastix_int_t *levels,
                   pastix_int_t  cblknum )
{
    /* cblknum level has already been computed */
    if ( levels[cblknum] != 0 ) {
        return levels[cblknum];
    }
    else {
        pastix_int_t father = treetab[cblknum];

        if ( father == -1 ) {
            return 1;
        }
        else {
            return compute_cblklevel( treetab, levels, father ) + 1;
        }
    }
}

void
orderComputeClif( const pastix_graph_t *graph,
                  SCOTCH_Graph         *sgraph,
                  Order                *order,
                  SCOTCH_Ordering      *sorder)
{
    SCOTCH_Graph    sn_sgraph;
    SCOTCH_Strat    sn_strat;
    FILE *file;
    pastix_int_t  *sn_parttab;
    pastix_int_t   i, vertnbr = order->rangtab[order->cblknbr - 1];
    pastix_int_t   sn_vertnbr = graph->n - vertnbr; /* This works for the last supernode */
    SCOTCH_graphBase( sgraph, 0 );
    orderBase( order, 0 );

    //clifton: node reordering somewhere around here

    /* Output graph before any changes */
    //DEBUG: Before partition, output entire graph, map
    file = fopen("before.grf","w");
    SCOTCH_graphSave(sgraph, file);
    fclose(file);

    /* { */
    /*     pastix_int_t i; */
    /*     fprintf(stderr, "cblknbr = %ld:\n", order->cblknbr); */
    /*     for (i=0; i<=order->cblknbr; i++) { */
    /*         fprintf(stderr, "%ld ", order->rangtab[i]); */
    /*         if (i % 8 == 7 ) */
    /*             fprintf(stderr, "\n"); */
    /*     } */
    /*     fprintf(stderr, "\n"); */
    /* } */

    file = fopen("before.map","w");
    SCOTCH_graphOrderSaveMap(sgraph, sorder, file);
    fclose(file);

    file = fopen("before.ord","w");
    SCOTCH_graphOrderSave(sgraph, sorder, file);
    fclose(file);

    fprintf(stderr,"Before Files complete\n");

    //try graphIsolate here
    //get all the nodes except the ones in the block we want to keep
    /**
     * Extract the subgraph with unknowns of the last supernode
     * For that, we need to generate the list of all vertices in original
     * numbering that are not in this supernode.
     */
     {
        pastix_int_t *sn_colptr, *sn_rows;

        pastix_int_t sn_id  = order->cblknbr - 1; /* set to whichever supernode should be extracted */
        pastix_int_t fnode = order->rangtab[sn_id];
        pastix_int_t lnode = order->rangtab[sn_id+1]-1;
        graphIsolateSupernode( graph->n,
			       graph->colptr,
			       graph->rows,
                               order->permtab,
                               order->peritab,
                               fnode,
                               lnode,
                               sn_vertnbr,
			       order->peritab + order->rangtab[order->cblknbr - 1],
			       &sn_colptr,
			       &sn_rows );
        /*{
        pastix_int_t *non_blk_vertices;
        pastix_int_t *sn_colptr, *sn_rows, *sn_perm, *sn_invp;

        MALLOC_INTERN(non_blk_vertices, vertnbr, pastix_int_t);

        for(i=0; i<vertnbr; i++) {
            non_blk_vertices[i] = order->peritab[i];
        }

        graphIsolate( graph->n,
                      graph->colptr,
                      graph->rows,
                      vertnbr,
                      non_blk_vertices,
                      &sn_colptr,
                      &sn_rows,
                      &sn_perm,
                      &sn_invp );
         */
        //csc_symgraph_int(blksize,ncol_ptr,nrows,NULL,&nblksize,&nncol,&nnrow,NULL,API_YES);
        //csc_noDiag(nncol[0],nblksize,nncol,nnrow,NULL);

        if(!SCOTCH_graphInit(&sn_sgraph))
        {
            SCOTCH_graphBuild( &sn_sgraph,
                               order->baseval,
                               sn_vertnbr,
                               sn_colptr,
                               NULL,
                               NULL,
                               NULL,
                               sn_colptr[ sn_vertnbr ] - order->baseval,
                               sn_rows,
                               NULL);
        }
        else
            fprintf(stderr,"Failed to build graph\n");

        fprintf(stderr,"Check: %d\n", SCOTCH_graphCheck( &sn_sgraph ));

        //Initialize partitioning strategy
        if(SCOTCH_stratInit(&sn_strat) != 0 )
            fprintf(stderr,"Failed to initialize partitioning strategy\n");

        /* Create two partitions of the supernode */
        MALLOC_INTERN( sn_parttab, sn_vertnbr, pastix_int_t );

        /* if(SCOTCH_graphPart( &sn_sgraph, 2, &sn_strat, sn_parttab ) != 0) */
        /*     fprintf(stderr,"Partitioning Failed\n"); */

        file = fopen("part.grf","w");
        SCOTCH_graphSave(&sn_sgraph, file);
        fclose(file);

        /* file = fopen("part.map","w"); */
        /* SCOTCH_graphMapSave( &sn_sgraph, &mappdat, file ); */
        /* fclose(file); */

        if(0)
        {
            pastix_int_t   partnbr = 2;
            SCOTCH_Arch    archdat;
            SCOTCH_Mapping mappdat;

            SCOTCH_archInit  (&archdat);
            SCOTCH_archCmplt (&archdat, partnbr);
            SCOTCH_graphMapInit (&sn_sgraph, &mappdat, &archdat, sn_parttab);

            SCOTCH_graphMapCompute (&sn_sgraph, &mappdat, &sn_strat);

            file = fopen("part.grf","w");
            SCOTCH_graphSave(&sn_sgraph, file);
            fclose(file);


            SCOTCH_graphMapExit (&sn_sgraph, &mappdat);
            SCOTCH_archExit (&archdat);
        }
        fprintf(stderr, "After partitioning\n");

        /* Extract the data from the xyz file */
        {
            FILE *fileout;
            pastix_int_t   n, dim;
            int rc;
            (void) rc;

            file = fopen( "before.xyz", "r" );
            fileout = fopen( "part.xyz", "w" );

            rc = fscanf( file, "%ld %ld", &dim, &n );
            assert( n == graph->n );
            assert( rc == 2 );

            if ( dim == 2 ) {
                pastix_int_t fnode = order->rangtab[ order->cblknbr-1 ];

                /* Write header */
                fprintf(fileout, "%ld %ld\n", dim, sn_vertnbr );

                for(i=0; i<n; i++) {
                    pastix_int_t v, iv;
                    double x, y;

                    rc = fscanf(file, "%ld %lf %lf", &v, &x, &y );
                    assert( rc == 3 );

                    /* If permutation in the last supernode, we keep it */
                    iv = order->permtab[i];
                    if ( iv >= fnode ) {
                        fprintf(fileout, "%ld %lf %lf\n", iv-fnode, x, y);
                    }
                }
            }
            else if (dim == 3) {
                pastix_int_t fnode = order->rangtab[ order->cblknbr-1 ];

                /* Write header */
                fprintf(fileout, "%ld %ld\n", dim, sn_vertnbr );

                for(i=0; i<n; i++) {
                    pastix_int_t v, iv;
                    double x, y, z;

                    rc = fscanf(file, "%ld %lf %lf %lf", &v, &x, &y, &z );
                    assert( rc == 4 );

                    /* If permutation in the last supernode, we keep it */
                    iv = order->permtab[i];
                    if ( iv >= fnode ) {
                        fprintf(fileout, "%ld %lf %lf %lf\n", iv-fnode, x, y, z);
                    }
                }
            }
            else {
                fprintf(stderr, "Read xyz for dim=%ld not implemented\n", dim );
            }

            fclose(file);
            fclose(fileout);


            if (dim == 3){
                int *colors = malloc((lnode-fnode+1) * sizeof(int));
                for (i=0; i<lnode-fnode+1; i++){
                    colors[i] = COLOR;
                }
                //memset(colors, 5, (lnode-fnode+1) * sizeof(int));

                pastix_int_t *levels;
                levels = calloc(order->cblknbr, sizeof(pastix_int_t));

                for (i=0; i<order->cblknbr; i++){
                    levels[i] = compute_cblklevel( order->treetab, levels, i );
                }


                pastix_int_t *colptr = graph->colptr;
                pastix_int_t *rows   = graph->rows;
                pastix_int_t *invp   = order->peritab;
                pastix_int_t *perm   = order->permtab;
                pastix_int_t baseval = colptr[0];
                pastix_int_t j, k;
                pastix_int_t ip;
                for (ip=fnode; ip<=lnode; ip++){
                    i = invp[ip];
                    for (j = colptr[i]-baseval; j < colptr[i+1]-baseval; j++){
                        pastix_int_t jp = perm[ rows[j]-baseval ];
                        pastix_int_t cblknbr = dichotomic_search( jp, order );
                        if (levels[cblknbr] == 2 && colors[ip-fnode] == COLOR){
                            colors[ip-fnode] = 0;
                        }
                        if (levels[cblknbr] == 3 && colors[ip-fnode] != 2){
                            colors[ip-fnode] = 1;
                        }
                        if (levels[cblknbr] == 4 && colors[ip-fnode] == COLOR){
                            colors[ip-fnode] = 2;
                        }


                        for(k = colptr[rows[j]-baseval]-baseval; k < colptr[rows[j]+1-baseval]-baseval; k++){
                            pastix_int_t kp = perm[ rows[k]-baseval ];
                            pastix_int_t cblknbr = dichotomic_search( kp, order );
                            if (levels[cblknbr] == 2  && colors[ip-fnode] == COLOR){
                                colors[ip-fnode] = 0;
                            }
                            if (levels[cblknbr] == 3 && colors[ip-fnode] != 2){
                                colors[ip-fnode] = 1;
                            }
                            if (levels[cblknbr] == 4 && colors[ip-fnode] == COLOR){
                                colors[ip-fnode] = 2;
                            }
                        }
                    }
                }

                pastix_int_t v, c, iv;
                (void) iv;
                (void) c;
                (void) v;
                fileout    = fopen("part.map" ,"w");
                /* fileout = fopen("part.map2","w"); */
                /* rc = fscanf(file, "%ld", &v); */

                /* printf("NB NODES %ld %ld\n", v, lnode-fnode+1); */
                /* assert ( v == (lnode-fnode+1) ); */

                /* Write header */
                fprintf(fileout, "%ld\n", lnode-fnode+1);

                for (i=fnode; i<=lnode; i++){
                    /* rc = fscanf(file, "%ld %ld", &v, &c); */
                    /* assert( rc == 2 ); */
                    fprintf(fileout, "%ld %d\n", i-fnode, colors[i-fnode]);
                }

                /* fclose(file); */
                fclose(fileout);
            }
        }

        /* Update the invp/perm arrays */
    /* pastix_int_t  *new_rangtab; */
    /* pastix_int_t   new_cblknbr; */

    /*     if (1) */
    /*     { */
    /*         int idxones[ sn_vertnbr]; */
    /*         int localoff[sn_vertnbr]; */
    /*         int offzero = 0; */
    /*         int offone = 0; */
    /*         int crdi; */

    /*         for(crdi = 0; crdi < sn_vertnbr; crdi++) */
    /*         { */
    /*             if(sn_parttab[crdi]) */
    /*             { */
    /*                 idxones[offone] = crdi; */
    /*                 localoff[crdi] = offone++; */
    /*             } */
    /*             else */
    /*                 localoff[crdi] = offzero++; */
    /*         } */
    /*         for(crdi = 0; crdi < offone; crdi++) */
    /*         { */
    /*             localoff[idxones[crdi]] += offzero; */
    /*         } */

    /*         for(crdi = 0; crdi < sn_vertnbr; crdi++) */
    /*         { */
    /*             int cidx = sn_invp[crdi]; */
    /*             int npm = vertnbr+localoff[crdi]; */
    /*             order->permtab[cidx] = npm; */
    /*             order->peritab[npm] = cidx; */
    /*         } */

    /*         /\* Update rangtab with a two partition of the supernode *\/ */
    /*         new_cblknbr = order->cblknbr + 1; /\* Add Nb partition - 1 *\/ */
    /*         MALLOC_INTERN( new_rangtab, new_cblknbr+1, pastix_int_t ); */

    /*         memcpy( new_rangtab, order->rangtab, (order->cblknbr+1) * sizeof(pastix_int_t) ); */
    /*         new_rangtab[order->cblknbr]   = new_rangtab[order->cblknbr-1] + offzero; */
    /*         new_rangtab[order->cblknbr+1] = order->vertnbr; */
    /*     } */
    }

    /* { */
    /*     pastix_int_t i; */
    /*     fprintf(stderr, "cblknbr = %ld:\n", new_cblknbr); */
    /*     for (i=0; i<=new_cblknbr; i++) { */
    /*         fprintf(stderr, "%ld ", new_rangtab[i]); */
    /*         if (i % 8 == 7 ) */
    /*             fprintf(stderr, "\n"); */
    /*     } */
    /*     fprintf(stderr, "\n"); */
    /* } */
}
