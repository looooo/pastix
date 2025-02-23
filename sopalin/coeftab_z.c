/**
 *
 * @file coeftab_z.c
 *
 * Precision dependent sequential routines to apply operation of the full matrix.
 *
 * @copyright 2015-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @author Esragul Korkmaz
 * @author Tony Delarue
 * @date 2024-07-05
 *
 * @precisions normal z -> s d c
 *
 **/
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define _GNU_SOURCE 1
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
#include "common.h"
#include "blend/solver.h"
#include "lapacke.h"
#include "sopalin/coeftab_z.h"
#include "pastix_zcores.h"

/**
 *******************************************************************************
 *
 * @brief Dump the solver matrix coefficients into a file in human readable
 * format.
 *
 * All non-zeroes coefficients are dumped in the format:
 *    i j val
 * with one value per row.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data instance to access the unique directory id in which
 *          output the files.
 *
 * @param[in] solvmtx
 *          The solver matrix to print.
 *
 * @param[in] prefix
 *          The filename where to store the output matrix.
 *
 *******************************************************************************/
void
coeftab_zdump( pastix_data_t      *pastix_data,
               const SolverMatrix *solvmtx,
               const char         *prefix )
{
    SolverCblk *cblk = solvmtx->cblktab;
    pastix_int_t itercblk;
    char  filename[256];
    FILE *stream = NULL;

    pastix_gendirectories( pastix_data );

    /*
     * TODO: there is a problem right here for now, because there are no
     * distinctions between L and U coeffcients in the final file
     */
    for (itercblk=0; itercblk<solvmtx->cblknbr; itercblk++, cblk++)
    {
        if ( cblk->cblktype & (CBLK_FANIN|CBLK_RECV) ) {
            continue;
        }
        if ( solvmtx->clustnum != cblk->ownerid ) {
            continue;
        }

        sprintf( filename, "%s_%ld.txt", prefix, (long)cblk->gcblknum );
        stream = pastix_fopenw( pastix_data->dir_global, filename, "w" );
        if ( stream == NULL ){
            continue;
        }

        cpucblk_zdump( PastixLCoef, cblk, stream );
        if ( NULL != cblk->ucoeftab ) {
            cpucblk_zdump( PastixUCoef, cblk, stream );
        }

        fclose( stream );
    }
}

/**
 *******************************************************************************
 *
 * @brief Dump a single column block into a FILE in a human readale format.
 *
 * All non-zeroes coefficients are dumped in the format:
 *    i j val
 * with one value per row.
 *
 * The filename is as follows : {L, U}cblk{Index of the cblk}_init.txt
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the matrix must be initialized.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[in] cblk
 *          The column block to dump into the file.
 *
 * @param[in] itercblk
 *          The index of the cblk to dump
 *
 * @param[inout] directory
 *          The pointer to the temporary directory where to store the output
 *          files.
 *
 *******************************************************************************/
void
cpucblk_zdumpfile( pastix_coefside_t side,
                   SolverCblk       *cblk,
                   pastix_int_t      itercblk,
                   const char       *directory )
{
    FILE *f = NULL;
    char *filename;
    int rc;

    /* Lower part */
    if ( side != PastixUCoef )
    {
        rc = asprintf( &filename, "Lcblk%05ld_init.txt", (long int) itercblk );
        f  = pastix_fopenw( directory, filename, "w" );
        if ( f != NULL ) {
            cpucblk_zdump( PastixLCoef, cblk, f );
            fclose( f );
        }
        free( filename );
    }

    /* Upper part */
    if ( side != PastixLCoef )
    {
        rc = asprintf( &filename, "Ucblk%05ld_init.txt", (long int) itercblk );
        f  = pastix_fopenw( directory, filename, "w" );
        if ( f != NULL ) {
            cpucblk_zdump( PastixUCoef, cblk, f );
            fclose( f );
        }
        free( filename );
    }
    (void)rc;
}

/**
 *******************************************************************************
 *
 * @brief Compare two solver matrices in full-rank format with the same data
 * distribution.
 *
 * The second solver matrix is overwritten by the difference of the two
 * matrices.  The frobenius norm of the difference of each column block is
 * computed and the functions returns 0 if the result for all the column blocks
 * of:
 *      || B_k - A_k || / ( || A_k || * eps )
 *
 * is below 10. Otherwise, an error message is printed and 1 is returned.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the cblk must be tested.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[in] solvA
 *          The solver matrix A.
 *
 * @param[inout] solvB
 *          The solver matrix B.
 *          On exit, B coefficient arrays are overwritten by the result of
 *          (B-A).
 *
 *******************************************************************************
 *
 * @return 0 if the test is passed, >= 0 otherwise.
 *
 *******************************************************************************/
int
coeftab_zdiff( pastix_coefside_t   side,
               const SolverMatrix *solvA,
               SolverMatrix       *solvB )
{
    SolverCblk *cblkA = solvA->cblktab;
    SolverCblk *cblkB = solvB->cblktab;
    pastix_int_t cblknum;
    int rc       = 0;
    int saved_rc = 0;

    for(cblknum=0; cblknum<solvA->cblknbr; cblknum++, cblkA++, cblkB++) {
        rc += cpucblk_zdiff( side, cblkA, cblkB );
        if ( rc != saved_rc ){
            fprintf(stderr, "CBLK %ld was not correctly compressed\n", (long)cblknum);
            saved_rc = rc;
        }
    }

    return rc;
}

/**
 *******************************************************************************
 *
 * @brief Compress all the cblks marked as valid for low-rank format.
 *
 * All the cblk in the top levels of the elimination tree marked as candidates
 * for compression are compressed if there is a gain to compress them. The
 * compression to low-rank format is parameterized by the input information
 * stored in the lowrank structure. On exit, all the cblks marked for
 * compression are stored through the low-rank structure, even if they are kept
 * in their full-rank form.
 *
 * @remark This routine is sequential
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The solver matrix of the problem to compress.
 *
 *******************************************************************************
 *
 * @return The memory gain resulting from the compression to low-rank format in
 * Bytes.
 *
 *******************************************************************************/
pastix_int_t
coeftab_zcompress( SolverMatrix *solvmtx )
{
    SolverCblk       *cblk = solvmtx->cblktab;
    pastix_coefside_t side = (solvmtx->factotype == PastixFactLU) ? PastixLUCoef : PastixLCoef;
    int               ilu_lvl;
    pastix_int_t      cblknum, gain = 0;

    ilu_lvl = solvmtx->lowrank.compress_preselect ? -1 : solvmtx->lowrank.ilu_lvl;

    for(cblknum=0; cblknum<solvmtx->cblknbr; cblknum++, cblk++) {
        if ( cblk->cblktype & CBLK_COMPRESSED ) {
            gain += cpucblk_zcompress( solvmtx, side, ilu_lvl, cblk );
        }
    }
    return gain;
}

/**
 *******************************************************************************
 *
 * @brief Uncompress all column block in low-rank format into full-rank format
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The solver matrix of the problem.
 *
 *******************************************************************************/
void
coeftab_zuncompress( SolverMatrix *solvmtx )
{
    SolverCblk  *cblk   = solvmtx->cblktab;
    pastix_int_t cblknum;
    pastix_coefside_t side = (solvmtx->factotype == PastixFactLU) ? PastixLUCoef : PastixLCoef;

    for(cblknum=0; cblknum<solvmtx->cblknbr; cblknum++, cblk++) {
        if (cblk->cblktype & CBLK_COMPRESSED) {
            cpucblk_zuncompress( side, cblk );
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the memory gain of the low-rank form over the full-rank form
 * for the entire matrix.
 *
 * This function returns the memory gain in bytes for the full matrix when
 * column blocks are stored in low-rank format compared to a full rank storage.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix of the problem.
 *
 * @param[in] iparm
 *          The integer parameter array
 *
 * @param[inout] dparm
 *          The double parameter array which is going to be updated.
 *
 *******************************************************************************/
void
coeftab_zmemory_lr( const SolverMatrix *solvmtx,
                    const pastix_int_t *iparm,
                    pastix_fixdbl_t    *dparm )
{
    pastix_coefside_t side = (solvmtx->factotype == PastixFactLU) ? PastixLUCoef : PastixLCoef;
    SolverCblk       *cblk = solvmtx->cblktab;
    const SolverBlok *blok;
    pastix_int_t i, cblknum;
    pastix_int_t    gain[MEMORY_STATS_SIZE] = { 0 };
    pastix_int_t    orig[MEMORY_STATS_SIZE] = { 0 };
    pastix_fixdbl_t memlr[MEMORY_STATS_SIZE] = { 0. };
    pastix_fixdbl_t memfr[MEMORY_STATS_SIZE] = { 0. };
    pastix_fixdbl_t totlr, totfr;

#if defined(PASTIX_SUPERNODE_STATS)
    pastix_int_t      last[3] = { 0 };
    pastix_fixdbl_t   memlast[4];
    const SolverBlok *solvblok = solvmtx->bloktab;

    for(i=0; i<solvmtx->bloknbr; i++, solvblok++ ) {
        const SolverCblk *lcblk = solvmtx->cblktab + solvblok->lcblknm;
        pastix_int_t      ncols = cblk_colnbr( lcblk );
        pastix_int_t      nrows = blok_rownbr( solvblok );
        pastix_int_t      size  = ncols * nrows;

        /* Skip remote data */
        if ( cblk->ownerid != solvmtx->clustnum ) {
            continue;
        }

        /* Let's skip recv and fanin for now */
        if ( lcblk->cblktype & (CBLK_RECV|CBLK_FANIN) ) {
            continue;
        }

        if ( !(lcblk->cblktype & CBLK_COMPRESSED) ) {
            if ( side != PastixLCoef ) {
                last[solvblok->inlast] += 2 * size;
            }
            else{
                last[solvblok->inlast] += size;
            }
        }
        else{
            if ( side != PastixUCoef ) {
                if ( solvblok->LRblock[0].rk >= 0 ) {
                    assert( solvblok->LRblock[0].rk <= core_get_rklimit( nrows, ncols ) );
                    assert( ((nrows+ncols) * solvblok->LRblock[0].rkmax) <= size );
                    last[solvblok->inlast] += ((nrows+ncols) * solvblok->LRblock[0].rkmax);
                }
                else {
                    last[solvblok->inlast] += size;
                }
            }

            if ( side != PastixLCoef ) {
                if ( solvblok->LRblock[1].rk >= 0 ) {
                    assert( solvblok->LRblock[1].rk <= core_get_rklimit( nrows, ncols ) );
                    assert( ((nrows+ncols) * solvblok->LRblock[1].rkmax) <= size );
                    last[solvblok->inlast] += ((nrows+ncols) * solvblok->LRblock[1].rkmax);
                }
                else {
                    last[solvblok->inlast] += size;
                }
            }
        }
    }
    for (i=0; i<3; i++) {
        memlast[i] = last[i] * pastix_size_of( PastixComplex64 );
    }
    memlast[3] = memlast[0] + memlast[1] + memlast[2];

    pastix_print( solvmtx->clustnum, 0,
                  "    Compression on LAST\n"
                  "      ------------------------------------------------\n"
                  "        A11                     %8.3g %co\n"
                  "        A12                     %8.3g %co\n"
                  "        A22                     %8.3g %co\n"
                  "        SUM                     %8.3g %co\n",
                  pastix_print_value(memlast[0]), pastix_print_unit(memlast[0]),
                  pastix_print_value(memlast[1]), pastix_print_unit(memlast[1]),
                  pastix_print_value(memlast[2]), pastix_print_unit(memlast[2]),
                  pastix_print_value(memlast[3]), pastix_print_unit(memlast[3]));
#endif

    for(cblknum=0; cblknum<solvmtx->cblknbr; cblknum++, cblk++) {
        pastix_int_t colnbr = cblk_colnbr( cblk );

        /* Skip remote data */
        if ( cblk->ownerid != solvmtx->clustnum ) {
            continue;
        }

        /* Let's skip recv and fanin for now */
        if ( cblk->cblktype & (CBLK_RECV|CBLK_FANIN) ) {
            continue;
        }

        if ( !(cblk->cblktype & CBLK_COMPRESSED) )
        {
            pastix_int_t in_height = 0;
            pastix_int_t off_height = cblk->stride;

            /* Compute the size of the original supernode diagonal block */
            blok = cblk->fblokptr;
            while( (blok < cblk[1].fblokptr) &&
                   ((solvmtx->cblktab + blok->fcblknm)->sndeidx == cblk->sndeidx) )
            {
                in_height += blok_rownbr( blok );
                blok++;
            }

            /* Size of the cblk outside the diagonal block */
            off_height -= in_height;

            orig[FR_InDiag]  += colnbr * in_height;
            orig[FR_OffDiag] += colnbr * off_height;
        }
        else {
            /* The gain for the diagonal block is always 0 */
            orig[LR_DInD] += colnbr * colnbr;
            cpucblk_zmemory( side, solvmtx, cblk, orig, gain );
        }
    }

    if ( side == PastixLUCoef ) {
        orig[FR_InDiag]  *= 2;
        orig[FR_OffDiag] *= 2;
        orig[LR_InDiag]  *= 2;
        orig[LR_OffDiag] *= 2;
        orig[LR_InSele]  *= 2;
    }

    totlr = 0.;
    totfr = 0.;

    for (i=0; i<MEMORY_STATS_SIZE; i++) {
        memlr[i] = (orig[i] - gain[i]) * pastix_size_of( PastixComplex64 );
        memfr[i] =  orig[i]            * pastix_size_of( PastixComplex64 );
        totlr += memlr[i];
        totfr += memfr[i];
    }

    if( iparm[IPARM_VERBOSE] > PastixVerboseNot ) {
        pastix_print( solvmtx->clustnum, 0,
                      "    Compression:\n"
                      "      ------------------------------------------------\n"
                      "      Full-rank supernodes\n"
                      "        Inside                                %8.3g %co\n"
                      "        Outside                               %8.3g %co\n"
                      "      Low-rank supernodes\n"
                      "        Diag in diag                          %8.3g %co\n"
                      "        Inside not selected     %8.3g %co / %8.3g %co\n"
                      "        Inside selected         %8.3g %co / %8.3g %co\n"
                      "        Outside                 %8.3g %co / %8.3g %co\n"
                      "      ------------------------------------------------\n"
                      "      Total                     %8.3g %co / %8.3g %co\n",
                      pastix_print_value(memfr[FR_InDiag] ), pastix_print_unit(memfr[FR_InDiag] ),
                      pastix_print_value(memfr[FR_OffDiag]), pastix_print_unit(memfr[FR_OffDiag]),

                      pastix_print_value(memfr[LR_DInD]),    pastix_print_unit(memfr[LR_DInD]),

                      pastix_print_value(memlr[LR_InDiag] ), pastix_print_unit(memlr[LR_InDiag] ),
                      pastix_print_value(memfr[LR_InDiag] ), pastix_print_unit(memfr[LR_InDiag] ),

                      pastix_print_value(memlr[LR_InSele] ), pastix_print_unit(memlr[LR_InSele]),
                      pastix_print_value(memfr[LR_InSele] ), pastix_print_unit(memfr[LR_InSele]),

                      pastix_print_value(memlr[LR_OffDiag]), pastix_print_unit(memlr[LR_OffDiag]),
                      pastix_print_value(memfr[LR_OffDiag]), pastix_print_unit(memfr[LR_OffDiag]),

                      pastix_print_value(totlr),             pastix_print_unit(totlr),
                      pastix_print_value(totfr),             pastix_print_unit(totfr) );
    }

    dparm[DPARM_MEM_FR] = totfr;
    dparm[DPARM_MEM_LR] = totlr;

    return;
}

/**
 *******************************************************************************
 *
 * @brief Compute the memory usage of the full-rank form for the entire matrix.
 *
 * This function returns the memory usage in bytes for the full matrix when
 * column blocks are stored in full-rank format.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix of the problem.
 *
 * @param[inout] dparm
 *          The double parameter array which is going to be updated.
 *
 * @param[in] iparm
 *          The integer parameter array
 *
 *******************************************************************************/
void
coeftab_zmemory_fr( const SolverMatrix *solvmtx,
                    const pastix_int_t *iparm,
                    pastix_fixdbl_t    *dparm )
{
    pastix_coefside_t side = (solvmtx->factotype == PastixFactLU) ? PastixLUCoef : PastixLCoef;
    const SolverCblk *cblk = solvmtx->cblktab;
    pastix_int_t      cblknum;
    pastix_fixdbl_t   totmem = 0.;

    for(cblknum=0; cblknum<solvmtx->cblknbr; cblknum++, cblk++) {
        pastix_int_t colnbr = cblk_colnbr( cblk );
        pastix_int_t rownbr = cblk->stride;

        /* Skip remote data */
        if ( cblk->ownerid != solvmtx->clustnum ) {
            continue;
        }

        /* Let's skip recv and fanin for now */
        if ( cblk->cblktype & (CBLK_RECV|CBLK_FANIN) ) {
            continue;
        }

        totmem += (double)colnbr * (double)rownbr;
    }

    if ( side == PastixLUCoef ) {
        totmem *= 2.;
    }

    totmem *= (double)pastix_size_of( PastixComplex64 );

    dparm[DPARM_MEM_FR] = totmem;

    if( iparm[IPARM_VERBOSE] > PastixVerboseNot ) {
        pastix_print( solvmtx->clustnum, 0,
                      "    Memory usage of coeftab                   %8.3g %co\n",
                      pastix_print_value(dparm[DPARM_MEM_FR]), pastix_print_unit(dparm[DPARM_MEM_FR]) );
    }

    return;
}

/**
 *******************************************************************************
 *
 * @brief Compute the memory usage for the entire matrix.
 *
 * This functions computes the memory usage and gain if the matrix is compressed
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix of the problem.
 *
 * @param[in] iparm
 *          The integer parameter array. Uses IPARM_COMPRESS_WHEN, IPARM_VERBOSE
 *
 * @param[inout] dparm
 *          The double parameter array which is going to be updated.
 *          Update DPARM_MEM_FR and DPARM_MEM_LR.
 *
 *******************************************************************************/
void
coeftab_zmemory( const SolverMatrix *solvmtx,
                 const pastix_int_t *iparm,
                 pastix_fixdbl_t    *dparm )
{
    if ( iparm[IPARM_COMPRESS_WHEN] != PastixCompressNever ) {
        coeftab_zmemory_lr( solvmtx, iparm, dparm );
    }
    else {
        coeftab_zmemory_fr( solvmtx, iparm, dparm );
    }
}

/**
 *******************************************************************************
 *
 * @brief Extract the Schur complement
 *
 * This routine is sequential and returns the full Schur complement
 * uncommpressed in Lapack format.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix structure describing the problem.
 *
 * @param[inout] S
 *          The pointer to the allocated matrix array that will store the Schur
 *          complement.
 *
 * @param[in] lds
 *          The leading dimension of the S array.
 *
 *******************************************************************************/
void
coeftab_zgetschur( const SolverMatrix *solvmtx,
                   pastix_complex64_t *S, pastix_int_t lds )
{
    SolverCblk *cblk = solvmtx->cblktab + solvmtx->cblkschur;
    pastix_complex64_t *localS;
    pastix_int_t itercblk, fcolnum, nbcol;
    pastix_int_t ret;
    int upper_part = (solvmtx->factotype == PastixFactLU);
    fcolnum = cblk->fcolnum;

    nbcol = solvmtx->nodenbr - fcolnum;
    assert( nbcol <= lds );

    /* Initialize the array to 0 */
    ret = LAPACKE_zlaset_work( LAPACK_COL_MAJOR, 'A', nbcol, nbcol, 0., 0., S, lds );
    assert( ret == 0 );

    for (itercblk=solvmtx->cblkschur; itercblk<solvmtx->cblknbr; itercblk++, cblk++)
    {
        assert( cblk->cblktype & CBLK_IN_SCHUR );
        assert( lds >= cblk->stride );

        localS = S + (cblk->fcolnum - fcolnum) * lds + (cblk->fcolnum - fcolnum);

        cpucblk_zgetschur( cblk, upper_part, localS, lds );
    }
    (void)ret;
}

/**
 *******************************************************************************
 *
 * @brief Extract the diagonal
 *
 * This routine is sequential and returns the full diagonal in the vector D,
 * such that:
 *     D[incD*i]= A(i, i)
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix structure describing the problem.
 *
 * @param[inout] D
 *          The pointer to the allocated vector array that will store the diagonal.
 *          D must be of size solvmtx->nodenbr * incD.
 *
 * @param[in] incD
 *          The increment bewteen two elements of D. incD > 0.
 *
 *******************************************************************************/
void
coeftab_zgetdiag( const SolverMatrix *solvmtx,
                  pastix_complex64_t *D, pastix_int_t incD )
{
    SolverCblk *cblk = solvmtx->cblktab;
    pastix_complex64_t *A;
    pastix_int_t lda, itercblk, nbcol, i;

    for (itercblk=0; itercblk<solvmtx->cblknbr; itercblk++, cblk++)
    {
        nbcol = cblk_colnbr( cblk );
        if ( cblk->cblktype & CBLK_COMPRESSED ) {
            assert( cblk->fblokptr->LRblock[0]->rk == -1 );
            A   = cblk->fblokptr->LRblock[0]->u;
            lda = cblk_colnbr( cblk ) + 1;
        }
        else {
            A = cblk->lcoeftab;

            if ( cblk->cblktype & CBLK_LAYOUT_2D ) {
                lda = cblk_colnbr( cblk ) + 1;
            }
            else {
                lda = cblk->stride + 1;
            }
        }

        for (i=0; i<nbcol; i++, D += incD, A += lda ) {
            *D = *A;
        }
    }
}
