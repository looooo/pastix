/**
 *
 * @file pzlarft_blgtrd.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Azzam Haidar
  @date 2011-05-15
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include <lapacke.h>

#define V(m)     &(V[(m)])
#define TAU(m)   &(TAU[(m)])
#define T(m)   &(T[(m)])
/***************************************************************************/
/**
 *  Parallel compute T2 from bulgechasing of Symetric matrix - static scheduling
 *  Lower case is supported
 **/
/***************************************************************************/
void plasma_pzlarft_blgtrd(plasma_context_t *plasma)
{
    int my_core_id = PLASMA_RANK;
    int cores_num  = plasma->world_size;
    /*===========================*/
    int N, NB,Vblksiz;
    PLASMA_Complex64_t *V;
    PLASMA_Complex64_t *T;
    PLASMA_Complex64_t *TAU;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    /*===========================
     *   local variables
     *===========================*/
    int LDT, LDV;
    int Vm, Vn, mt, nt;
    int myrow, mycol, blkj, blki;
    int firstrow;
    int blkid,vpos,taupos,tpos;
    int blkpercore,blkcnt, myid;


    plasma_unpack_args_8(N, NB, Vblksiz, V, T, TAU, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;

    /* Quick return */
    if (N == 0){
        return;
    }
    if (NB == 0){
        return;
    }
    if (NB == 1){
        return;
    }

    findVTsiz(N, NB, Vblksiz, &blkcnt, &LDV);
    blkpercore = blkcnt/cores_num;
    blkpercore = blkpercore==0 ? 1:blkpercore;
    LDT        = Vblksiz;    
    LDV        = NB+Vblksiz-1;

    /*========================================
     * compute the T's in parallel.
     * The Ts are independent so each core pick
     * a T and compute it. The loop is based on 
     * the version 113 of the pzunmqr_blgtrd.c
     * which go over the losange block_column 
     * by block column. but it is not important 
     * here the order because Ts are independent.
     * ========================================
    */ 
    nt  = plasma_ceildiv((N-1),Vblksiz);
    for (blkj=nt-1; blkj>=0; blkj--) {
        /* the index of the first row on the top of block (blkj) */ 
        firstrow = blkj * Vblksiz + 1;
        /*find the number of tile for this block */
        if( blkj == nt-1 )
            mt = plasma_ceildiv( N -  firstrow,    NB);
        else
            mt = plasma_ceildiv( N - (firstrow+1), NB);
        /*loop over the tiles find the size of the Vs and apply it */
        for (blki=mt; blki>0; blki--) {
            /*calculate the size of each losange of Vs= (Vm,Vn)*/
            myrow     = firstrow + (mt-blki)*NB;
            mycol     = blkj*Vblksiz;
            Vm = min( NB+Vblksiz-1, N-myrow);
            if( ( blkj == nt-1 ) && ( blki == mt ) ){
                Vn = min (Vblksiz, Vm);
            } else {
                Vn = min (Vblksiz, Vm-1);
            }
            /*calculate the pointer to the Vs and the Ts.
             * Note that Vs and Ts have special storage done
             * by the bulgechasing function*/
            findVTpos(N,NB,Vblksiz,mycol,myrow, &vpos, &taupos, &tpos, &blkid);
            myid = blkid/blkpercore;
            if( my_core_id==(myid%cores_num) ){
                if( ( Vm > 0 ) && ( Vn > 0 ) ){
                    LAPACKE_zlarft_work(LAPACK_COL_MAJOR, 
                                  lapack_const(PlasmaForward), 
                                  lapack_const(PlasmaColumnwise),
                                  Vm, Vn, V(vpos), LDV, TAU(taupos), T(tpos), LDT);
                }
            }
        }
    }
}
#undef V
#undef TAU
#undef T
/***************************************************************************/



