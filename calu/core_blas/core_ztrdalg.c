/**
 *
 * @file core_ztrdalg.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Azzam Haidar
 * @date 2011-05-15
 * @precisions normal z -> c d s
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 *  CORE_ztrdalg is a part of the tridiagonal reduction algorithm (bulgechasing)
 *  It correspond to a local driver of the kernels that should be executed on a
 *  single core.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *         @arg PlasmaLower:
 *         @arg PlasmaUpper:
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[in] NB
 *          The size of the Bandwidth of the matrix A,
 *          which correspond to the tile size. NB >= 0.
 *
 * @param[in] pA
 *          A pointer to the descriptor of the matrix A.
 *
 * @param[out] V
 *          PLASMA_Complex64_t array, dimension (N).
 *          The scalar elementary reflectors are written in this
 *          array. So it is used as a workspace for V at each step
 *          of the bulge chasing algorithm.
 *
 * @param[out] TAU
 *          PLASMA_Complex64_t array, dimension (N).
 *          The scalar factors of the elementary reflectors are written
 *          in thisarray. So it is used as a workspace for TAU at each step
 *          of the bulge chasing algorithm.
 *
 * @param[in] i
 *          Integer that refer to the current sweep. (outer loop).
 *
 * @param[in] j
 *          Integer that refer to the sweep to chase.(inner loop).
 *
 * @param[in] m
 *          Integer that refer to a sweep step, to ensure order dependencies.
 *
 * @param[in] grsiz
 *          Integer that refer to the size of a group.
 *          group mean the number of kernel that should be executed sequentially
 *          on the same core.
 *          group size is a trade-off between locality (cache reuse) and parallelism.
 *          a small group size increase parallelism while a large group size increase
 *          cache reuse.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ztrdalg = PCORE_ztrdalg
#define CORE_ztrdalg PCORE_ztrdalg
#endif
void CORE_ztrdalg(PLASMA_enum        uplo,
                        int n,
                        int nb,
                        PLASMA_Complex64_t *A,
                        int lda,
                        PLASMA_Complex64_t *V,
                        PLASMA_Complex64_t *TAU,
                        int Vblksiz, int wantz, 
                        int grsiz, int lcsweep, int id, int blksweep,
                        PLASMA_Complex64_t *work)
{
  int i, blkid, st, ed, sweepid;
  int nbtiles = plasma_ceildiv(n,nb);
  int KDM1  =  nb-1;


  /* code for all tiles */
  for (i = 0; i < grsiz ; i++) {
     blkid = id+i;
     st    = min(blkid*nb+lcsweep+1, n-1);
     ed    = min(st+KDM1, n-1);
     sweepid = blksweep*nb + lcsweep;
     /*printf("  COUCOU voici n %5d   nb %5d  sweepid %d   st %5d  ed %5d lcsweep %5d   id %5d   blkid %5d\n",n, nb,sweepid, st, ed,lcsweep, id, blkid);*/


     if((st == ed) && (sweepid != n-2) ) /* quick return in case of last tile */
        return;

     if(blkid==blksweep){
        CORE_zhbtype1cb(n, nb, A, lda, V, TAU, st, ed, sweepid, Vblksiz, wantz, work);     
     }else{
        CORE_zhbtype3cb(n, nb, A, lda, V, TAU, st, ed, sweepid, Vblksiz, wantz, work);     
     }
     /*printf("n %d  nb %d sweepid %d  myid %d   st %d  ed %d \n",n,nb,sweepid,blkid*2+1,st,ed);*/

     if(id!=(nbtiles-1)) {
         CORE_zhbtype2cb(n, nb, A, lda, V, TAU, st, ed, sweepid, Vblksiz, wantz, work);
         /*printf("n %d  nb %d sweepid %d  myid %d   st %d  ed %d \n",n,nb,sweepid,blkid*2+1+1,st,ed);*/
     }


  }

}

/***************************************************************************//**
 *
 **/
#define A(m_)    (A + (lda * nb * (m_))) 
void QUARK_CORE_ztrdalg(Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum uplo,
                        int n,
                        int nb,
                        PLASMA_Complex64_t *A,
                        int lda,
                        PLASMA_Complex64_t *V,
                        PLASMA_Complex64_t *TAU,
                        int Vblksiz, int wantz, 
                        int grsiz, int lcsweep, int id, int blksweep)
{
  Quark_Task *MYTASK;
  int ii, cur_id,  nbtiles    = plasma_ceildiv(n,nb);

           /*printf("coucou from quark function id %d    lcsweep %d    blksweep %d   grsiz %d  nbtiles %d\n", id, lcsweep, blksweep, grsiz, nbtiles);*/
           MYTASK = QUARK_Task_Init( quark, CORE_ztrdalg_quark,   task_flags);
           QUARK_Task_Pack_Arg(quark, MYTASK,  sizeof(PLASMA_enum),          &uplo,               VALUE );
           QUARK_Task_Pack_Arg(quark, MYTASK,  sizeof(int),                     &n,               VALUE );
           QUARK_Task_Pack_Arg(quark, MYTASK,  sizeof(int),                    &nb,               VALUE );
           QUARK_Task_Pack_Arg(quark, MYTASK,  sizeof(PLASMA_Complex64_t),       A,               NODEP );
           QUARK_Task_Pack_Arg(quark, MYTASK,  sizeof(int),                   &lda,               VALUE );
           QUARK_Task_Pack_Arg(quark, MYTASK,  sizeof(PLASMA_Complex64_t),       V,               NODEP );
           QUARK_Task_Pack_Arg(quark, MYTASK,  sizeof(PLASMA_Complex64_t),     TAU,               NODEP );
           QUARK_Task_Pack_Arg(quark, MYTASK,  sizeof(int),               &Vblksiz,               VALUE );
           QUARK_Task_Pack_Arg(quark, MYTASK,  sizeof(int),                 &wantz,               VALUE );
           QUARK_Task_Pack_Arg(quark, MYTASK,  sizeof(int),                 &grsiz,               VALUE );
           QUARK_Task_Pack_Arg(quark, MYTASK,  sizeof(int),               &lcsweep,               VALUE );
           QUARK_Task_Pack_Arg(quark, MYTASK,  sizeof(int),                    &id,               VALUE );
           QUARK_Task_Pack_Arg(quark, MYTASK,  sizeof(int),               &blksweep,              VALUE );
           QUARK_Task_Pack_Arg(quark, MYTASK,  sizeof(PLASMA_Complex64_t)*nb,     NULL,           SCRATCH);

           QUARK_Task_Pack_Arg(quark, MYTASK, sizeof(PLASMA_Complex64_t),    A(id),           INOUT );
           if( id<(nbtiles-1) )
              QUARK_Task_Pack_Arg(quark, MYTASK, sizeof(PLASMA_Complex64_t),    A(id+1),           INOUT );

           cur_id = id;
           for (ii = 1; ii < grsiz ; ii++) {
                cur_id = cur_id+1;
                if( cur_id<(nbtiles-1) )
                   QUARK_Task_Pack_Arg(quark, MYTASK, sizeof(PLASMA_Complex64_t),    A(cur_id+1),           INOUT );
           }

           QUARK_Insert_Task_Packed(quark, MYTASK);
}
#undef A
/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ztrdalg_quark = PCORE_ztrdalg_quark
#define CORE_ztrdalg_quark PCORE_ztrdalg_quark
#endif
void CORE_ztrdalg_quark(Quark *quark)
{
    PLASMA_enum         uplo;
    int                    n;
    int                   nb;
    PLASMA_Complex64_t    *A;
    int                  lda;
    PLASMA_Complex64_t    *V;
    PLASMA_Complex64_t  *TAU;
    int              Vblksiz;
    int                wantz;
    int                grsiz;
    int              lcsweep;
    int                   id;
    int             blksweep;
    PLASMA_Complex64_t *work;

    quark_unpack_args_14(quark, uplo, n, nb, A, lda, V, TAU, Vblksiz, wantz, grsiz, lcsweep, id, blksweep, work);
    CORE_ztrdalg(uplo, n, nb, A, lda, V, TAU, Vblksiz, wantz, grsiz, lcsweep, id, blksweep, work);
}
