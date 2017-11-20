/**
 *
 * @file pastix_starpu_model.c
 *
 * PaStiX zgetrf StarPU wrapper.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Ian Masliah
 * @date 2017-06-24
 *
 * @{
 *
 **/

#include "common.h"
#include "solver.h"
#include "sopalin_data.h"
#include "pastix_zcores.h"
#include "pastix_starpu.h"
#include "codelets.h"
#if !defined(PASTIX_WITH_STARPU)
#error "This file should not be compiled if Starpu is not enabled"
#endif
#include <starpu.h>

double cblk_gemmsp_cost( struct starpu_task *task, struct starpu_perfmodel_arch *arch 
                        , unsigned nimpl)
{
  pastix_coefside_t sideA;
  pastix_coefside_t sideB;
  pastix_trans_t    trans;
  SolverCblk       *cblk;
  SolverBlok       *blok;
  SolverCblk       *fcblk;
  sopalin_data_t   *sopalin_data;

  starpu_codelet_unpack_args(task->cl_arg, &sideA, &sideB, &trans, &cblk, &blok, &fcblk, &sopalin_data);
  int shift = (sideA == PastixUCoef) ? 1 : 0;
  int K = cblk_colnbr( cblk );
  int N = blok_rownbr( blok );
  int M = cblk->stride - (cblk->cblktype & CBLK_LAYOUT_2D ? (blok + shift)->coefind / K : (blok + shift)->coefind);

  double cost = 0.;

  if (arch->devices->type == STARPU_CPU_WORKER)
  {
    cost = modelsGetCost3Param(&sopalin_data->cpu_coefs[PastixKernelGEMMCblk2d2d], M, N, K); 
  }
  else if(arch->devices->type == STARPU_CUDA_WORKER)
  {
    cost = modelsGetCost3Param(&sopalin_data->gpu_coefs[PastixKernelGEMMCblk2d2d], M, N, K);
  }

  (void)nimpl;

  return cost*10e6;
}

double blok_gemmsp_cost( struct starpu_task *task, struct starpu_perfmodel_arch *arch 
                               , unsigned nimpl)
{
  pastix_coefside_t sideA;
  pastix_coefside_t sideB;
  pastix_trans_t    trans;
  SolverCblk       *cblk;
  SolverBlok       *lblk;
  SolverBlok       *blokA,*blokB;
  SolverCblk       *fcblk;
  sopalin_data_t   *sopalin_data;
  pastix_int_t blok_mk, blok_nk, blok_mn;
  starpu_codelet_unpack_args(task->cl_arg, &sopalin_data,&sideA, &sideB, &trans, &cblk,  &fcblk, &blok_mk, &blok_nk, &blok_mn , &fcblk);
 
  int K = cblk_colnbr( cblk );
  int N = 0;
  int M = 0;

  blokA = cblk->fblokptr + blok_mk;
  blokB = cblk->fblokptr + blok_nk;
  lblk = cblk[1].fblokptr;

  M = blok_rownbr(blokA); 
  while( (blokA < lblk) &&
        (blokA[0].fcblknm == blokA[1].fcblknm) &&
        (blokA[0].lcblknm == blokA[1].lcblknm) )
  {
      blokA++;
      M += blok_rownbr(blokA);
  }
 
  N = blok_rownbr(blokB); 
  while( (blokB < lblk) &&
        (blokB[0].fcblknm == blokB[1].fcblknm) &&
        (blokB[0].lcblknm == blokB[1].lcblknm) )
  {
      blokB++;
      N += blok_rownbr(blokB);
  }

  double cost = 0.;

  if (arch->devices->type == STARPU_CPU_WORKER)
    cost =  modelsGetCost3Param(&sopalin_data->cpu_coefs[PastixKernelGEMMBlok2d2d],M,N,K);  
  else if(arch->devices->type == STARPU_CUDA_WORKER)
    cost = modelsGetCost3Param(&sopalin_data->gpu_coefs[PastixKernelGEMMBlok2d2d],M,N,K);

  (void)nimpl;

  return cost*10e6;
}

double cblk_getrfsp1d_cost( struct starpu_task *task, struct starpu_perfmodel_arch *arch 
                        , unsigned nimpl)
{
  sopalin_data_t *sopalin_data;
  SolverCblk *cblk;

  starpu_codelet_unpack_args(task->cl_arg, &cblk, &sopalin_data);

  double cost = 0.;

  if (arch->devices->type == STARPU_CPU_WORKER)
    cost =  modelsGetCost1Param(&sopalin_data->cpu_coefs[PastixKernelGETRF],cblk_colnbr(cblk));  
  else if(arch->devices->type == STARPU_CUDA_WORKER)
    cost = modelsGetCost1Param(&sopalin_data->gpu_coefs[PastixKernelGETRF],cblk_colnbr(cblk));

  (void)nimpl;

  return cost*10e6;
}

double blok_getrfsp1d_cost( struct starpu_task *task, struct starpu_perfmodel_arch *arch 
                               , unsigned nimpl)
{
  sopalin_data_t *sopalin_data;
  SolverCblk *cblk;

  starpu_codelet_unpack_args(task->cl_arg, &cblk, &sopalin_data);

  double cost = 0.;

  if (arch->devices->type == STARPU_CPU_WORKER)
    cost =  modelsGetCost1Param(&sopalin_data->cpu_coefs[PastixKernelGETRF],cblk_colnbr(cblk));  
  else if(arch->devices->type == STARPU_CUDA_WORKER)
    cost = modelsGetCost1Param(&sopalin_data->gpu_coefs[PastixKernelGETRF],cblk_colnbr(cblk));

  (void)nimpl;

  return cost*10e6;
}

double cblk_potrfsp1d_cost( struct starpu_task *task, struct starpu_perfmodel_arch *arch 
                        , unsigned nimpl)
{
  sopalin_data_t *sopalin_data;
  SolverCblk *cblk;

  starpu_codelet_unpack_args(task->cl_arg, &cblk, &sopalin_data);

  double cost = 0.;

  if (arch->devices->type == STARPU_CPU_WORKER)
    cost =  modelsGetCost1Param(&sopalin_data->cpu_coefs[PastixKernelPOTRF],cblk_colnbr(cblk));  
  else if(arch->devices->type == STARPU_CUDA_WORKER)
    cost = modelsGetCost1Param(&sopalin_data->gpu_coefs[PastixKernelPOTRF],cblk_colnbr(cblk));

  (void)nimpl;

  return cost*10e6;
}

double blok_potrfsp1d_cost( struct starpu_task *task, struct starpu_perfmodel_arch *arch 
                               , unsigned nimpl)
{
  sopalin_data_t *sopalin_data;
  SolverCblk *cblk;

  starpu_codelet_unpack_args(task->cl_arg, &cblk, &sopalin_data);

  double cost = 0.;

  if (arch->devices->type == STARPU_CPU_WORKER)
    cost =  modelsGetCost1Param(&sopalin_data->cpu_coefs[PastixKernelPOTRF],cblk_colnbr(cblk));  
  else if(arch->devices->type == STARPU_CUDA_WORKER)
    cost = modelsGetCost1Param(&sopalin_data->gpu_coefs[PastixKernelPOTRF],cblk_colnbr(cblk));

  (void)nimpl;

  return cost*10e6;
}

double blok_trsmsp_cost( struct starpu_task *task, struct starpu_perfmodel_arch *arch 
                               , unsigned nimpl)
{
  pastix_coefside_t coef;
  pastix_side_t     side;
  pastix_uplo_t     uplo;
  pastix_trans_t    trans;
  pastix_diag_t     diag;
  SolverCblk       *cblk;
  SolverBlok       *blok;
  SolverBlok       *lblk;
  pastix_int_t      blok_m;
  sopalin_data_t   *sopalin_data;

  starpu_codelet_unpack_args(task->cl_arg, &coef, &side, &uplo, &trans, &diag,
                             &cblk, &blok_m, &sopalin_data);

  int M = 0;
  blok = cblk->fblokptr + blok_m;
  lblk = cblk[1].fblokptr;

  M = blok_rownbr(blok); 
  while( (blok < lblk) &&
        (blok[0].fcblknm == blok[1].fcblknm) &&
        (blok[0].lcblknm == blok[1].lcblknm) )
  {
      blok++;
      M += blok_rownbr(blok);
  }
  double cost = 0.;

  if (arch->devices->type == STARPU_CPU_WORKER)
    cost =  modelsGetCost2Param(&sopalin_data->cpu_coefs[PastixKernelTRSMBlok2d],M,cblk_colnbr(cblk));  
  else if(arch->devices->type == STARPU_CUDA_WORKER)
    cost = modelsGetCost2Param(&sopalin_data->gpu_coefs[PastixKernelTRSMBlok2d],M,cblk_colnbr(cblk));

  (void)nimpl;

  return cost*10e6;
}
