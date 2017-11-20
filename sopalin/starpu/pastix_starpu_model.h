/**
 *
 * @file starpu_model.h
 *
 * Pastix StarPU model function
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Ian Masliah
 * @date 2013-06-24
 *
 **/

#ifndef _pastix_starpu_model_h_
#define _pastix_starpu_model_h_

double cblk_gemmsp_cost( struct starpu_task *task, struct starpu_perfmodel_arch *arch 
                        , unsigned nimpl);

double blok_gemmsp_cost( struct starpu_task *task, struct starpu_perfmodel_arch *arch 
                               , unsigned nimpl);

double cblk_getrfsp1d_cost( struct starpu_task *task, struct starpu_perfmodel_arch *arch 
                        , unsigned nimpl);

double blok_getrfsp1d_cost( struct starpu_task *task, struct starpu_perfmodel_arch *arch 
                               , unsigned nimpl);

double cblk_potrfsp1d_cost( struct starpu_task *task, struct starpu_perfmodel_arch *arch 
                        , unsigned nimpl);

double blok_potrfsp1d_cost( struct starpu_task *task, struct starpu_perfmodel_arch *arch 
                               , unsigned nimpl);

double blok_trsmsp_cost( struct starpu_task *task, struct starpu_perfmodel_arch *arch 
                               , unsigned nimpl);


#endif /* _pastix_starpu_model_h_ */