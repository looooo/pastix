/**
 *
 * @file cand.h
 *
 * PaStiX analyse headers for candidate array functions.
 *
 * @copyright 1998-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Pascal Henon
 * @date 2013-06-24
 *
 * @addtogroup blend_dev_cand
 * @{
 *    This module contains all subroutines to initialize the candidates array
 *    for each supernode, as well as supernode properties that are defined by
 *    level such as 2D layouts and 2D tasks.
 *
 **/
#ifndef _CAND_H_
#define _CAND_H_

/**
 * @brief Processor candidate group to own a column blok
 */
typedef struct cand_s {
    double       costlevel; /**< Cost from root to node                              */
    pastix_int_t treelevel; /**< Level of the cblk in the elimination tree (depth from the root) */
    pastix_int_t fcandnum;  /**< first processor number of this candidate group      */
    pastix_int_t lcandnum;  /**< last processor number of this candidate group       */
    pastix_int_t fccandnum; /**< first cluster number of the cluster candidate group */
    pastix_int_t lccandnum; /**< last cluster number of the cluster candidate group  */
    pastix_int_t cluster;   /**< Cluster id on which the task will be executed       */
    int8_t       cblktype;  /**< type of the distribution                            */
} Cand;

void candInit           (       Cand          *candtab,
                                pastix_int_t   cblknbr );

void candSetSubCandidate(       Cand          *candtab,
                          const EliminTree    *etree,
                                pastix_int_t   rootnum,
                                pastix_int_t   procnum );

int  candCheck          ( const Cand            *candtab,
                          const symbol_matrix_t *symbmtx );

void candSetClusterCand (       Cand          *candtab,
                                pastix_int_t   cblknbr,
                          const pastix_int_t  *core2clust,
                                pastix_int_t   coresnbr );

void candSave           ( const Cand          *candtab,
                                pastix_int_t   cblknbr,
                                char         **directory );

void candBuild          (       pastix_int_t           level_tasks2d,
                                pastix_int_t           width_tasks2d,
                                pastix_compress_when_t lr_when,
                                pastix_int_t           lr_width,
                                Cand            *candtab,
                                EliminTree      *etree,
                          const symbol_matrix_t *symbmtx,
                          const CostMatrix      *costmtx );
void candUpdate         ( Cand                  *candtab,
                          EliminTree            *etree,
                          const symbol_matrix_t *symbmtx,
                          const CostMatrix      *costmtx );


void candGenDot          ( const EliminTree *etree,
                           const Cand       *candtab,
                           FILE             *stream );
void candGenDotLevel     ( const EliminTree *etree,
                           const Cand       *candtab,
                           FILE             *stream,
                           pastix_int_t      level );
void candGenCompressedDot( const EliminTree *etree,
                           const Cand       *candtab,
                           FILE             *stream);

#endif

/**
 * @}
 */
