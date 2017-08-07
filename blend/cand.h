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

int  candCheck          ( const Cand          *candtab,
                          const SymbolMatrix  *symbmtx );

void candSetClusterCand (       Cand          *candtab,
                                pastix_int_t   cblknbr,
                          const pastix_int_t  *core2clust,
                                pastix_int_t   coresnbr );

void candSave           (       pastix_data_t *pastix,
                          const Cand          *candtab,
                                pastix_int_t   cblknbr );

void candBuild          (       pastix_int_t   autolevel,
                                pastix_int_t   level2D,
                                pastix_int_t   ratiolimit,
                                Cand          *candtab,
                                EliminTree    *etree,
                          const SymbolMatrix  *symbmtx,
                          const CostMatrix    *costmtx );

#endif

/**
 * @}
 */
