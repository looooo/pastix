#ifndef SMART_CBLK_SPLIT_H
#define SMART_CBLK_SPLIT_H

int smart_cblk_split(BlendCtrl      * ctrl,
		     SymbolMatrix   * symbmtx, 
		     PASTIX_INT       cblknum,
		     PASTIX_INT       procnbr,
		     PASTIX_INT       blas_min_col,
		     PASTIX_INT       blas_max_col,
		     PASTIX_INT     * nseq,
		     PASTIX_INT    ** seq);

#endif /* SMART_CBLK_SPLIT_H */
