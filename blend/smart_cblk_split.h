#ifndef SMART_CBLK_SPLIT_H
#define SMART_CBLK_SPLIT_H

int smart_cblk_split(const BlendCtrl      * ctrl,
		     const SymbolMatrix   * symbmtx, 
		     pastix_int_t       cblknum,
		     pastix_int_t       procnbr,
		     pastix_int_t       blas_min_col,
		     pastix_int_t       blas_max_col,
		     pastix_int_t     * nseq,
		     pastix_int_t    ** seq);

#endif /* SMART_CBLK_SPLIT_H */
