/************************************************************/
/**                                                        **/
/**   NAME       : kass.h                                  **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Compute a block structure of the factor **/
/**                obtained by a ILU(k) factorization      **/
/**                                                        **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 30/01/2006      **/
/**                                 to                     **/
/**                                                        **/
/************************************************************/
#ifndef _KASS_H_
#define _KASS_H_

struct kass_csr_s {
    pastix_int_t   n;
    pastix_int_t  *nnz;
    pastix_int_t **rows;
};
typedef struct kass_csr_s kass_csr_t;

void kass_csrInit( pastix_int_t n, kass_csr_t *csr );
void kass_csrClean( kass_csr_t *csr );
pastix_int_t kass_csrGetNNZ( kass_csr_t *csr );
int  kass_csrGenPA( const pastix_graph_t *graphA,
                    const pastix_int_t   *perm,
                          kass_csr_t     *graphPA );
void kass_csrCompact( kass_csr_t *csr );

void
kassBuildSymbol(      kass_csr_t   *P,
                      pastix_int_t  cblknbr,
                const pastix_int_t *rangtab,
                      SymbolMatrix *symbmtx);
void
kassPatchSymbol( SymbolMatrix *symbmtx );

pastix_int_t SF_Direct(const kass_csr_t   *graphA,
                             pastix_int_t  cblknbr,
                       const pastix_int_t *rangtab,
                             pastix_int_t *treetab,
                             kass_csr_t   *graphL );

void amalgamate2(double rat,
                 kass_csr_t    *graphL,
                 pastix_int_t   nnzL,
                 pastix_int_t   snodenbr,
                 pastix_int_t  *snodetab,
                 pastix_int_t  *treetab,
                 pastix_int_t  *cblknbr,
                 pastix_int_t **rangtab,
                 pastix_int_t  *nodetab,
                 MPI_Comm pastix_comm );

int kass2(int             ilu,
          int             levelk,
          int             rat,
          SymbolMatrix   *symbmtx,
          pastix_graph_t *graph,
          Order          *orderptr,
          MPI_Comm        pastix_comm);

void ifax(pastix_int_t n,
          pastix_int_t *ia,
          pastix_int_t *ja,
          pastix_int_t levelk,
          pastix_int_t  cblknbr,
          pastix_int_t *rangtab,
          pastix_int_t *perm,
          pastix_int_t *iperm,
          SymbolMatrix *symbmtx);

#endif
