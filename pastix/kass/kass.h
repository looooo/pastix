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

void kass(int            levelk, 
	  int            rat, 
	  SymbolMatrix * symbptr, 
	  pastix_int_t            baseval,
	  pastix_int_t            vertnbr, 
	  pastix_int_t            edgenbr, 
	  pastix_int_t          * verttab,
	  pastix_int_t          * edgetab, 
	  Order        * orderptr, 
	  MPI_Comm       pastix_comm);

int kass2(int            ilu,
          int            levelk,
          int            rat,
          SymbolMatrix * symbmtx,
          pastix_csc_t * csc,
          Order        * orderptr,
          MPI_Comm       pastix_comm);

/* void kass(int alpha, int rat, SymbolMatrix * symbptr, Graph * graphptr, Order * orderptr, MPI_Comm pastix_comm); */

void ifax(pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_int_t levelk, pastix_int_t  cblknbr, pastix_int_t *rangtab, pastix_int_t *perm, pastix_int_t *iperm, SymbolMatrix *symbmtx);

