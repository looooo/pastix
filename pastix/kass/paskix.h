/************************************************************/
/**                                                        **/
/**   NAME       :  paskix.h                              **/
/**                                                        **/
/**   AUTHOR     : Pascal HENON                            **/
/**                                                        **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 28 Nov 2006     **/
/**                                                        **/
/**  Some functions on CSC matrices                        **/
/**                                                        **/
/************************************************************/

/*
**  The function prototypes.
*/
void CSC_Fnum2Cnum(pastix_int_t *ja, pastix_int_t *ia, pastix_int_t n);
void CSC_Cnum2Fnum(pastix_int_t *ja, pastix_int_t *ia, pastix_int_t n);

void CSC_rowPerm(pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_float_t *a, pastix_int_t *rperm);
void CSC_colPerm(pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_float_t *a, pastix_int_t *cperm);
void CSC_colScale(pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_float_t *a, pastix_float_t *dcol);
void CSC_rowScale(pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_float_t *a, pastix_float_t *drow);
void CSC_sort(pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_float_t *a);
