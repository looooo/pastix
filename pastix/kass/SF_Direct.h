/************************************************************/
/**                                                        **/
/**   NAME       :  SF_Direct.h                            **/
/**                                                        **/
/**   AUTHOR     : Pascal HENON                            **/
/**                                                        **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : Sept 2006       **/
/**                                                        **/
/**                                                        **/
/************************************************************/

/*
**  The function prototypes.
*/

#ifndef SF_DIRECT_
#define SF_DIRECT_
pastix_int_t SF_Direct(csptr A, pastix_int_t cblknbr, const pastix_int_t *rangtab, pastix_int_t *treetab, csptr P);
#endif
