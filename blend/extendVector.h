/************************************************************/
/**                                                        **/
/**   NAME       : extendVector.h                          **/
/**                                                        **/
/**   AUTHOR     : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Vector that can extend its size         **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 27 jul 1998     **/
/**                                 to     03 sep 1998     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

typedef struct ExtendVectorINT_ {
    pastix_int_t          vecsize;
    pastix_int_t          eltnbr;          /*+ number of elements +*/
    pastix_int_t    *     inttab;          /*+ array of pastix_int_t       +*/
} ExtendVectorINT;


/*
**  The function prototypes.
*/

#ifndef EXTENDVECTOR
#define static
#endif

pastix_int_t                     *extendint_Init    (ExtendVectorINT *, pastix_int_t);
void                     extendint_Exit    (ExtendVectorINT *);
void                     extendint_Add     (ExtendVectorINT *, pastix_int_t);
pastix_int_t                      extendint_Size    (ExtendVectorINT *);
pastix_int_t                      extendint_Read    (ExtendVectorINT *, pastix_int_t);
void                     extendint_Clear   (ExtendVectorINT *);
void                     extendint_ToSize  (pastix_int_t, ExtendVectorINT *);
void                     extendint_incr    (ExtendVectorINT *);
#undef static
