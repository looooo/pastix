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
    PASTIX_INT          vecsize;
    PASTIX_INT          eltnbr;          /*+ number of elements +*/
    PASTIX_INT    *     inttab;          /*+ array of PASTIX_INT       +*/
} ExtendVectorINT;


/*
**  The function prototypes.
*/

#ifndef EXTENDVECTOR
#define static
#endif

PASTIX_INT                     *extendint_Init    (ExtendVectorINT *, PASTIX_INT);
void                     extendint_Exit    (ExtendVectorINT *);
void                     extendint_Add     (ExtendVectorINT *, PASTIX_INT);
PASTIX_INT                      extendint_Size    (ExtendVectorINT *);
PASTIX_INT                      extendint_Read    (ExtendVectorINT *, PASTIX_INT);
void                     extendint_Clear   (ExtendVectorINT *);
void                     extendint_ToSize  (PASTIX_INT, ExtendVectorINT *);
void                     extendint_incr    (ExtendVectorINT *);
#undef static
