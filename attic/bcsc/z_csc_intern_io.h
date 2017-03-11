/**
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/
/*
  File: z_csc_intern_io.h

  Functions to save or load internal CSC in binary or ascii mode.
  
*/

#ifndef Z_CSC_INTERN_IO_H
#define Z_CSC_INTERN_IO_H
/*
  Function: z_CscSave

  Writes on disk an internal CSCd in text format.

  Format is :
  
  > CSC_FNBR(cscptr)
  > CSC_COLNBR(cscptr,iter)    ! iter = 0 to CSC_FNBR(cscptr) - 1
  > CSC_COL(cscptr,iter,iter2) ! iter2 = 0 to CSC_COLNBR(cscptr,iter)
  > ...
  > CSC_ROW(cscptr,iter) ! For all rows and values (iter)
  > CSC_VAL(cscptr,iter)

  Parameters :
    cscprt - the internal CSCd structure to save.
    stream - the FILE to write into, open in write mode.
*/
pastix_int_t z_CscSave(const z_CscMatrix * const cscptr,
                       FILE            * const stream);

/*
  Function: z_CscBSave

  Writes on disk an internal CSCd in binary format.

  Parameters :
    cscprt - the internal CSCd structure to save.
    stream - the FILE to write into, open in write mode.
*/
pastix_int_t z_CscBSave(const z_CscMatrix * const cscptr, 
	     FILE            * const stream);

/* 
   Function: z_CscLoad

   Reads an internal CSCd from disk.

   Format is :
   
   > CSC_FNBR(cscptr)
   > CSC_COLNBR(cscptr,iter)    ! iter = 0 to CSC_FNBR(cscptr) - 1
   > CSC_COL(cscptr,iter,iter2) ! iter2 = 0 to CSC_COLNBR(cscptr,iter)
   > ...
   > CSC_ROW(cscptr,iter) ! For all rows and values (iter)
   > CSC_VAL(cscptr,iter)

   Parameters :
     cscprt - the internal CSCd structure to load.
     stream - the FILE to write into, open in read mode. 
*/
pastix_int_t z_CscLoad(z_CscMatrix * cscptr, 
	    FILE      * stream);

/*
  Function: z_CscBLoad
  
  Loads an internal CSCd from a file saved in binary mode.

  Parameters :
    cscprt - the internal CSCd structure to load.
    stream - the FILE to write into, open in read mode. 
*/
pastix_int_t z_CscBLoad(z_CscMatrix * cscptr, 
	     FILE      * stream);

#endif /* Z_CSC_INTERN_IO_H */
