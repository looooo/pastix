/*
  File: csc_intern_io.h

  Functions to save or load internal CSC in binary or ascii mode.
  
*/

#ifndef CSC_INTERN_IO_H
#define CSC_INTERN_IO_H
/*
  Function: CscSave

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
PASTIX_INT CscSave(const CscMatrix * const cscptr, 
	    FILE            * const stream);

/*
  Function: CscBSave

  Writes on disk an internal CSCd in binary format.

  Parameters :
    cscprt - the internal CSCd structure to save.
    stream - the FILE to write into, open in write mode.
*/
PASTIX_INT CscBSave(const CscMatrix * const cscptr, 
	     FILE            * const stream);

/* 
   Function: CscLoad

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
PASTIX_INT CscLoad(CscMatrix * cscptr, 
	    FILE      * stream);

/*
  Function: CscBLoad
  
  Loads an internal CSCd from a file saved in binary mode.

  Parameters :
    cscprt - the internal CSCd structure to load.
    stream - the FILE to write into, open in read mode. 
*/
PASTIX_INT CscBLoad(CscMatrix * cscptr, 
	     FILE      * stream);

#endif /* CSC_INTERN_IO_H */
