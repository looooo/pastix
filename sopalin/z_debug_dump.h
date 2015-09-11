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
  File: z_debug_dump.h
  
  Functions to dump informations on disk.
*/

#ifndef DEBUG_DUMP_H
#define DEBUG_DUMP_H

/*
  Function: z_dump1
  
  Dumps ord->permtab on disk.

  Format: 
    > i -> permtab[i]

  parameters:
    ord    - Order structure to print permtab from.
    stream - FILE *, opened in write mode, in which permtab will be writen.
    colnbr - Number of elements in permtab.
 */
void z_dump1(Order *ord,
	   FILE  *stream, 
	   pastix_int_t    colnbr);


/*
  Function: z_dump2

  Prints internal CSCd, in (i,j,v) format, in a file.

  Parameters:
    datacode - z_SolverMatrix.
    stream   - FILE * opened in write mode.
  
 */
void z_dump2(const z_SolverMatrix * datacode,
           z_CscMatrix          * cscmtx,
	   pastix_complex64_t              * trandcsc,
	   FILE               *stream);


/*
  Function: z_dump3

  Prints z_solver matrix informations, in (i,j,v) format, in a file.

  Parameters:
    datacode - z_SolverMatrix.
    stream   - FILE * opened in write mode.
*/
void z_dump3(const z_SolverMatrix *datacode, 
	   FILE               *stream);

/*
  Function: z_dump3_LU

  Prints z_solver matrix informations, in (i,j,v) format, in a file, 
  for LU decomposition.

  Parameters:
    datacode - z_SolverMatrix.
    streamL  - FILE * opened in write mode.
    streamU  - FILE * opened in write mode.
*/
void z_dump3_LU(const z_SolverMatrix * datacode, 
	      FILE               * streamL, 
	      FILE               * streamU);


/*
  Function: z_dump4
  
  Writes column blocks and blocs dimension in a file.
  
  Parameters:
    datacode - z_SolverMatrix containing informations about blocs
    stream   - FILE * opened in write mode.
*/
void z_dump4(const z_SolverMatrix *datacode, 
	   FILE               *stream);


/*
  Function: z_dump5

  Writes right-hand-side memeber in a file.

  Parameters:
    datacode - z_SolverMatrix containing right-hand-side member.
    stream   - FILE * opened in write mode.
*/
void z_dump5(const z_SolverMatrix *datacode, 
	   FILE               *stream);


/*
  Function: z_dump6

  Prints diagonal blocks in the folowing format :
  > ** block diag <cblknbr> **
  > <line1> [<value1> <value2> ... ]
  > <line2> [...                   ]
  
  Prints one file dor L and one for U.

  Parameters:
    datacode - z_SolverMatrix.
    streamL  - FILE * into which L diagonal blocs will be writen.
    streamU  - FILE * into which U diagonal blocs will be writen.
*/
void z_dump6(const z_SolverMatrix *datacode, 
	   FILE               *streamL, 
	   FILE               *streamU);


/*
  Function: z_dump7
  
  Writes a vector in the folowing format :
  > <line1> <value[line1]>
  > <line2> <value[line1]>

  Parameters: 
    v      - vector to write.
    stream - FILE * opened in write mode.
    nbr    - Size of the vector v.
*/
void z_dump7(pastix_complex64_t *v, 
	   FILE  *stream, 
	   pastix_int_t    colnbr);

#endif /* DEBUG_DUMP_H */
