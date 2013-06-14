/*
  File: debug_dump.h
  
  Functions to dump informations on disk.
*/

#ifndef DEBUG_DUMP_H
#define DEBUG_DUMP_H

/*
  Function: dump1
  
  Dumps ord->permtab on disk.

  Format: 
    > i -> permtab[i]

  parameters:
    ord    - Order structure to print permtab from.
    stream - FILE *, opened in write mode, in which permtab will be writen.
    colnbr - Number of elements in permtab.
 */
void dump1(Order *ord,
	   FILE  *stream, 
	   pastix_int_t    colnbr);


/*
  Function: dump2

  Prints internal CSCd, in (i,j,v) format, in a file.

  Parameters:
    datacode - SolverMatrix.
    stream   - FILE * opened in write mode.
  
 */
void dump2(const SolverMatrix * datacode,
           CscMatrix          * cscmtx,
	   pastix_float_t              * trandcsc,
	   FILE               *stream);


/*
  Function: dump3

  Prints solver matrix informations, in (i,j,v) format, in a file.

  Parameters:
    datacode - SolverMatrix.
    stream   - FILE * opened in write mode.
*/
void dump3(const SolverMatrix *datacode, 
	   FILE               *stream);

/*
  Function: dump3_LU

  Prints solver matrix informations, in (i,j,v) format, in a file, 
  for LU decomposition.

  Parameters:
    datacode - SolverMatrix.
    streamL  - FILE * opened in write mode.
    streamU  - FILE * opened in write mode.
*/
void dump3_LU(const SolverMatrix * datacode, 
	      FILE               * streamL, 
	      FILE               * streamU);


/*
  Function: dump4
  
  Writes column blocks and blocs dimension in a file.
  
  Parameters:
    datacode - SolverMatrix containing informations about blocs
    stream   - FILE * opened in write mode.
*/
void dump4(const SolverMatrix *datacode, 
	   FILE               *stream);


/*
  Function: dump5

  Writes right-hand-side memeber in a file.

  Parameters:
    datacode - SolverMatrix containing right-hand-side member.
    stream   - FILE * opened in write mode.
*/
void dump5(const SolverMatrix *datacode, 
	   FILE               *stream);


/*
  Function: dump6

  Prints diagonal blocks in the folowing format :
  > ** block diag <cblknbr> **
  > <line1> [<value1> <value2> ... ]
  > <line2> [...                   ]
  
  Prints one file dor L and one for U.

  Parameters:
    datacode - SolverMatrix.
    streamL  - FILE * into which L diagonal blocs will be writen.
    streamU  - FILE * into which U diagonal blocs will be writen.
*/
void dump6(const SolverMatrix *datacode, 
	   FILE               *streamL, 
	   FILE               *streamU);


/*
  Function: dump7
  
  Writes a vector in the folowing format :
  > <line1> <value[line1]>
  > <line2> <value[line1]>

  Parameters: 
    v      - vector to write.
    stream - FILE * opened in write mode.
    nbr    - Size of the vector v.
*/
void dump7(pastix_float_t *v, 
	   FILE  *stream, 
	   pastix_int_t    colnbr);

#endif /* DEBUG_DUMP_H */
