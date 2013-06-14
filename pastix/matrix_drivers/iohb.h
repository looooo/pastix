#ifndef IOHB_H
#define IOHB_H

#include<stdio.h>
#include<stdlib.h>
#ifdef X_ARCHi686_mac
#  include<malloc/malloc.h>
#else /* X_ARCHi686_mac */
#  ifdef __FreeBSD__
#    include            <stdlib.h>
#  else /* not __FreeBSD__ */
#    include            <malloc.h>
#  endif /* not __FreeBSD__ */
#endif /* X_ARCHi686_mac */

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __INTEL_COMPILER
/* Ignore icc remark : "operands are evaluated in unspecified order"*/
#pragma warning(disable:981)
/* Ignore icc remark : "external function definition with no prior declaration" */
#pragma warning(disable:1418)
/* Ignore icc remark : "external declaration in primary source file" */
#pragma warning(disable:1419)
/* Ignore icc remark : "parameter "arg" was never referenced" */
#pragma warning(disable:869)
/* Ignore icc remark : "floating-point equality and inequality comparisons are unreliable" */
#pragma warning(disable:1572)
#endif

int readHB_info(const char* filename, int* M, int* N, int* nz, char** Type, 
                                                      int* Nrhs);

int readHB_header(FILE* in_file, char* Title, char* Key, char* Type, 
                    int* Nrow, int* Ncol, int* Nnzero, int* Nrhs,
                    char* Ptrfmt, char* Indfmt, char* Valfmt, char* Rhsfmt, 
                    int* Ptrcrd, int* Indcrd, int* Valcrd, int* Rhscrd, 
                    char *Rhstype);

int readHB_mat_double(const char* filename, int colptr[], int rowind[], 
                                                                 double val[]);

int readHB_newmat_double(const char* filename, int* M, int* N, int* nonzeros, 
                         int** colptr, int** rowind, double** val);

int readHB_aux_double(const char* filename, const char AuxType, double b[]);

int readHB_newaux_double(const char* filename, const char AuxType, double** b);

int writeHB_mat_double(const char* filename, int M, int N, 
                        int nz, const int colptr[], const int rowind[], 
                        const double val[], int Nrhs, const double rhs[], 
                        const double guess[], const double exact[],
                        const char* Title, const char* Key, const char* Type, 
                        char* Ptrfmt, char* Indfmt, char* Valfmt, char* Rhsfmt,
                        const char* Rhstype);

int readHB_mat_char(const char* filename, int colptr[], int rowind[], 
                                           char val[], char* Valfmt);

int readHB_newmat_char(const char* filename, int* M, int* N, int* nonzeros, int** colptr, 
                          int** rowind, char** val, char** Valfmt);

int readHB_aux_char(const char* filename, const char AuxType, char b[]);

int readHB_newaux_char(const char* filename, const char AuxType, char** b, char** Rhsfmt);

int writeHB_mat_char(const char* filename, int M, int N, 
                        int nz, const int colptr[], const int rowind[], 
                        const char val[], int Nrhs, const char rhs[], 
                        const char guess[], const char exact[], 
                        const char* Title, const char* Key, const char* Type, 
                        char* Ptrfmt, char* Indfmt, char* Valfmt, char* Rhsfmt,
                        const char* Rhstype);

int ParseIfmt(char* fmt, int* perline, int* width);

int ParseRfmt(char* fmt, int* perline, int* width, int* prec, char* flag);

void IOHBTerminate(char* message);
#ifdef __cplusplus
}
#endif

#endif
