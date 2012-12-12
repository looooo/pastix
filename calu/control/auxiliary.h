/**
 *
 * @file control/auxiliary.h
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Jakub Kurzak
 * @author Piotr Luszczek
 * @author Emmanuel Agullo
 * @date 2010-11-15
 *
 **/
#ifndef _PLASMA_AUXILIARY_H_
#define _PLASMA_AUXILIARY_H_

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  Internal routines
 **/
void plasma_warning(const char *func_name, char* msg_text);
void plasma_error(const char *func_name, char* msg_text);
void plasma_fatal_error(const char *func_name, char* msg_text);
int plasma_element_size(int type);
void plasma_memcpy(void *dst, void *src, PLASMA_size size, int type);
void plasma_memzero(void *memptr, PLASMA_size size, int type);
void plasma_memset_int(int *mem, int size, int value);
int plasma_rank(plasma_context_t *plasma);
int plasma_tune(PLASMA_enum func, int M, int N, int NRHS);

/***************************************************************************//**
 *  User routines
 **/
int PLASMA_Version(int *ver_major, int *ver_minor, int *ver_micro);

#ifdef __cplusplus
}
#endif

#endif
