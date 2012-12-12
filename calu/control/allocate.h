/**
 *
 * @file allocate.h
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Jakub Kurzak
 * @date 2010-11-15
 *
 **/
#ifndef _PLASMA_ALLOCATE_H_
#define _PLASMA_ALLOCATE_H_

#ifdef __cplusplus
extern "C" {
#endif

void *plasma_shared_alloc(plasma_context_t *plasma, size_t size, int type);
void plasma_shared_free(plasma_context_t *plasma, void *ptr);
void *plasma_private_alloc(plasma_context_t *plasma, size_t size, int type);
void plasma_private_free(plasma_context_t *plasma, void *ptr);

#ifdef __cplusplus
}
#endif

#endif
