/**
 *
 * @file plasmawinthread.h
 *
 *  This file handles the mapping from pthreads calls to windows threads
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Piotr Luszczek
 * @date 2010-11-15
 *
 **/
#ifndef PLASMWINTHREAD_H
#define PLASMWINTHREAD_H

#include <windows.h>

/*
typedef struct pthread_s {
  HANDLE Hth;
  unsigned IDth;
  void *(*Fth) (void *);
} pthread_t;
*/
typedef struct pthread_s {
  HANDLE hThread;
  unsigned int uThId;
} pthread_t;

typedef HANDLE pthread_mutex_t;
typedef int pthread_mutexattr_t;
typedef int pthread_attr_t;
typedef int pthread_condattr_t;

typedef struct pthread_cond_s {
  HANDLE hSem;
  HANDLE hEvt;
  CRITICAL_SECTION cs;
  int waitCount; /* waiting thread counter */
} pthread_cond_t;

typedef int pthread_attr_t;

#define PTHREAD_MUTEX_INITIALIZER ((pthread_mutex_t) -1)

#define PTHREAD_SCOPE_SYSTEM 1

#define PLASMA_DLLPORT
#define PLASMA_CDECL __cdecl

PLASMA_DLLPORT pthread_t PLASMA_CDECL pthread_self(void);
PLASMA_DLLPORT int PLASMA_CDECL pthread_mutex_init(pthread_mutex_t *mutex, const pthread_mutexattr_t * attr);
PLASMA_DLLPORT int PLASMA_CDECL pthread_mutex_destroy(pthread_mutex_t *mutex);
PLASMA_DLLPORT int PLASMA_CDECL pthread_mutex_lock(pthread_mutex_t *mutex);
PLASMA_DLLPORT int PLASMA_CDECL pthread_mutex_trylock(pthread_mutex_t *mutex);
PLASMA_DLLPORT int PLASMA_CDECL pthread_mutex_unlock(pthread_mutex_t *mutex);
PLASMA_DLLPORT int PLASMA_CDECL pthread_attr_init(pthread_attr_t *attr);
PLASMA_DLLPORT int PLASMA_CDECL pthread_attr_destroy(pthread_attr_t *attr);
PLASMA_DLLPORT int PLASMA_CDECL pthread_attr_setscope(pthread_attr_t *attr, int scope);
PLASMA_DLLPORT int PLASMA_CDECL pthread_create(pthread_t *tid, const pthread_attr_t *attr, void *(*start) (void *), void *arg);
PLASMA_DLLPORT int PLASMA_CDECL pthread_cond_init(pthread_cond_t *cond, const pthread_condattr_t *attr);
PLASMA_DLLPORT int PLASMA_CDECL pthread_cond_destroy(pthread_cond_t *cond);
PLASMA_DLLPORT int PLASMA_CDECL pthread_cond_wait(pthread_cond_t *cond, pthread_mutex_t *mutex);
PLASMA_DLLPORT int PLASMA_CDECL pthread_cond_broadcast(pthread_cond_t *cond);
PLASMA_DLLPORT int PLASMA_CDECL pthread_join(pthread_t thread, void **value_ptr);
PLASMA_DLLPORT int PLASMA_CDECL pthread_equal(pthread_t thread1, pthread_t thread2);

PLASMA_DLLPORT int PLASMA_CDECL pthread_setconcurrency (int);

PLASMA_DLLPORT unsigned int PLASMA_CDECL pthread_self_id(void);

#endif
