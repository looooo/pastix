/**
 *
 * @file core_zplghe.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Piotr Luszczek
 * @author Pierre Lemarinier
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c
 *
 **/
#include "common.h"

#define COMPLEX
#undef REAL

/*
 Rnd64seed is a global variable but it doesn't spoil thread safety. All matrix
 generating threads only read Rnd64seed. It is safe to set Rnd64seed before
 and after any calls to create_tile(). The only problem can be caused if
 Rnd64seed is changed during the matrix generation time.
 */

//static unsigned long long int Rnd64seed = 100;
#define Rnd64_A 6364136223846793005ULL
#define Rnd64_C 1ULL
#define RndF_Mul 5.4210108624275222e-20f
#define RndD_Mul 5.4210108624275222e-20

#ifdef COMPLEX
#define NBELEM   2
#else
#define NBELEM   1
#endif

static unsigned long long int
Rnd64_jump(unsigned long long int n, unsigned long long int seed ) {
  unsigned long long int a_k, c_k, ran;
  int i;

  a_k = Rnd64_A;
  c_k = Rnd64_C;

  ran = seed;
  for (i = 0; n; n >>= 1, i++) {
    if (n & 1)
      ran = a_k * ran + c_k;
    c_k *= (a_k + 1);
    a_k *= a_k;
  }

  return ran;
}

#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zplghe = PCORE_zplghe
#define CORE_zplghe PCORE_zplghe
#endif
void CORE_zplghe( double bump, int m, int n, PLASMA_Complex64_t *A, int lda,
                  int bigM, int m0, int n0, unsigned long long int seed )
{
    PLASMA_Complex64_t *tmp = A;
    int64_t i, j;
    unsigned long long int ran, jump;

    jump = (unsigned long long int)m0 + (unsigned long long int)n0 * (unsigned long long int)bigM;

    /*
     * Tile diagonal
     */
    if ( m0 == n0 ) {
        for (j = 0; j < n; j++) {
            ran = Rnd64_jump( NBELEM * jump, seed );

            for (i = j; i < m; i++) {
                *tmp = 0.5f - ran * RndF_Mul;
                ran  = Rnd64_A * ran + Rnd64_C;
#ifdef COMPLEX
                *tmp += I*(0.5f - ran * RndF_Mul);
                ran   = Rnd64_A * ran + Rnd64_C;
#endif
                tmp++;
            }
            tmp  += (lda - i + j + 1);
            jump += bigM + 1;
        }

        for (j = 0; j < n; j++) {
#ifdef COMPLEX
            A[j+j*lda] += bump - I*cimag( A[j+j*lda] );
#else
            A[j+j*lda] += bump;
#endif

            for (i=0; i<j; i++) {
                A[lda*j+i] = conj( A[lda*i+j] );
            }
        }
    }
    /*
     * Lower part
     */
    else if ( m0 > n0 ) {
        for (j = 0; j < n; j++) {
            ran = Rnd64_jump( NBELEM * jump, seed );

            for (i = 0; i < m; i++) {
                *tmp = 0.5f - ran * RndF_Mul;
                ran  = Rnd64_A * ran + Rnd64_C;
#ifdef COMPLEX
                *tmp += I*(0.5f - ran * RndF_Mul);
                ran   = Rnd64_A * ran + Rnd64_C;
#endif
                tmp++;
            }
            tmp  += (lda - i);
            jump += bigM;
        }
    }
    /*
     * Upper part
     */
    else if ( m0 < n0 ) {
        /* Overwrite jump */
        jump = n0 + m0 * bigM;

        for (i = 0; i < m; i++) {
            ran = Rnd64_jump( NBELEM * jump, seed );

            for (j = 0; j < n; j++) {
                A[j*lda+i] = 0.5f - ran * RndF_Mul;
                ran = Rnd64_A * ran + Rnd64_C;
#ifdef COMPLEX
                A[j*lda+i] -= I*(0.5f - ran * RndF_Mul);
                ran = Rnd64_A * ran + Rnd64_C;
#endif
            }
            jump += bigM;
        }
    }
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zplghe( Quark *quark, Quark_Task_Flags *task_flags,
                        double bump, int m, int n, PLASMA_Complex64_t *A, int lda,
                        int bigM, int m0, int n0, unsigned long long int seed )
{
    DAG_CORE_PLGHE;
    QUARK_Insert_Task(quark, CORE_zplghe_quark, task_flags,
        sizeof(double),                   &bump, VALUE,
        sizeof(int),                      &m,    VALUE,
        sizeof(int),                      &n,    VALUE,
        sizeof(PLASMA_Complex64_t)*lda*n, A,         OUTPUT,
        sizeof(int),                      &lda,  VALUE,
        sizeof(int),                      &bigM, VALUE,
        sizeof(int),                      &m0,   VALUE,
        sizeof(int),                      &n0,   VALUE,
        sizeof(unsigned long long int),   &seed, VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zplghe_quark = PCORE_zplghe_quark
#define CORE_zplghe_quark PCORE_zplghe_quark
#endif
void CORE_zplghe_quark(Quark *quark)
{
    double bump;
    int m;
    int n;
    PLASMA_Complex64_t *A;
    int lda;
    int bigM;
    int m0;
    int n0;
    unsigned long long int seed;

    quark_unpack_args_9( quark, bump, m, n, A, lda, bigM, m0, n0, seed );
    CORE_zplghe( bump, m, n, A, lda, bigM, m0, n0, seed );
}

