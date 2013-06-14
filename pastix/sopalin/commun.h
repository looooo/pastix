#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <inttypes.h>
#include <assert.h>

#include "common_pastix.h"
#include "sopalin_compute.h"

#define SIZEEND  128
#define SIZEPAS  4
#define SIZEINIT 4

#define ITER 100

#define SIZEMAX 1000

pastix_float_t vecbuf1[SIZEMAX];
pastix_float_t vecbuf2[SIZEMAX];
pastix_float_t matbuf1[SIZEMAX][SIZEMAX];
pastix_float_t matbuf2[SIZEMAX][SIZEMAX];
pastix_float_t matbuf3[SIZEMAX][SIZEMAX];
pastix_float_t matppf1[SIZEMAX*SIZEMAX];
pastix_float_t matppf2[SIZEMAX*SIZEMAX];

void init_commun();
void init_vecteur(pastix_int_t n,pastix_float_t *v,pastix_float_t new[SIZEMAX],pastix_int_t inc);
void end_vecteur(pastix_int_t n,pastix_float_t *v,pastix_float_t new[SIZEMAX],pastix_int_t inc);
void init_matrice(char *t,pastix_int_t n,pastix_int_t m,pastix_float_t *a,pastix_float_t new[SIZEMAX][SIZEMAX],pastix_int_t lda);
void end_matrice(char *t,pastix_int_t n,pastix_int_t m,pastix_float_t *a,pastix_float_t new[SIZEMAX][SIZEMAX],pastix_int_t lda);

void init_vecteur(pastix_int_t n,pastix_float_t *v,pastix_float_t new[SIZEMAX],pastix_int_t inc)
{
  pastix_int_t i;
  for (i=0;i<n;i++)
    new[i]=v[i*inc];
}

void end_vecteur(pastix_int_t n,pastix_float_t *v,pastix_float_t new[SIZEMAX],pastix_int_t inc)
{
  pastix_int_t i;
  for (i=0;i<n;i++)
    v[i*inc]=new[i];
}

void init_matrice(char *t,pastix_int_t n,pastix_int_t m,pastix_float_t *a,pastix_float_t new[SIZEMAX][SIZEMAX],pastix_int_t lda)
{
  pastix_int_t i,j;
  for (i=0;i<n;i++)
    for (j=0;j<m;j++)
      if (*t=='N')
	new[j][i]=a[i+j*lda];
      else
	new[i][j]=a[i+j*lda];
}

void end_matrice(char *t,pastix_int_t n,pastix_int_t m,pastix_float_t *a,pastix_float_t new[SIZEMAX][SIZEMAX],pastix_int_t lda)
{
  pastix_int_t i,j;
  for (i=0;i<n;i++)
    for (j=0;j<m;j++)
      if (*t=='N')
	a[i+j*lda]=new[j][i];
      else
	a[i+j*lda]=new[i][j];
}

static Clock test_clk;

#define TEST_CLOCK_INIT {clockInit(&test_clk);clockStart(&test_clk);}
#define TEST_CLOCK_STOP {clockStop(&test_clk);}
#define TEST_CLOCK_GET  clockVal(&test_clk)
