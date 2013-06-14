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

PASTIX_FLOAT vecbuf1[SIZEMAX];
PASTIX_FLOAT vecbuf2[SIZEMAX];
PASTIX_FLOAT matbuf1[SIZEMAX][SIZEMAX];
PASTIX_FLOAT matbuf2[SIZEMAX][SIZEMAX];
PASTIX_FLOAT matbuf3[SIZEMAX][SIZEMAX];
PASTIX_FLOAT matppf1[SIZEMAX*SIZEMAX];
PASTIX_FLOAT matppf2[SIZEMAX*SIZEMAX];

void init_commun();
void init_vecteur(PASTIX_INT n,PASTIX_FLOAT *v,PASTIX_FLOAT new[SIZEMAX],PASTIX_INT inc);
void end_vecteur(PASTIX_INT n,PASTIX_FLOAT *v,PASTIX_FLOAT new[SIZEMAX],PASTIX_INT inc);
void init_matrice(char *t,PASTIX_INT n,PASTIX_INT m,PASTIX_FLOAT *a,PASTIX_FLOAT new[SIZEMAX][SIZEMAX],PASTIX_INT lda);
void end_matrice(char *t,PASTIX_INT n,PASTIX_INT m,PASTIX_FLOAT *a,PASTIX_FLOAT new[SIZEMAX][SIZEMAX],PASTIX_INT lda);

void init_vecteur(PASTIX_INT n,PASTIX_FLOAT *v,PASTIX_FLOAT new[SIZEMAX],PASTIX_INT inc)
{
  PASTIX_INT i;
  for (i=0;i<n;i++)
    new[i]=v[i*inc];
}

void end_vecteur(PASTIX_INT n,PASTIX_FLOAT *v,PASTIX_FLOAT new[SIZEMAX],PASTIX_INT inc)
{
  PASTIX_INT i;
  for (i=0;i<n;i++)
    v[i*inc]=new[i];
}

void init_matrice(char *t,PASTIX_INT n,PASTIX_INT m,PASTIX_FLOAT *a,PASTIX_FLOAT new[SIZEMAX][SIZEMAX],PASTIX_INT lda)
{
  PASTIX_INT i,j;
  for (i=0;i<n;i++)
    for (j=0;j<m;j++)
      if (*t=='N')
	new[j][i]=a[i+j*lda];
      else
	new[i][j]=a[i+j*lda];
}

void end_matrice(char *t,PASTIX_INT n,PASTIX_INT m,PASTIX_FLOAT *a,PASTIX_FLOAT new[SIZEMAX][SIZEMAX],PASTIX_INT lda)
{
  PASTIX_INT i,j;
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
