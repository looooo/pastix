/**
 *
 * @file testing_zmain.h
 *
 * @precisions normal z -> c d s
 *
 **/
#ifndef TESTING_ZMAIN_H
#define TESTING_ZMAIN_H

#include "common.h"
#include <cblas.h>
#include <lapacke.h>
#include "blend/solver.h"
#include "pastix_zcores.h"

#define USAGE(name, args, details)                                      \
    printf(" Proper Usage is : ./testing_zmain " name " " args " with\n" \
           "   - " name "  : name of function to test\n"                \
           details); \

extern int IONE;
extern int ISEED[4];

int testing_zrradd(int argc, char **argv);
int testing_zgradd(int argc, char **argv);
int testing_zlrm2(int argc, char **argv);
int testing_zlrmm(int argc, char **argv);
int testing_zlrmge(int argc, char **argv);
int testing_zge2lr(int argc, char **argv);

#endif /* TESTING_ZMAIN_H */
