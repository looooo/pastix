/**
 *
 * @file z_rradd_tests.c
 *
 * Tests and validate the Xrradd routine.
 *
 * @version 5.1.0
 * @author Gregoire Pichon
 * @date 2016-11-24
 *
 * @precisions normal z -> c d s
 *
 **/
#define _GNU_SOURCE
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "cblas.h"
#include "lapacke.h"
#include <common.h>

#define PRINT_RES(_ret_)                        \
    if(_ret_ == -1) {                           \
        printf("UNDEFINED\n");                  \
    }                                           \
    else if(_ret_ > 0) {                        \
        printf("FAILED(%d)\n", _ret_);          \
        err++;                                  \
    }                                           \
    else {                                      \
        printf("SUCCESS\n");                    \
    }

int
z_rradd_test(){
    return 0;
}

int main (int argc, char **argv)
{
    int err = 0;
    int ret;
    pastix_int_t i;

    for (i=0; i<10; i++){
        printf("   -- Test RRADD M=%ld ", i);

        ret = z_rradd_test();
        PRINT_RES(ret);
    }


    if( err == 0 ) {
        printf(" -- All tests PASSED --\n");
        return EXIT_SUCCESS;
    }
    else
    {
        printf(" -- %d tests FAILED --\n", err);
        return EXIT_FAILURE;
    }

}
