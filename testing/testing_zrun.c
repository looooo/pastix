/**
 *
 * @file testing_zmain.c
 *
 * @precisions normal z -> c d s
 *
 **/

#include <pastix.h>
#include <testing_zmain.h>

int   IONE     = 1;
int   ISEED[4] = {0,0,0,1};   /* initial seed for zlarnv() */

int main (int argc, char **argv)
{
    int info = -1;

    /* Check for number of arguments*/
    (void) argv;
    if ( argc != 1 ) {
        printf(" Proper Usage is : ./testing_zun without parameters\n ");
        exit(1);
    }

    int nb_params = 9;
    char **params;

    MALLOC_INTERN(params, nb_params, char*);
    int i;
    for (i=0; i<nb_params; i++){
        MALLOC_INTERN(params[i], 32, char);
    }


    /* TESTING RRADD */
    {
        double tolerance;
        pastix_int_t MA;
        pastix_int_t NA;
        int RA;
        pastix_int_t MB;
        pastix_int_t NB;
        int RB;
        pastix_int_t offx;
        pastix_int_t offy;

        int nb_params = 9;

        int tol = 0;
        for (tol = 1; tol < 100000000; tol*=10){
            tolerance = 1. / tol;
            MA = NA = MB = NB = 100;
            RA = RB = 10;
            offx = offy = 0;

            sprintf(params[0], "%f", tolerance);
            sprintf(params[1], "%ld", MA);
            sprintf(params[2], "%ld", NA);
            sprintf(params[3], "%d", RA);
            sprintf(params[4], "%ld", MB);
            sprintf(params[5], "%ld", NB);
            sprintf(params[6], "%d", RB);
            sprintf(params[7], "%ld", offx);
            sprintf(params[8], "%ld", offy);

            info = testing_zrradd(nb_params, params);
            if (info != 1){
                printf("ERROR in testing_zrradd with parameters: \n");
                printf("Tolerance = %f\n", tolerance);
                printf("MA = %ld, NA = %ld, RA = %d\n", MA, NA, RA);
                printf("MB = %ld, NB = %ld, RB = %d\n", MB, NB, RB);
                printf("Offx = %ld Offf = %ld\n", offx, offy);
                exit(1);
            }
        }
    }


    /* TESTING GRADD */

    /* TESTING LRMM */

    /* TESTING GE2LR */

    for (i=0; i<nb_params; i++){
        memFree_null(params[i]);
    }
    memFree_null(params);


    printf("\n\n ALL TESTS WERE SUCCESSFULLY PERFORMED\n\n");
    return EXIT_SUCCESS;
}
