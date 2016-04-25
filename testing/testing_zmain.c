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
    char func[32];

    /* Check for number of arguments*/
    if ( argc < 2) {
        printf(" Proper Usage is : ./ztesting FUNC ...\n"
               "   - FUNC   : name of function to test\n");
        exit(1);
    }

    sscanf( argv[1], "%s", func   );

    argc -= 2;
    argv += 2;

    /*
     * LR operations
     */
    if ( strcmp(func, "RRADD") == 0 ) {
        info = testing_zrradd( argc, argv );
    } else {
        fprintf(stderr, "Function unknown\n");
    }

    if ( info == -1 ) {
        printf( "TESTING %s FAILED : incorrect number of arguments\n", func);
    } else if ( info == -2 ) {
        printf( "TESTING %s FAILED : not enough memory\n", func);
    } else if ( info == -3 ) {
        printf( "TESTING %s FAILED : invalid operation\n", func);
    }

    return EXIT_SUCCESS;
}
