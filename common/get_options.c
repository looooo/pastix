/**
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 */
#include "common.h"
#include <ctype.h>
#if defined(HAVE_GETOPT_H)
#include <getopt.h>
#endif  /* defined(HAVE_GETOPT_H) */
#include <string.h>

/**
 *******************************************************************************
 *
 * getfilename - Sets filename to source if source doesn't starts with '-'.
 * Otherwise, filename is set to defaultname.
 *
 *******************************************************************************
 *
 * @param[out] filename
 *          string to set to correct filename.
 *
 * @param[in] source
 *          possible source for filename.
 *
 * @param[in] defaultname
 *          default filename.
 *
 * Returns:
 *   0 if set to default.
 *   1 if set to source.
 */
static inline int
getfilename(char **filename, char *source, char *defaultname)
{
    if (source == NULL || source[0] == '-')
    {
        *filename = (char *) malloc((strlen(defaultname)+1)*sizeof(char));
        strcpy(*filename,defaultname);
        return 0;
    }
    *filename = (char *) malloc((strlen(source)+1)*sizeof(char));
    strcpy(*filename,source);
    return 1;
}

static inline void
pastix_usage(void)
{
    fprintf(stderr,
            "Matrix input (mandatory):\n"
            " -0 --rsa          : RSA/Matrix Market Fortran driver (only real)\n"
            " -1 --hb           : Harwell Boeing C driver\n"
            " -2 --ijv          : IJV coordinate C driver\n"
            " -3 --mm           : Matrix Market C driver\n"
            " -9 --lap          : Generate a Laplacian (5-points stencil)\n"
            " -x --xlap         : Generate an extended Laplacian (9-points stencil)\n"
            " -G --graph        : SCOTCH Graph file\n"
            "\n"
            "Architecture arguments:\n"
            " -t --threads      : Number of threads per node\n"
            " -g --gpus         : Number of gpus per node\n"
            " -s --sched        : Set the default scheduler (default: 1)\n"
            "                     0: Sequential, 1: Static, 2:PaRSEC\n"
            "\n"
            "Optional arguments:\n"
            " -f --fact                     : Choose factorization method (default: LU)\n"
            "                                 0: Cholesky, 1: LDL^[th], 2: LU, 3:LL^t, 4:LDL^t\n"
            "                                 3 and 4 are for complex matrices only\n"
            " -o --ord                      : Choose between ordering libraries [scotch|ptscotch|metis]\n"
            " -i --iparm <IPARM_ID> <value> : set any given integer parameter\n"
            " -d --dparm <DPARM_ID> <value> : set any given floating parameter\n"
            "\n"
            " -v --verbose[=lvl] : extra verbose output\n"
            " -h --help          : this message\n"
            "\n"
            );
}

#define GETOPT_STRING "0:1:2:3:9:x:G:t:g:d:f:o:i:d:v::h"

#if defined(HAVE_GETOPT_LONG)
static struct option long_options[] =
{
    {"0",           required_argument,  0, '0'},
    {"rsa",         required_argument,  0, '0'},
    {"1",           required_argument,  0, '1'},
    {"hb",          required_argument,  0, '1'},
    {"2",           required_argument,  0, '2'},
    {"ijv",         required_argument,  0, '2'},
    {"3",           required_argument,  0, '3'},
    {"mm",          required_argument,  0, '3'},
    {"9",           required_argument,  0, '9'},
    {"lap",         required_argument,  0, '9'},
    {"x",           required_argument,  0, 'x'},
    {"xlap",        required_argument,  0, 'x'},
    {"G",           required_argument,  0, 'G'},
    {"graph",       required_argument,  0, 'G'},

    {"threads",     required_argument,  0, 't'},
    {"t",           required_argument,  0, 't'},
    {"gpus",        required_argument,  0, 'g'},
    {"g",           required_argument,  0, 'g'},
    {"sched",       required_argument,  0, 's'},
    {"s",           required_argument,  0, 's'},

    {"ord",         required_argument,  0, 'o'},
    {"o",           required_argument,  0, 'o'},
    {"fact",        required_argument,  0, 'f'},
    {"f",           required_argument,  0, 'f'},
    {"iparm",       required_argument,  0, 'i'},
    {"i",           required_argument,  0, 'i'},
    {"dparm",       required_argument,  0, 'd'},
    {"d",           required_argument,  0, 'd'},

    {"verbose",     optional_argument,  0, 'v'},
    {"v",           optional_argument,  0, 'v'},
    {"help",        no_argument,        0, 'h'},
    {"h",           no_argument,        0, 'h'},
    {0, 0, 0, 0}
};
#endif  /* defined(HAVE_GETOPT_LONG) */

void
pastix_getOptions( int argc, char **argv,
                   pastix_int_t *iparam, double *dparam,
                   pastix_driver_t *driver, char **filename )
{
    int opt = 0;
    int c;
    (void)dparam;

    if (argc == 1) {
        pastix_usage(); exit(0);
    }

    do
    {
#if defined(HAVE_GETOPT_LONG)
        c = getopt_long_only(argc, argv, "",
                             long_options, &opt);
#else
        c = getopt(argc, argv, GETOPT_STRING);
        (void) opt;
#endif  /* defined(HAVE_GETOPT_LONG) */

        //       printf("%c: %s = %s\n", c, long_options[opt].name, optarg);
        switch(c)
        {
        case '0':
            *driver = PastixDriverRSA;
            getfilename( filename, optarg, "rsaname" );
            break;

        case '1':
            *driver = PastixDriverHB;
            getfilename( filename, optarg, "hbname" );
            break;

        case '2':
            *driver = PastixDriverIJV;
            getfilename( filename, optarg, "ijvname" );
            break;

        case '3':
            *driver = PastixDriverMM;
            getfilename( filename, optarg, "mmname" );
            break;

        case '9':
            *driver = PastixDriverLaplacian;
            getfilename( filename, optarg, "d:1000" );
            break;

        case 'x':
            *driver = PastixDriverXLaplacian;
            getfilename( filename, optarg, "d:1000" );
            break;

        case 'G':
            *driver = PastixDriverGraph;
            getfilename( filename, optarg, "graphname" );
            break;

        case 't': iparam[IPARM_THREAD_NBR] = atoi(optarg); break;
        case 'g': iparam[IPARM_GPU_NBR] = atoi(optarg); break;

        case 'o':
            if (strcmp(optarg, "scotch") == 0)
            {
                iparam[IPARM_ORDERING] = PastixOrderScotch;
            }
            else if (strcmp(optarg, "metis") == 0)
            {
                iparam[IPARM_ORDERING] = PastixOrderMetis;
            }
            else if (strcmp(optarg, "ptscotch") == 0)
            {
                iparam[IPARM_ORDERING] = PastixOrderPtScotch;
            }
            else if (strcmp(optarg, "personal") == 0)
            {
                iparam[IPARM_ORDERING] = PastixOrderPersonal;
            }
            else {
                fprintf(stderr, "Wrong values (ord=%s)!!!\nPossible values for ordering are: scotch, metis and ptscotch (Default scotch is chosen)\n", optarg);
            }
            break;

        case 'f': {
            int factotype = atoi( optarg );
            if ( (factotype >= 0) && (factotype <= 3)){
                iparam[IPARM_FACTORIZATION] = factotype;
            }
        }
            break;

        case 's': {
            int schedtype = atoi( optarg );
            if ( (schedtype >= 0) && (schedtype <= 3)){
                iparam[IPARM_SCHEDULER] = schedtype;
            }
        }
            break;

        case 'i':
        {
            int iparm_idx, iparm_val;

            /* Get iparm index */
            iparm_idx = iparm_to_int( optarg );
            if ( iparm_idx == -1 ) {
                fprintf(stderr, "Couldn't find value for %s\n", optarg );
                goto unknown_option;
            }

            /* Get iparm value */
            iparm_val = api_to_int( argv[optind] );
            if ( iparm_val == -1 ){
                fprintf(stderr, "Couldn't find value for %s\n", argv[optind] );
                goto unknown_option;
            }
            iparam[iparm_idx] = iparm_val;
        }
        break;

        case 'd':
        {
            int dparm_idx;
            double dparm_val;

            /* Get iparm index */
            dparm_idx = dparm_to_int( optarg );
            if ( dparm_idx == -1 )
                goto unknown_option;

            /* Get iparm value */
            dparm_val = atof( argv[optind] );
            dparam[dparm_idx] = dparm_val;
        }
        break;

        case 'v':
            if(optarg)  iparam[IPARM_VERBOSE] = atoi(optarg);
            else        iparam[IPARM_VERBOSE] = 2;
            break;

        case 'h': pastix_usage(); exit(0);

        case '?': /* getopt_long already printed an error message. */
            exit(1);
        default:
            break;
        }
    } while(-1 != c);

    return;

  unknown_option:
    fprintf(stderr, "ERROR: Unknown option\n");
    pastix_usage(); exit(0);
}
