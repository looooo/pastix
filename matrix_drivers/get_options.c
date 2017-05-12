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
#include "drivers.h"

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


/* int api_iparmreader(char * filename, pastix_int_t *iparmtab); */
/* int api_dparmreader(char * filename, double       *dparmtab); */

/*   Function: str_tolower */

/*   Rewrites *string* in lower case. */

/*   Parameters: */
/*   string - string to rewrite in lower case. */
static inline void
str_tolower(char * string)
{
    int j = 0;
    while (string[j] != '\0')
    {
        string[j] = (char)tolower(string[j]);
        j++;
    }
    return;
}

static inline void
pastix_ex_usage(void)
{
    fprintf(stderr,
            "Matrix input (mandatory):\n"
            " -0 --rsa          : RSA format (use Fortran, only real)\n"
            " -1 --hb           : Harwell Boeing (RSA Driver in C, support real/complex)\n"
            " -2 --ccc          : CCC format\n"
            " -3 --rcc          : RCC format\n"
            " -4 --olaf         : OLAF format\n"
            " -5 --peer         : PEER format\n"
            " -7 --ijv          : IJV 3 files format\n"
            " -8 --mm           : Matrix Market format\n"
            "    --dmm          : Matrix Market distributed format\n"
            " -9 --lap          : Generate a random 2D Laplacian of specified size\n"
            "    --fdup         : BRGM (Fabrice Dupros)\n"
            "    --fdupd        : BRGM (Fabrice Dupros) distributed\n"
            "    --petsc_s      : PETSc symmetric\n"
            "    --petsc_h      : PETSc hermitian\n"
            "    --petsc_u      : PETSc unsymmetric\n"
            " -G --graph        : SCOTCH Graph file\n"
            "\n"
            "Architecture arguments:\n"
            " -t --threads      : Number of threads\n"
            " -g --gpus         : Number of gpus\n"
            "\n"
            "Optional arguments:\n"
            " -o --ord                      : Choose between ordering libraries [scotch|ptscotch|metis]\n"
            " -i --iparm <IPARM_ID> <value> : set an integer parameter\n"
            " -d --dparm <DPARM_ID> <value> : set a floating parameter\n"
            "\n"
            " -v --verbose      : extra verbose output\n"
            " -h --help         : this message\n"
            "\n"
            );
}

#define GETOPT_STRING "0:1:2:3:9:x:G:t:g:o:i:d:f:s:v::h"

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

    {"ord",         required_argument,  0, 'o'},
    {"o",           required_argument,  0, 'o'},
    {"fact",        required_argument,  0, 'f'},
    {"f",           required_argument,  0, 'f'},
    {"sched",       required_argument,  0, 's'},
    {"s",           required_argument,  0, 's'},
    {"iparm",       required_argument,  0, 'i'},
    {"i",           required_argument,  0, 'i'},
    {"dparm",       required_argument,  0, 'd'},
    {"d",           required_argument,  0, 'd'},

    {"verbose",     optional_argument,  0, 'v'},
    {"v",           optional_argument,  0, 'v'},
    {"help",        no_argument,        0, 'h'},
    {"h",           no_argument,        0, 'h'},
    /*{"iparmfile",   no_argument,        0, 'i'},*/
    {"iparm",       no_argument,        0, 'i'},
    {"i",           no_argument,        0, 'i'},
    {0, 0, 0, 0}
};
#endif  /* defined(HAVE_GETOPT_LONG) */

void pastix_ex_getoptions(int argc, char **argv,
                          pastix_int_t *iparam, double *dparam,
                          pastix_driver_t *driver, char **filename )
{
    int opt = 0;
    int c;
    (void)dparam;

    if (argc == 1) {
        pastix_ex_usage(); exit(0);
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
                iparam[IPARM_ORDERING] = PastixOrderPtscotch;
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

        case 'h': pastix_ex_usage(); exit(0);

        case '?': /* getopt_long already printed an error message. */
            exit(1);
        default:
            break;
        }
    } while(-1 != c);

    return;

  unknown_option:
    fprintf(stderr, "ERROR: Unknown option\n");
    pastix_ex_usage(); exit(0);
}
