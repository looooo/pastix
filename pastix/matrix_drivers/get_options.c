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
#if defined(HAVE_GETOPT_H)
#include <getopt.h>
#endif  /* defined(HAVE_GETOPT_H) */
#include "drivers.h"

/**
 * Function: getfilename
 *
 * Sets filename to source if source doesn't starts with '-'.
 * Otherwise, filename is set to defaultname.
 *
 * Parameters:
 *   filename    - string to set to correct filename.
 *   source      - possible source for filename.
 *   defaultname - default filename.
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

/* /\* */
/*   Function: str_tolower */

/*   Rewrites *string* in lower case. */

/*   Parameters: */
/*     string - string to rewrite in lower case. */
/* *\/ */
/* int str_tolower(char * string) */
/* { */
/*   int j = 0; */
/*   while (string[j] != '\0') */
/*     { */
/*       string[j] = (char)tolower(string[j]); */
/*       j++; */
/*     } */
/*   return EXIT_SUCCESS; */
/* } */

/* /\* */
/*   Function: global_usage */

/*   Print usage corresponding to all pastix exemples. */

/*   Parameters: */
/*     mpi_comm - MPI communicator. */
/*     argv     - program argument */

/* *\/ */
/* void global_usage(MPI_Comm mpi_comm, char ** argv) */
/* { */
/*   int rank; */
/*   (void)mpi_comm; */

/*   MPI_Comm_rank(mpi_comm, &rank); */
/*   if (rank == 0) */
/*     { */
/*       fprintf(stdout, "Usage : %s [option] \n",argv[0]); */
/*       fprintf(stdout, "\toptions : \n"); */
/*       fprintf(stdout, "\t\t -rsa     [filename]          driver RSA (use Fortran, double only) \n"); */
/*       fprintf(stdout, "\t\t -hb      [filename]          driver Harwell Boeing (Similar RSA but in C with complex support\n"); */
/*       fprintf(stdout, "\t\t -ccc     [filename]          driver CCC\n"); */
/*       fprintf(stdout, "\t\t -rcc     [filename]          driver RCC\n"); */
/*       fprintf(stdout, "\t\t -olaf    [filename]          driver OLAF\n"); */
/*       fprintf(stdout, "\t\t -peer    [filename]          driver PEER\n"); */
/*       fprintf(stdout, "\t\t -petsc_s [filename]          driver PETSc symmetric\n"); */
/*       fprintf(stdout, "\t\t -petsc_h [filename]          driver PETSc hermitian\n"); */
/*       fprintf(stdout, "\t\t -petsc_u [filename]          driver PETSc unsymmetric\n"); */
/*       fprintf(stdout, "\t\t -3files  [filename]          driver IJV 3files \n"); */
/*       fprintf(stdout, "\t\t -mm      [filename]          driver Matrix Market\n"); */
/*       fprintf(stdout, "\t\t -dmm     [filename]          driver Matrix Market (distributed)\n"); */
/* #ifdef FDUPROS */
/*       fprintf(stdout, "\t\t -fdup    [filename]          driver from Fabrice Dupros\n"); */
/*       fprintf(stdout, "\t\t -fdupd   [filename]          driver from Fabrice Dupros, distributed\n"); */
/* #endif */
/*       fprintf(stdout, "\t\t -ord     <scotch|metis>      select ordering library\n"); */
/*       fprintf(stdout, "\t\t -lap     <integer>           generate a laplacian of size <integer>\n"); */
/*       fprintf(stdout, "\t\t -incomp  <integer> <integer> incomplete factorization, with the given level of fill [1-5],\n"); */
/*       fprintf(stdout, "\t\t                              and amalgamation [10-70]\n"); */
/*       fprintf(stdout, "\t\t -ooc     <integer>           Memory limit in Mo/percent depending on compilation options\n"); */
/*       fprintf(stdout, "\t\t -kass    <integer>           kass, with the given amalgamation\n"); */
/*       fprintf(stdout, "\t\t -t       <integer>           define thread number\n"); */
/*       fprintf(stdout, "\t\t -v       <integer>           define verbose level (1,2 or 3)\n"); */
/*       fprintf(stdout, "\t\t -iparm   <IPARM_ID> <value>  set an integer parameter\n"); */
/*       fprintf(stdout, "\t\t -dparm   <DPARM_ID> <value>  set a floating parameter\n"); */

/*       /\*       fprintf(stdout, "\t\t b         driver \"Fabrice Dupros\"\n"); *\/ */
/*       fprintf(stdout, "\t\t -h                          print this help\n"); */
/*     } */
/* } */
/* /\* */
/*   Function: get_options */

/*   Get options from argv. */

/*   Parameters: */
/*   argc          - number of arguments. */
/*   argv          - argument tabular. */
/*   driver_type   - type of driver (output, -1 if not set). */
/*   filename      - Matrix filename (output). */
/*   nbmatrices    - number of matrices in arguments. */
/*   nbthread      - number of thread (output, 1 if not set). */
/*   verbose       - verbose level 1,2 or 3 */
/*   ordering      - ordering to choose (see <API_ORDER>). */
/*   incomplete    - indicate if -incomp is present */
/*   level_of_fill - Level of fill for incomplete factorization. */
/*   amalgamation  - Amalgamation for kass. */
/*   ooc           - Out-of-core limite (Mo or percent depending on compilation option) */
/*   size          - Size of the matrix (generated matrix only) */
/* *\/ */
/* int get_options(int              argc, */
/*              char           **argv, */
/*              pastix_driver_t  **driver_type, */
/*              char          ***filename, */
/*              int             *nbmatrices, */
/*              int             *nbthread, */
/*              int             *verbose, */
/*              int             *ordering, */
/*              int             *incomplete, */
/*              int             *level_of_fill, */
/*              int             *amalgamation, */
/*              int             *ooc, */
/*              pastix_int_t    *size) */
/* { */

/*   int i = 1; */
/*   int maxmatrices = 10; */

/*   (*driver_type) = (pastix_driver_t*)malloc(maxmatrices*sizeof(pastix_driver_t)); */
/*   (*filename)    = (char **       )malloc(maxmatrices*sizeof(char*)); */
/*   *nbmatrices    = 0; */
/*   *nbthread      = 1; */
/*   *verbose       = 1; */
/*   *size          = 0; */
/*   *ordering      = API_ORDER_SCOTCH; */
/*   *incomplete    = API_NO; */
/*   *level_of_fill = 0; */
/*   *amalgamation  = 5; */
/*   *ooc           = 2000; */

/*   if (argc == 1) */
/*     goto usage; */
/*   while(i < argc) */
/*     { */
/*       if (argv[i][0] == '-') */
/*      { */

/*        switch (argv[i][1]) { */

/*        case 'c': */
/*        case 'C': */
/*          str_tolower(argv[i]); */
/*          if (strcmp(argv[i], "-chb") == 0) */
/*            { */
/*              (*driver_type)[(*nbmatrices)] = CHB; */
/*              i+=getfilename(&(*filename)[(*nbmatrices)], */
/*                             (i+1<argc)?argv[i+1]:NULL, "rsaname"); */
/*              (*nbmatrices)++; */
/*            } */
/*          else */
/*            { */
/*              if (strcmp(argv[i], "-ccc") == 0) */
/*                { */
/*                  (*driver_type)[(*nbmatrices)] = CCC; */
/*                  i+= getfilename(&(*filename)[(*nbmatrices)], */
/*                                  (i+1<argc)?argv[i+1]:NULL, "dirname"); */
/*                  (*nbmatrices)++; */
/*                } */
/*              else */
/*                { */
/*                  if (strcmp(argv[i], "-cscd") == 0) */
/*                    { */
/*                      (*driver_type)[(*nbmatrices)] = CSCD; */
/*                      i+= getfilename(&(*filename)[(*nbmatrices)], */
/*                                      (i+1<argc)?argv[i+1]:NULL, "dirname"); */
/*                      (*nbmatrices)++; */
/*                    } */
/*                  else */
/*                    goto unknown_option; */
/*                } */
/*            } */
/*          break; */
/*        case 'd': */
/*        case 'D': */
/*          str_tolower(argv[i]); */
/*          if (strcmp(argv[i], "-dparmfile") == 0) */
/*            { */
/*              i++; */
/*            } */
/*          else if (strcmp(argv[i], "-dparm") == 0) */
/*            { */
/*              i+=2; */
/*            } */
/*          else if (strcmp(argv[i],"-dmm") == 0 || */
/*                   strcmp(argv[i],"distributedmatrixmarket") == 0) */
/*            { */
/*              (*driver_type)[(*nbmatrices)] = MMD; */
/*              i += getfilename(&(*filename)[(*nbmatrices)], */
/*                               (i+1<argc)?argv[i+1]:NULL, "mmname"); */
/*              (*nbmatrices)++; */
/*            } */
/*          else */
/*            goto unknown_option; */
/*          break; */
/* #ifdef FDUPROS */
/*        case 'f': */
/*        case 'F': */
/*          { */
/*            str_tolower(argv[i]); */
/*            if (strcmp(argv[i], "-fdup") == 0) */
/*              { */
/*                (*driver_type)[(*nbmatrices)] = FDUP; */
/*                i+=getfilename(&(*filename)[(*nbmatrices)], */
/*                               (i+1<argc)?argv[i+1]:NULL, "dirname"); */
/*                (*nbmatrices)++; */
/*              } */
/*            else */
/*              { */
/*                if (strcmp(argv[i], "-fdupd") == 0) */
/*                  { */
/*                    (*driver_type)[(*nbmatrices)] = FDUP_DIST; */
/*                    i+=getfilename(&(*filename)[(*nbmatrices)], */
/*                                   (i+1<argc)?argv[i+1]:NULL, "dirname"); */
/*                    (*nbmatrices)++; */
/*                  } */
/*                else */
/*                  { */
/*                    goto unknown_option; */
/*                  } */
/*              } */
/*          } */
/*          break; */
/* #endif */

/*        case 'h': */
/*        case 'H': */
/*          str_tolower(argv[i]); */
/*          if (strcmp(argv[i],"-h") ==0 || strcmp(argv[i],"-help") ==0) */
/*            goto usage; */
/*          else */
/*            if (strcmp(argv[i],"-hb") == 0 || */
/*                strcmp(argv[i],"-harwell-boeing") == 0 || */
/*                strcmp(argv[i],"-harwellboeing") == 0) */
/*              { */
/*                (*driver_type)[(*nbmatrices)] = HB; */
/*                i+=getfilename(&(*filename)[(*nbmatrices)], */
/*                               (i+1<argc)?argv[i+1]:NULL, "rsaname"); */
/*                (*nbmatrices)++; */
/*              } */
/*          break; */
/*        case 'i': */
/*        case 'I': */
/*          str_tolower(argv[i]); */
/*          if (strcmp(argv[i], "-incomp") == 0) */
/*            { */
/*              *incomplete    = API_YES; */
/*              *level_of_fill = atoi(argv[i+1]); */
/*              i++; */
/*              *amalgamation  = atoi(argv[i+1]); */
/*              i++; */
/*            } */
/*          else if (strcmp(argv[i], "-iparmfile") == 0) */
/*            { */
/*              i++; */
/*            } */
/*          else if (strcmp(argv[i], "-iparm") == 0) */
/*            { */
/*              i+=2; */
/*            } */
/*          else */
/*            goto unknown_option; */
/*          break; */

/*        case 'k': */
/*        case 'K': */
/*          str_tolower(argv[i]); */
/*          if (strcmp(argv[i], "-kass") == 0) */
/*            { */
/*              *level_of_fill = -1; */
/*              *amalgamation  = atoi(argv[i+1]); */
/*              i++; */
/*            } */
/*          else */
/*            goto unknown_option; */
/*          break; */
/*        case 'l': */
/*        case 'L': */
/*          str_tolower(argv[i]); */
/*          if (strcmp(argv[i], "-lap") == 0) */
/*            { */
/*              (*driver_type)[(*nbmatrices)] = LAPLACIAN; */
/*              *size = atoi(argv[i+1]); */
/*              (*filename)[(*nbmatrices)] = NULL; */
/*              if (0 == *size) */
/*                goto unknown_option; */
/*              i++; */
/*              (*nbmatrices)++; */
/*            } */
/*          else */
/*            goto unknown_option; */
/*          break; */

/*        case 'm': */
/*        case 'M': */
/*          str_tolower(argv[i]); */
/*          if (strcmp(argv[i],"-mm") == 0 || */
/*              strcmp(argv[i],"matrixmarket") == 0) */
/*            { */
/*              (*driver_type)[(*nbmatrices)] = MM; */
/*              i += getfilename(&(*filename)[(*nbmatrices)], */
/*                               (i+1<argc)?argv[i+1]:NULL, "mmname"); */
/*              (*nbmatrices)++; */
/*            } */
/*          else */
/*            goto unknown_option; */
/*          break; */

/*        case 'o': */
/*        case 'O': */
/*          str_tolower(argv[i]); */
/*          if (strcmp(argv[i],"-olaf") == 0) */
/*            { */
/*              (*driver_type)[(*nbmatrices)] = OLAF; */
/*              i+= getfilename(&(*filename)[(*nbmatrices)], */
/*                              (i+1<argc)?argv[i+1]:NULL, "olafcsr"); */
/*              (*nbmatrices)++; */
/*            } */
/*          else */
/*            { */
/*              if (strcmp(argv[i],"-ord") == 0) */
/*                { */
/*                  if (EXIT_FAILURE == */
/*                      getordering(ordering,(i+1<argc)?argv[i+1]:NULL)) */
/*                    goto usage; */
/*                  else */
/*                    i++; */
/*                } */
/*              else */
/*                { */
/*                  if (strcmp(argv[i], "-ooc") == 0) */
/*                    { */
/*                      *ooc = atoi(argv[i+1]); */
/*                      if (0 == *ooc) */
/*                        goto unknown_option; */
/*                      i++; */
/*                    } */
/*                  else */
/*                    goto unknown_option; */
/*                } */
/*            } */
/*          break; */

/*        case 'p': */
/*        case 'P': */
/*          str_tolower(argv[i]); */
/*          if (strcmp(argv[i],"-peer") == 0) */
/*            { */
/*              (*driver_type)[(*nbmatrices)] = PEER; */
/*              i+= getfilename(&(*filename)[(*nbmatrices)], */
/*                              (i+1<argc)?argv[i+1]:NULL, "rsaname"); */
/*              (*nbmatrices)++; */
/*            } */
/*          else */
/*            { */
/*              if ( strcmp(argv[i], "-petsc_s") == 0 ) */
/*                { */
/*                  (*driver_type)[(*nbmatrices)] = PETSCS; */
/*                  i+= getfilename(&(*filename)[(*nbmatrices)], */
/*                                  (i+1<argc)?argv[i+1]:NULL, "PETSCFILE"); */
/*                  (*nbmatrices)++; */
/*                } */
/*              else */
/*                { */
/*                  if ( strcmp(argv[i], "-petsc_u") == 0 ) */
/*                    { */
/*                      (*driver_type)[(*nbmatrices)] = PETSCU; */
/*                      i+= getfilename(&(*filename)[(*nbmatrices)], */
/*                                      (i+1<argc)?argv[i+1]:NULL, "PETSCFILE"); */
/*                      (*nbmatrices)++; */
/*                    } */
/*                else */
/*                  { */
/*                    if ( strcmp(argv[i], "-petsc_h") == 0 ) */
/*                      { */
/*                        (*driver_type)[(*nbmatrices)] = PETSCH; */
/*                        i+= getfilename(&(*filename)[(*nbmatrices)], */
/*                                        (i+1<argc)?argv[i+1]:NULL, "PETSCFILE"); */
/*                        (*nbmatrices)++; */
/*                      } */
/*                    else { */
/*                      goto unknown_option; */
/*                    } */
/*                  } */
/*                } */
/*            } */
/*          break; */

/*        case 'r': */
/*        case 'R': */
/*          str_tolower(argv[i]); */
/*          if (strcmp(argv[i],"-rsa") == 0) */
/*            { */
/*              (*driver_type)[(*nbmatrices)] = RSA; */
/*              i+= getfilename(&(*filename)[(*nbmatrices)], */
/*                              (i+1<argc)?argv[i+1]:NULL, "rsaname"); */
/*              (*nbmatrices)++; */
/*            } */
/*          else */
/*            { */
/*              if (strcmp(argv[i],"-rcc") == 0) */
/*                { */
/*                  (*driver_type)[(*nbmatrices)] = RCC; */
/*                  i+= getfilename(&(*filename)[(*nbmatrices)], */
/*                                  (i+1<argc)?argv[i+1]:NULL, "dirname"); */
/*                  (*nbmatrices)++; */
/*                } */
/*              else */
/*                goto unknown_option; */
/*            } */
/*          break; */

/*        case 't': */
/*        case 'T': */
/*          str_tolower(argv[i]); */
/*          if (strcmp(argv[i], "-t") == 0) */
/*            { */
/*              *nbthread = atoi(argv[i+1]); */
/*              if (0 == *nbthread) */
/*                goto unknown_option; */
/*              i++; */
/*            } */
/*          else */
/*            goto unknown_option; */
/*          break; */

/*        case 'v': */
/*          str_tolower(argv[i]); */
/*          if (strcmp(argv[i], "-v") == 0) */
/*            { */
/*              *verbose = atoi(argv[i+1]); */
/*              if (0 == *verbose) */
/*                goto unknown_option; */
/*              i++; */
/*              (*verbose)--; */
/*            } */
/*          else */
/*            goto unknown_option; */
/*          break; */
/*        case '3': */
/*          str_tolower(argv[i]); */
/*          if (strcmp(argv[i],"-3files") == 0) */
/*            { */
/*              (*driver_type)[(*nbmatrices)] = THREEFILES; */
/*              i+= getfilename(&(*filename)[(*nbmatrices)], */
/*                              (i+1<argc)?argv[i+1]:NULL, "dirname"); */
/*              (*nbmatrices)++; */
/*            } */
/*          else */
/*            goto unknown_option; */
/*          break; */

/*        default: */
/*        unknown_option: */
/*          fprintf(stderr, */
/*                  "ERROR: main: unprocessed option (\"%s\")\n", argv[i]); */
/*        usage: */
/*          global_usage(MPI_COMM_WORLD,argv); */
/*          free(*driver_type); */
/*          free(*filename); */
/*          MPI_Finalize(); */
/*          return EXIT_FAILURE; */

/*        } */
/*      } */
/*       if (maxmatrices == (*nbmatrices)) */
/*      { */
/*        maxmatrices *=2; */
/*        (*driver_type) = (pastix_driver_t*)realloc((*driver_type), */
/*                                                 maxmatrices * */
/*                                                 sizeof(pastix_driver_t)); */
/*        (*filename)    = (char **       )realloc((*filename), */
/*                                                 maxmatrices*sizeof(char*)); */
/*      } */
/*       i++; */
/*     } */

/*   /\\* default driver *\\/ */
/*   if ((*nbmatrices) == 0) */
/*     { */
/*       (*driver_type)[0] = RSA; */
/*       (*filename)[0] = malloc((strlen("rsaname")+1)*sizeof(char)); */
/*       strcpy((*filename)[0],"rsaname"); */
/*       (*nbmatrices) ++; */
/*     } */
/*   return EXIT_SUCCESS; */
/* } */

/* /\\* */
/*   Function: d_get_idparm */

/*   Get options from argv. */

/*   Parameters: */
/*   argc          - number of arguments. */
/*   argv          - argument tabular. */
/*   iparm         - type of driver (output, -1 if not set). */
/*   dparm         - type of driver (output, -1 if not set). */
/* *\\/ */
/* int d_get_idparm(int            argc, */
/*                  char         **argv, */
/*                  pastix_int_t  *iparm, */
/*                  double        *dparm) */
/* { */
/*   int i             = 1; */

/*   while(i < argc) */
/*     { */
/*       if (argv[i][0] == '-') */
/*      { */
/*        switch (argv[i][1]) { */

/*        case 'd': */
/*        case 'D': */
/*          str_tolower(argv[i]); */
/*          if (strcmp(argv[i], "-dparmfile") == 0) */
/*            { */
/*              i++; */
/*              api_dparmreader(argv[i], dparm); */
/*            } */
/*          else if (strcmp(argv[i], "-dparm") == 0) */
/*            { */
/*              int    dparm_idx; */
/*              double value; */
/*              char * endptr; */
/*              i++; */
/*              dparm_idx = (int)strtol(argv[i], &endptr, 10); */
/*              if (endptr == argv[i]) */
/*                { */
/*                  if( 1 == api_str_to_int(argv[i], &dparm_idx)) */
/*                    goto unknown_option; */
/*                } */
/*              i++; */
/*              value = (double)strtod(argv[i], &endptr); */
/*              if (endptr == argv[i]) */
/*                goto unknown_option; */
/*              dparm[dparm_idx] = value; */

/*            } */
/*          break; */
/*        case 'i': */
/*        case 'I': */
/*          str_tolower(argv[i]); */
/*          if (strcmp(argv[i], "-iparmfile") == 0) */
/*            { */
/*              i++; */
/*              api_iparmreader(argv[i], iparm); */
/*            } */
/*          else if (strcmp(argv[i], "-iparm") == 0) */
/*            { */
/*              int iparm_idx, value; */
/*              char * endptr; */
/*              i++; */
/*              iparm_idx = (int)strtol(argv[i], &endptr, 10); */
/*              if (endptr == argv[i]) */
/*                { */
/*                  if( 1 == api_str_to_int(argv[i], &iparm_idx)) */
/*                    goto unknown_option; */
/*                } */
/*              i++; */
/*              value = (int)strtol(argv[i], &endptr, 10); */
/*              if (endptr == argv[i]) */
/*                if( 1 == api_str_to_int(argv[i], &value)) */
/*                  goto unknown_option; */
/*              iparm[iparm_idx] = value; */
/*            } */
/*          break; */
/*        } */
/*      } */
/*       i++; */
/*     } */

/*   return EXIT_SUCCESS; */
/*  unknown_option: */
/*   fprintf(stderr, */
/*        "ERROR: main: unprocessed option (\"%s\")\n", argv[i]); */
/*   global_usage(MPI_COMM_WORLD,argv); */
/*   MPI_Finalize(); */
/*   return EXIT_FAILURE; */

/* } */


/* /\\* */
/*   Function: s_get_idparm */

/*   Get options from argv. */

/*   Parameters: */
/*   argc          - number of arguments. */
/*   argv          - argument tabular. */
/*   iparm         - type of driver (output, -1 if not set). */
/*   dparm         - type of driver (output, -1 if not set). */
/* *\\/ */
/* int s_get_idparm(int            argc, */
/*                  char         **argv, */
/*                  pastix_int_t  *iparm, */
/*                  float        *dparm) { */
/*     int i, ret=EXIT_SUCCESS; */
/*     double *mydparm = malloc(DPARM_SIZE*sizeof(double)); */
/*     for (i = 0; i < DPARM_SIZE; i++) */
/*         mydparm[i] = (double)(dparm[i]); */
/*     ret = d_get_idparm(argc, argv, iparm, mydparm); */
/*     if (ret != EXIT_SUCCESS) return ret; */
/*     for (i = 0; i < DPARM_SIZE; i++) */
/*         dparm[i] = (float)(mydparm[i]); */
/*     return ret; */
/* } */



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

#define GETOPT_STRING "0:1:2:3:4:5:6:7:8:9:G:t:g:o:i:d:v::h"

#if defined(HAVE_GETOPT_LONG)
static struct option long_options[] =
{
    {"0",           required_argument,  0, '0'},
    {"rsa",         required_argument,  0, '0'},
    {"1",           required_argument,  0, '1'},
    {"hb",          required_argument,  0, '1'},
    {"2",           required_argument,  0, '2'},
    {"ccc",         required_argument,  0, '2'},
    {"3",           required_argument,  0, '3'},
    {"rcc",         required_argument,  0, '3'},
    {"4",           required_argument,  0, '4'},
    {"olaf",        required_argument,  0, '4'},
    {"5",           required_argument,  0, '5'},
    {"peer",        required_argument,  0, '5'},
    {"7",           required_argument,  0, '7'},
    {"ijv",         required_argument,  0, '7'},
    {"8",           required_argument,  0, '8'},
    {"mm",          required_argument,  0, '8'},
    {"9",           required_argument,  0, '9'},
    {"lap",         required_argument,  0, '9'},
    {"dmm",         required_argument,  0, 'A'},
    {"fdup",        required_argument,  0, 'B'},
    {"fdupd",       required_argument,  0, 'C'},
    {"petsc_s",     required_argument,  0, 'D'},
    {"petsc_h",     required_argument,  0, 'E'},
    {"petsc_u",     required_argument,  0, 'F'},
    {"G",           required_argument,  0, 'G'},
    {"graph",       required_argument,  0, 'G'},

    {"threads",     required_argument,  0, 't'},
    {"t",           required_argument,  0, 't'},
    {"gpus",        required_argument,  0, 'g'},
    {"g",           required_argument,  0, 'g'},

    {"ord",         required_argument,  0, 'o'},
    {"o",           required_argument,  0, 'o'},
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

void pastix_ex_getoptions(int argc, char **argv,
                          pastix_int_t *iparam, double *dparam,
                          pastix_driver_t *driver, char **filename )
{
    int opt = 0;
    int c;

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
            *driver = PastixDriverCCC;
            getfilename( filename, optarg, "cccname" );
            break;

        case '3':
            *driver = PastixDriverRCC;
            getfilename( filename, optarg, "rccname" );
            break;

        case '4':
            *driver = PastixDriverOlaf;
            getfilename( filename, optarg, "olafname" );
            break;

        case '5':
            *driver = PastixDriverPeer;
            getfilename( filename, optarg, "peername" );
            break;

        case '7':
            *driver = PastixDriverIJV;
            getfilename( filename, optarg, "ijvname" );
            break;

        case '8':
            *driver = PastixDriverMM;
            getfilename( filename, optarg, "mmname" );
            break;

        case '9':
            *driver = PastixDriverLaplacian;
            getfilename( filename, optarg, "d:1000" );
            break;

        case 'A':
            *driver = PastixDriverDMM;
            getfilename( filename, optarg, "dmmname" );
            break;

        case 'B':
            *driver = PastixDriverBRGM;
            getfilename( filename, optarg, "brgmname" );
            break;

        case 'C':
            *driver = PastixDriverBRGMD;
            getfilename( filename, optarg, "brgmdname" );
            break;

        case 'D':
            *driver = PastixDriverPetscS;
            getfilename( filename, optarg, "petscname" );
            break;

        case 'E':
            *driver = PastixDriverPetscH;
            getfilename( filename, optarg, "petscname" );
            break;

        case 'F':
            *driver = PastixDriverPetscU;
            getfilename( filename, optarg, "petscname" );
            break;

        case 'G':
            *driver = PastixDriverGraph;
            getfilename( filename, optarg, "graphname" );
            break;

        case 't': iparam[IPARM_THREAD_NBR] = atoi(optarg); break;
        case 'g': iparam[IPARM_GPUS_NBR] = atoi(optarg); break;

        case 'o':
            if (strcmp(optarg, "scotch") == 0)
            {
                iparam[IPARM_ORDERING] = API_ORDER_SCOTCH;
            }
            else if (strcmp(optarg, "metis") == 0)
            {
                iparam[IPARM_ORDERING] = API_ORDER_METIS;
            }
            else if (strcmp(optarg, "ptscotch") == 0)
            {
                iparam[IPARM_ORDERING] = API_ORDER_PTSCOTCH;
            }
            else {
                fprintf(stderr, "Wrong values (ord=%s)!!!\nPossible values for ordering are: scotch, metis and ptscotch (Default scotch is chosen)\n", optarg);
            }
            break;

        case 'i':
            break;
        case 'd':
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

    //    int verbose = iparam[IPARM_VERBOSE];
}
