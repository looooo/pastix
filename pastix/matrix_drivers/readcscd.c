/**
 * @file readcscd.c
 *
 *  $COPYRIGHTS$
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 **/
#include "common.h"
#include "drivers.h"

/**
 * ******************************************************************************
 *
 * @ingroup pastix_csc_driver
 *
 * readCSCD - Read a real matrix in CSCD format.
 *
 *******************************************************************************
 *
 * @param[in] dirname
 *          Path to the directory containing matrix.
 *
 * @param[in] csc
 *          At exit, contains the matrix in csc format.
 *
 * @param[out] Rhs
 *          At exit, contains the right hand side.
 *
 * @param[in] pastix_comm
 *          MPI communicator used in PaStiX.
 *
 *******************************************************************************/
void
readCSCD( const char          *dirname,
          pastix_csc_t        *csc,
          void               **rhs,
          MPI_Comm             pastix_comm )
{
  const int          nbreltperline = 4; /* nbr of elt per line */
  FILE              *infile;
  char               line[BUFSIZ], file[BUFSIZ];
  int                myrank, nbproc, tmpint;
  long               tempint1,   tempint2,   tempint3,   tempint4;
  double             tempfloat1, tempfloat2, tempfloat3, tempfloat4;
  double            *rhs_temp    = NULL;
  double            *values      = NULL;
  int                vertloc, edgeloc;
  int                iterelt;
  int               *vectsize    = NULL;
  int               *vectsizercv = NULL;
  int                offset      = 1;
  char              *filename;
  int                i;
  (void)pastix_comm;


  MPI_Comm_rank( pastix_comm, &myrank);
  MPI_Comm_size( pastix_comm, &nbproc);

  filename = (char*)malloc(sizeof(char)*(strlen(dirname)+40));
  sprintf(filename,"%s/main",dirname);

  if (myrank == 0)
    {
      infile = fopen(filename, "r");
      if (infile==NULL)
  {
    fprintf(stderr,"cannot load %s\n", filename);
    exit(EXIT_FAILURE);
  }
      fgets(line, BUFSIZ, infile);
      sscanf(line, "%d", &tmpint); /* Read number of filename */
      fprintf(stdout, "Nombre de fichier %d\n", tmpint);
      fclose(infile);

      if (nbproc != tmpint)
  {
    if (myrank == 0)
      fprintf(stderr, "Veuillez fournir un communicateur MPI de %d processus\nActuellement, le communicateur contient %d processus\n", tmpint, nbproc);
    exit(EXIT_FAILURE);
  }
    }

  infile = fopen(filename, "r");
  if (infile == NULL)
    {
      fprintf(stderr,"cannot load %s\n", filename);
      exit(EXIT_FAILURE);
    }
  fgets(line, BUFSIZ, infile);

  for (i=0; i<=myrank; i++)
    {
      fgets(line, BUFSIZ, infile);
      sscanf(line, "%s", file);
    }
  fclose(infile);

  sprintf(filename,"%s/%s",dirname,file);
  infile = fopen(filename, "r");
  if (infile==NULL)
    {
      fprintf(stderr,"[P%d] cannot load %s\n", myrank, filename);
      exit(EXIT_FAILURE);
    }

  fgets(line, BUFSIZ, infile);
  sscanf(line, "%ld %ld", &tempint1,&tempint2);
  vertloc = tempint2;
/*   fprintf(stderr, "[P%d] rowlocal %ld\n",myrank, (long) vertloc); */
  fgets(line, BUFSIZ, infile);
  sscanf(line, "%ld", &tempint1);
  edgeloc = tempint1;
/*   fprintf(stderr, "[P%d] nzlocal %ld\n", myrank, (long) edgeloc); */

  csc->colptr      = (int *)   malloc((vertloc+1) * sizeof(int));
  if (   (csc->colptr)   == NULL)
    fprintf(stderr, "[P%d] z_cscdRead : Not enough memory for csc->colptr\n",myrank);
  csc->loc2glob = (int *)   malloc(  vertloc   * sizeof(int));
  if ( (csc->loc2glob) == NULL)
    fprintf(stderr, "[P%d] z_cscdRead : Not enough memory for csc->loc2glob\n",myrank);
  csc->rows      = (int *)   malloc(  edgeloc   * sizeof(int));
  if (  (csc->rows)   == NULL)
    fprintf(stderr, "[P%d] z_cscdRead : Not enough memory for csc->rows\n",myrank);
  csc->avals   = (double *) malloc(  edgeloc   * sizeof(double));
  values       = (double *) malloc(  edgeloc   * sizeof(double));
  if (  (csc->avals) == NULL || values == NULL )
    fprintf(stderr, "[P%d] z_cscdRead : Not enough memory for csc->avals\n",myrank);
  memset(values, 0, (edgeloc)*sizeof(double));
  *rhs     = (double *) malloc(  vertloc   * sizeof(double));
  rhs_temp = (double *) malloc(  vertloc   * sizeof(double));
  if (   (*rhs)  == NULL || rhs_temp  == NULL)
    fprintf(stderr, "[P%d] z_cscdRead : Not enough memory for *rhs\n",myrank);
  memset(*rhs, 0, (vertloc)*sizeof(double));

  /* Recuperation de Loc2glb*/
  for (iterelt=0; iterelt<vertloc+1-nbreltperline;iterelt+=nbreltperline )
    {
      fgets(line,BUFSIZ,infile);
      sscanf(line,"%ld %ld %ld %ld",
       &tempint1, &tempint2, &tempint3, &tempint4);
      (csc->loc2glob)[iterelt]   = (int)tempint1;
      (csc->loc2glob)[iterelt+1] = (int)tempint2;
      (csc->loc2glob)[iterelt+2] = (int)tempint3;
      (csc->loc2glob)[iterelt+3] = (int)tempint4;
    }
  switch (vertloc-iterelt)
    {
    case 1:
      fgets(line,BUFSIZ,infile);
      sscanf(line,"%ld",&tempint1);
      (csc->loc2glob)[iterelt] += (int)tempint1;
      iterelt++;
      break;
    case 2:
      fgets(line,BUFSIZ,infile);
      sscanf(line,"%ld %ld", &tempint1, &tempint2);
      (csc->loc2glob)[iterelt]   = (int)tempint1;
      (csc->loc2glob)[iterelt+1] = (int)tempint2;
      iterelt+=2;
      break;
    case 3:
      fgets(line,BUFSIZ,infile);
      sscanf(line,"%ld %ld %ld", &tempint1, &tempint2, &tempint3);
      (csc->loc2glob)[iterelt]   = (int)tempint1;
      (csc->loc2glob)[iterelt+1] = (int)tempint2;
      (csc->loc2glob)[iterelt+2] = (int)tempint3;
      iterelt+=3;
      break;
    default:
      break;
    }

  /* Recuperation de Colptr dans un style tres particulier... */
  for (iterelt=0; iterelt<vertloc+1+1-nbreltperline;iterelt+=nbreltperline )
    {
      fgets(line,BUFSIZ,infile);
      if (4 != sscanf(line,"%ld %ld %ld %ld", &tempint1, &tempint2, &tempint3, &tempint4))
  {
    fprintf(stderr, "ERROR: reading colptr\n");
    exit(1);
  }
      (csc->colptr)[iterelt]   = (int)tempint1;
      (csc->colptr)[iterelt+1] = (int)tempint2;
      (csc->colptr)[iterelt+2] = (int)tempint3;
      (csc->colptr)[iterelt+3] = (int)tempint4;
    }

  switch (vertloc-iterelt+1)
    {
    case 1:
      fgets(line,BUFSIZ,infile);
      if (1 != sscanf(line,"%ld",&tempint1))
  {
    fprintf(stderr, "ERROR: reading colptr\n");
    exit(1);
  }
      (csc->colptr)[iterelt] += (int)tempint1;
      iterelt++;
      break;
    case 2:
      fgets(line,BUFSIZ,infile);
      if (2 != sscanf(line,"%ld %ld", &tempint1, &tempint2))
  {
    fprintf(stderr, "ERROR: reading colptr\n");
    exit(1);
  }
      (csc->colptr)[iterelt]   = (int)tempint1;
      (csc->colptr)[iterelt+1] = (int)tempint2;
      iterelt+=2;
      break;
    case 3:
      fgets(line,BUFSIZ,infile);
      if (3 != sscanf(line,"%ld %ld %ld", &tempint1, &tempint2, &tempint3))
  {
    fprintf(stderr, "ERROR: reading colptr\n");
    exit(1);
  }
      (csc->colptr)[iterelt]   = (int)tempint1;
      (csc->colptr)[iterelt+1] = (int)tempint2;
      (csc->colptr)[iterelt+2] = (int)tempint3;
      iterelt+=3;
      break;
    default:
      break;
    }
  fprintf(stdout, "iterelt %ld, vertloc %ld\n", (long)iterelt, (long)vertloc);
  vectsize    = malloc(sizeof(int)*nbproc);
  vectsize    = memset(vectsize, 0, sizeof(int)*nbproc);
  vectsizercv = malloc(sizeof(int)*nbproc);
  vectsizercv = memset(vectsizercv, 0, sizeof(int)*nbproc);

  vectsize[myrank] = vertloc;

  if (vectsize    == NULL) fprintf(stderr, "[P%d] Erreur : alloc vectsize\n",    myrank);
  if (vectsizercv == NULL) fprintf(stderr, "[P%d] Erreur : alloc vectsizercv\n", myrank);

  MPI_Allreduce(vectsize, vectsizercv, nbproc, PASTIX_MPI_INT, MPI_SUM, pastix_comm);
  for (i=0; i<myrank; i++)
    offset += vectsizercv[i];

  /* Recuperation de ROW*/
  for (iterelt=0; iterelt<edgeloc+1-nbreltperline; iterelt+=4)
    {
      fgets(line,BUFSIZ,infile);
      sscanf(line,"%ld %ld %ld %ld", &tempint1,&tempint2,&tempint3,&tempint4);
      (csc->rows)[iterelt]   = (int)tempint1;
      (csc->rows)[iterelt+1] = (int)tempint2;
      (csc->rows)[iterelt+2] = (int)tempint3;
      (csc->rows)[iterelt+3] = (int)tempint4;

    }
  switch (edgeloc-iterelt)
    {
    case 1:
      fgets(line,BUFSIZ,infile);
      sscanf(line,"%ld",&tempint1);
      (csc->rows)[iterelt] = (int)tempint1;
      iterelt++;
      break;
    case 2:
      fgets(line,BUFSIZ,infile);
      sscanf(line,"%ld %ld", &tempint1, &tempint2);
      (csc->rows)[iterelt]   = (int)tempint1;
      (csc->rows)[iterelt+1] = (int)tempint2;
      iterelt+=2;
      break;
    case 3:
      fgets(line,BUFSIZ,infile);
      sscanf(line,"%ld %ld %ld", &tempint1, &tempint2, &tempint3);
      (csc->rows)[iterelt]   = (int)tempint1;
      (csc->rows)[iterelt+1] = (int)tempint2;
      (csc->rows)[iterelt+2] = (int)tempint3;
      iterelt+=3;
      break;
    }

  /* Recuperation de Aval */
  for (iterelt=0; iterelt<edgeloc+1-nbreltperline; iterelt+=4)
    {
      fgets(line,BUFSIZ,infile);
      sscanf(line,"%lf %lf %lf %lf",&tempfloat1,&tempfloat2,&tempfloat3,&tempfloat4);
      (values)[iterelt]   = (double)tempfloat1;
      (values)[iterelt+1] = (double)tempfloat2;
      (values)[iterelt+2] = (double)tempfloat3;
      (values)[iterelt+3] = (double)tempfloat4;
    }

  switch (edgeloc-iterelt)
    {
    case 1:
      fgets(line,BUFSIZ,infile);
      sscanf(line,"%lf",&tempfloat1);
      (values)[iterelt] = (double)tempfloat1;
      iterelt++;
      break;
    case 2:
      fgets(line,BUFSIZ,infile);
      sscanf(line,"%lf %lf",&tempfloat1,&tempfloat2);
      (values)[iterelt]   = (double)tempfloat1;
      (values)[iterelt+1] = (double)tempfloat2;
      iterelt+=2;
      break;
    case 3:
      fgets(line,BUFSIZ,infile);
      sscanf(line,"%lf %lf %lf",&tempfloat1,&tempfloat2,&tempfloat3);
      (values)[iterelt]   = (double)tempfloat1;
      (values)[iterelt+1] = (double)tempfloat2;
      (values)[iterelt+2] = (double)tempfloat3;
      iterelt+=3;
      break;
    }

  /* Recuperation de Rhs, second membre... */
  for (iterelt=0; iterelt<vertloc+1-nbreltperline; iterelt+=4)
    {
      fgets(line,BUFSIZ,infile);
      sscanf(line,"%lf %lf %lf %lf",&tempfloat1,&tempfloat2,&tempfloat3,&tempfloat4);
      (rhs_temp)[iterelt]   = (double)tempfloat1;
      (rhs_temp)[iterelt+1] = (double)tempfloat2;
      (rhs_temp)[iterelt+2] = (double)tempfloat3;
      (rhs_temp)[iterelt+3] = (double)tempfloat4;
    }

  switch (vertloc-iterelt)
    {
    case 1:
      fgets(line,BUFSIZ,infile);
      sscanf(line,"%lf",&tempfloat1);
      (rhs_temp)[iterelt] = (double)tempfloat1;
      iterelt++;
      break;
    case 2:
      fgets(line,BUFSIZ,infile);
      sscanf(line,"%lf %lf",&tempfloat1,&tempfloat2);
      (rhs_temp)[iterelt]   = (double)tempfloat1;
      (rhs_temp)[iterelt+1] = (double)tempfloat2;
      iterelt++;
      break;
    case 3:
      fgets(line,BUFSIZ,infile);
      sscanf(line,"%lf %lf %lf",&tempfloat1,&tempfloat2,&tempfloat3);
      (rhs_temp)[iterelt]   = (double)tempfloat1;
      (rhs_temp)[iterelt+1] = (double)tempfloat2;
      (rhs_temp)[iterelt+2] = (double)tempfloat3;
      iterelt++;
      break;
    default:
      break;
    }

  memcpy(csc->avals, values, edgeloc*sizeof(double));
  memFree_null(values);
  memcpy(*rhs, rhs_temp, vertloc*sizeof(double));
  memFree_null(rhs_temp);
  csc->n = vertloc;
  csc->flttype = PastixDouble;
  csc->mtxtype = PastixGeneral;
  csc->fmttype = PastixCSC;
  fclose(infile);
  free(filename);
  free(vectsize);
  free(vectsizercv);
}
