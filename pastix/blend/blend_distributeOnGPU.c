#include "blend_distributeOnGPU.h"
#include "flops.h"

#define symbol_get_cblk_stride(  _datacode_, _cblknum_ ) (SOLV_STRIDE((_cblknum_)))
#define symbol_get_cblk_width(   _datacode_, _cblknum_ ) (SYMB_LCOLNUM((_cblknum_))-SYMB_FCOLNUM((_cblknum_))+1)
#define symbol_get_blok_coefind( _datacode_, _bloknum_ ) (SOLV_COEFIND((_bloknum_)))
#define symbol_get_blok_height(  _datacode_, _bloknum_ ) (BLOK_ROWNBR((_bloknum_)))
//#include "fifo.h"

typedef struct _cblk_elem_ {
  z_SolverCblk       *cblk;
  double criterium;
  int cblk_size;
  int updates;
  float flop; /* in Mflop */
}cblk_elem;

typedef struct _gpu_elem_ {
  int current_avail_mem;
  int gpu_total_mem;
  int id;
}gpu_elem;


int blend_distributeOnGPU(z_SolverMatrix * datacode,
                          double         maxMem,
                          int            pageSize,
                          int            criterium,
                          enum API_GPU_CRITERIUM nGPUs,
                          enum API_FLOAT floatType,
                          enum API_FACT  factType) {
  z_SolverCblk   *cblktab = datacode->cblktab;
  z_UpDownVector * updovect = &(datacode->updovct);
  PASTIX_INT   cblknbr, cblknum;
  cblk_elem * cblkInfo;
  gpu_elem  * gpuInfo;
  pastix_queue_t cblkQueue;
  pastix_queue_t gpuQueue;
  PASTIX_INT gcblk2list   = -1;
  PASTIX_INT ndevices;
  int i;
  PASTIX_INT j, bloknum, m, k, n;
  float flop;
  gpu_elem *new_gpu;
  cblk_elem *new_cblk;
  gpu_elem *gpu_cursor;
  cblk_elem *cblk_cursor;
  int type_sze;
  size_t unit_size;
  ndevices = nGPUs;

  switch (floatType) {
  case API_REALSINGLE:
    type_sze=sizeof(float);
    break;
  case API_REALDOUBLE:
    type_sze=sizeof(double);
    break;
  case API_COMPLEXSINGLE:
    type_sze=2*sizeof(float);
    break;
  case API_COMPLEXDOUBLE:
    type_sze=2*sizeof(double);
    break;
  default:
    errorPrint("Unkwnown type");
    return FLOAT_TYPE_ERR;
  }
  unit_size = pageSize/type_sze;

  int* devices_cblk = malloc((1+ndevices)*sizeof(int));
  int* devices_gemm = malloc((1+ndevices)*sizeof(int));
  memset(devices_cblk,0, (1+ndevices)*sizeof(int));
  memset(devices_gemm,0, (1+ndevices)*sizeof(int));

  /*
   * Compute the GPU distribution
   */
  MALLOC_INTERN(cblkInfo, SYMB_CBLKNBR, cblk_elem);
  MALLOC_INTERN(gpuInfo, ndevices, gpu_elem);
  pqueueInit(&cblkQueue, SYMB_CBLKNBR);
  pqueueInit(&gpuQueue, ndevices);

  /* Sort the cblk according to criterium */
  cblknbr = SYMB_CBLKNBR;
  //fprintf(stdout,"start loop on  cblk \n");
  for(cblknum = 0; cblknum < cblknbr; cblknum++){
    gcblk2list = UPDOWN_GCBLK2LIST(UPDOWN_LOC2GLOB(cblknum));

    /* if it's a leaf, gpuid = cpu */
    if ( gcblk2list == -1 || ndevices <=0) {
      cblktab[cblknum].gpuid = -1;
      devices_cblk[0]++;
    } else {
      size_t cblksize;
      int    nbpages;
      PASTIX_INT updates = UPDOWN_LISTPTR(gcblk2list+1)
        -                 UPDOWN_LISTPTR(gcblk2list);

      cblktab[cblknum].gpuid = -2;
      /* If we are here, we have at least one update */
      assert(updates > 0);
      new_cblk = &(cblkInfo[cblknum]);
      new_cblk->updates = updates;
      new_cblk->cblk = cblktab + cblknum;

      /* FLOP */
      int cblknum2;
      flop = 0.0;
      //fprintf(stdout,"start flops computation for cblk %ld \n", cblknum);
      for(j=UPDOWN_LISTPTR(gcblk2list); j<UPDOWN_LISTPTR(gcblk2list+1);j++){
        bloknum = UPDOWN_LISTBLOK( j );
        cblknum2 = updovect->listcblk[j];
        //cblknum2 = sparse_matrix_get_lcblknum(datacode, bloknum);
        m = symbol_get_cblk_stride( datacode, cblknum2 ) - symbol_get_blok_coefind(datacode, bloknum);
        k = symbol_get_cblk_width( datacode, cblknum2 );
        n = symbol_get_blok_height( datacode, bloknum );
        flop += (float)FLOPS_SGEMM(m,n,k)/(float)(2*1e6);
      }
      //fprintf(stdout,"stop flops computation for cblk %ld \n", cblknum);

      /* if (floatType == API_COMPLEXSINGLE || floatType == API_COMPLEXDOUBLE) { */
      /*   flop *= 4.0; */
      /* } */

      /* if (factType == API_FACT_LU){ *\/ */
      /*   flop *= 2.0; */
      /* } */

      new_cblk->flop = flop;
      /* check for int overflow */
      assert(flop < 4294967296.);

      /* Amount of memory to push for the cblk */
      cblksize = symbol_get_cblk_width( datacode, cblknum ) * symbol_get_cblk_stride( datacode, cblknum );
      nbpages = ( cblksize + unit_size - 1) / unit_size;

      if (factType == API_FACT_LU){
          nbpages *= 2.0;
      }

      new_cblk->cblk_size = nbpages;
      /* Sort criterium */
      switch(criterium){
      case API_GPU_CRITERION_UPDATES:
        new_cblk->criterium = (double)-new_cblk->updates;
        break;
      case API_GPU_CRITERION_CBLKSIZE:
        new_cblk->criterium = (double)-new_cblk->cblk_size;
        break;
      case API_GPU_CRITERION_FLOPS:
        new_cblk->criterium = (double)-new_cblk->flop;
        break;
      default :
        new_cblk->criterium = (double)cblknum;
      }
      pqueuePush1(&cblkQueue, cblknum, new_cblk->criterium);
    }
  }
  // fprintf(stdout,"end loop on  cblk \n");

  /* Sort the GPUs according to available memory */
  for(i = 0; i < ndevices; i++) {
    new_gpu = gpuInfo+i;
    //new_gpu->gpu_device = gpu_enabled_devices[i];
    new_gpu->id = (PASTIX_INT)i;
    new_gpu->gpu_total_mem     = (int)(GPU_MAX_FILL * maxMem /pageSize);
    //new_gpu->gpu_total_mem     =  new_gpu->gpu_device->memory->max_segment);
    new_gpu->current_avail_mem = new_gpu->gpu_total_mem;
    pqueuePush1(&gpuQueue, i, -new_gpu->current_avail_mem);
  }

    /* Create association between cblks and gpus */
  if(ndevices > 0){
    while (! pqueueSize(&cblkQueue) == 0 ) {
      pastix_int_t cblknum = pqueuePop(&cblkQueue);
      pastix_int_t gpuid;
      pastix_int_t gpu_add;
      gpu_add = 0;
      gpu_cursor = NULL;

      /* Get the cblk with the highest criterium */
      cblk_cursor = cblkInfo + cblknum;

      /* The cblk has a predicted GPU */
      gpuid = cblk_cursor->cblk->gpuid;

      if (gpuid != -2) {
        assert( gpuid != -1 );
        gpu_cursor = gpuInfo + gpuid;
        if(gpu_cursor != NULL){
          if(gpu_cursor->current_avail_mem < cblk_cursor->cblk_size){
            gpu_cursor = NULL;
          }
        }
      }

      /* Get the gpu with the highest criterium */
      if (gpu_cursor == NULL) {
        gpuid = pqueuePop(&gpuQueue);
        gpu_add = 1;
        gpu_cursor  = gpuInfo + gpuid;
      }

      //fprintf(stdout, "Pop: cblk(%ld) size=%ld, nbpages=%d, updates=%d\n",
      //(cblk_cursor->cblk) - cblktab,
      //(cblk_cursor->cblk->lcolnum - cblk_cursor->cblk->fcolnum + 1) * cblk_cursor->cblk->stride,
      //cblk_cursor->cblk_size, cblk_cursor->updates );
      //fprintf(stdout, "Pop: gpu(%d) available mem=%d\n",
      //gpu_cursor->id, gpu_cursor->current_avail_mem);

      assert(gpu_cursor != NULL);
      if( (gpu_cursor->current_avail_mem >= cblk_cursor->cblk_size)
          && (gpu_cursor->current_avail_mem > 0)
          && (cblk_cursor->updates >= GPU_MIN_UPDATES)
          && (cblk_cursor->cblk_size >= GPU_MIN_NBPAGES)
          && (cblk_cursor->flop >= GPU_MIN_FLOP)
          ) {
        /* association */
        gpu_cursor->current_avail_mem -= cblk_cursor->cblk_size;
        cblk_cursor->cblk->gpuid = gpu_cursor->id;

        /* Add prediction to contibutors */
        for(j=UPDOWN_LISTPTR(gcblk2list); j<UPDOWN_LISTPTR(gcblk2list+1);j++){
          bloknum = UPDOWN_LISTBLOK( j );
          cblknum = updovect->listcblk[j];
          //cblknum = sparse_matrix_get_lcblknum(datacode, bloknum);
          if(cblktab[cblknum].gpuid == -2 && gcblk2list != -1){
            cblktab[cblknum].gpuid = gpu_cursor->id;
          }
        }
      } else {
        /* send to CPU */
        cblk_cursor->cblk->gpuid = -1;
      }

      if (gpu_add==1)
        pqueuePush1(&gpuQueue, gpuid, -gpu_cursor->current_avail_mem);

      devices_cblk[cblk_cursor->cblk->gpuid+1]++;
      devices_gemm[cblk_cursor->cblk->gpuid+1]+=cblk_cursor->updates;
    }
  }

  fprintf(stdout,"cpu :  %d cblk, %d gemm\n", devices_cblk[0], devices_gemm[0]);

  for (i = 0 ; i < ndevices; i++) {
    gpu_cursor  = gpuInfo + i;

    fprintf(stdout,"gpu %d:  %d cblk, %d gemm, memory : %d / %d \n",
            i, devices_cblk[i+1], devices_gemm[i+1],
            gpu_cursor->current_avail_mem, gpu_cursor->gpu_total_mem );
  }

  pqueueExit(&cblkQueue);
  pqueueExit(&gpuQueue);
  return cblknbr;
}
