// Interface to Metis k-way partitioning and reordering the block rows according to
// the partitions, P. Kumar, March 2012

// #include "defs_for_metis.h"
// #include "protos_for_metis.h"

#include "print_part.h"


void Interface_Metis_KWay_Part(int npes, MATCOO *ACOO,
					    idx_t *part ){

  MATCSC     ACSC2;
  MATCSC     ACSC;
  MATCSC     C;
  idx_t    *xadj, *adjncy, *vwgt, *adjwgt;
  idx_t        *perm;
	idx_t      *iperm;
  idx_t        options[METIS_NOPTIONS];
  idx_t        edgecut;
  //   int        i;
  int        wgtflag, numflag;

  ACSC.M = ACOO->M;
  ACSC.N = ACOO->N;
  ACSC.nnz = ACOO->nnz;

  ACSC.I = (int*)malloc((ACOO->nnz)*sizeof(int));
  ACSC.J = (int*)malloc((ACOO->M+1)*sizeof(int));
  ACSC.V = (double*)malloc((ACOO->nnz)*sizeof(double));

  coocsc3( ACOO->M, ACOO->nnz, ACOO->V, ACOO->I, ACOO->J, ACSC.V,
	   ACSC.I, ACSC.J );

  ACSC2.M   = ACSC.M;
  ACSC2.N   = ACSC.N;
  ACSC2.nnz = ACSC.nnz;
  ACSC2.I   = (int*)malloc((ACSC.nnz)*sizeof(int));
  ACSC2.J   = (int*)malloc((ACSC.M+1)*sizeof(int));
  ACSC2.V   = (double*)malloc((ACSC.nnz)*sizeof(double));

  coocsc3( ACOO->M, ACOO->nnz, ACOO->V, ACOO->J, ACOO->I, ACSC2.V,
	   ACSC2.I, ACSC2.J );

  Brahm_Add_Csc(&ACSC, &ACSC2, &C, 1.0, 1.0);
  C.nnz = C.J[C.N];

  free(ACSC.I); free(ACSC.J); free(ACSC.V);
  free(ACSC2.I); free(ACSC2.J); free(ACSC2.V);

  xadj = (idx_t*) calloc (ACOO->M+1, sizeof(idx_t));
  adjncy = (idx_t*) calloc (ACOO->nnz, sizeof(idx_t));
  vwgt = (idx_t*) calloc (ACOO->M+1, sizeof(idx_t));
  adjwgt = (idx_t*) calloc (ACOO->nnz, sizeof(idx_t));
  //perm = (idx_t*) calloc (ACOO->M, sizeof(idx_t));
  //iperm = (idx_t*) calloc (ACOO->M, sizeof(idx_t));

  csc2metisgraph3 (C.M, C.J[C.N], C.I, C.J, C.V,
		   xadj, adjncy, vwgt, adjwgt);

  free(C.I); free(C.J); free(C.V);

  //options[0]=0; // default options for metis
	METIS_SetDefaultOptions(options);
	//options[METIS_OPTION_NCUTS] = 1;

  /* added on 19th march 2012, P. Kumar, Leuven */

  wgtflag = 0; // metis_partgraphkway needs pointers to these flags, hence extra variables :(
  numflag = 0;

  idx_t ncon;
	ncon = 1;
	idx_t nvtxs = (idx_t)ACOO->M;
	idx_t nparts = (idx_t) npes;

  //printf("OK before k part\n");
	//printf("nparts = %ld\n", nparts);

  METIS_PartGraphKway(&nvtxs, &ncon, xadj, adjncy, NULL, NULL, NULL, &nparts, NULL, NULL, options, &edgecut, part);
  //   for(i=0;i<10;i++) printf("%d\n", part[i]);

 //printf("OK after k part\n");

  //   METIS_NodeNDP(ACOO->M, xadj, adjncy, npes,options, perm, iperm, part) ;

  free(xadj); free(adjncy);

// free(iperm); free(perm);
free(adjwgt);
free(vwgt);

  //  Brahm_Permute_Coo( ACOO, perm); free(perm);
}