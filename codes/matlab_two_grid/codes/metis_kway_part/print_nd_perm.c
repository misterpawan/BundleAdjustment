/* Print the part array to be used in GAMG (Graph based Algebraic Multigrid) 

   Author : Pawan Kumar
   Date   : 1st December 2012
   Place  : Leuven
   email  : kumar.work@rediffmail.com

  This interface was written out of frustration where the mex interface to metis 
  needed gcc-4.2 to compile. The gcc-4.2 is no longer supported by some linux 
  distributation, it is annoying to do a separate install of gcc-4.2 manually.
  Moreover, when we only need output of a vector, why should we write the whole 
  mex interface to the often complicated input datastructures. Imagine classes, 
  array of structs or struct of arrays. It is simply core convenient to write 
  the output vector in a file and read the necessary input from a file. This 
  is programming language independent, platform independent if the external program 
  is platform independent too, for example C is.
  This code is run inside MATLAB using "system" command, the code writes the part 
  array which is then read by MATLAB. The author hopes that it could be useful 
  to those who need metis partitions in MATLAB and do not worry about metis interface!
  
  Currently, only interface to the 2-way ND is provided, but it could be modified 
  to interface with the other functions of METIS.

  The input file is supposed to be in .mtx format
  
 */

#include "print_part.h"

int main(int argc, char *argv[]){
  
  int        err;
  MATCOO     ACOO;
  MATCSC     ACSC, ACSC2, C;
  idx_t      *part, *perm, *iperm;
  idx_t     options[METIS_NOPTIONS];
  int        myoptions[18];
  idx_t        *sizes=0;
  idx_t    *xadj, *adjncy, *vwgt, *adjwgt, ncon, nn;
  int        i;
  idx_t    numflag, edgecut;
  idx_t      scheme, wgtflag;
  clock_t start, end;
  double cpu_time_used;

  if(argc != 22) 
    {
      printf("Usage: %s [matrix] [nparts] [wgtflag] [ptype] [objtype] [ctype] [iptype] [rtype] [ncuts] [nseps] [numb] [niter] [seed] [minconn] [contig] [compress] [ccorder] [pfactor] [ufactor] [dbglvl] [default]\n", argv[0]); 
      exit(1);
    } 

  FILE *matFile = fopen( argv[1], "r" );
  int nparts    = atoi(argv[2]);
  wgtflag       = atoi(argv[3]);
  
//  set_options(argv, options, myoptions);

  err   = readMatrixCoo(matFile, &ACOO);

  ACSC.M   = ACOO.M;
  ACSC.N   = ACOO.N;
  ACSC.nnz = ACOO.nnz;
  
  ACSC.I = (int*)malloc((ACOO.nnz)*sizeof(int));
  ACSC.J = (int*)malloc((ACOO.M+1)*sizeof(int));
  ACSC.V = (double*)malloc((ACOO.nnz)*sizeof(double));
  
  coocsc( ACOO.M, ACOO.nnz, ACOO.V, ACOO.I, ACOO.J, ACSC.V,
	   ACSC.I, ACSC.J );
  
  ACSC2.M   = ACSC.M;
  ACSC2.N   = ACSC.N;
  ACSC2.nnz = ACSC.nnz;
  ACSC2.I   = (int*)malloc((ACSC.nnz)*sizeof(int));
  ACSC2.J   = (int*)malloc((ACSC.M+1)*sizeof(int));
  ACSC2.V   = (double*)malloc((ACSC.nnz)*sizeof(double));

  coocsc( ACOO.M, ACOO.nnz, ACOO.V, ACOO.J, ACOO.I, ACSC2.V,
	   ACSC2.I, ACSC2.J ); 

  Brahm_Add_Csc(&ACSC, &ACSC2, &C, 1.0, 1.0);
  C.nnz = C.J[C.N]; 

  free(ACSC.I); free(ACSC.J); free(ACSC.V);
  free(ACSC2.I); free(ACSC2.J); free(ACSC2.V);

  xadj = (idx_t*) calloc (ACOO.M+1, sizeof(idx_t));
  adjncy = (idx_t*) calloc (ACOO.nnz, sizeof(idx_t));
  vwgt = (idx_t*) calloc (ACOO.M+1, sizeof(idx_t));
  adjwgt = (idx_t*) calloc (ACOO.nnz, sizeof(idx_t)); 
  part = (idx_t*) calloc (ACOO.M, sizeof(idx_t));
  perm = (idx_t*) calloc (ACOO.M, sizeof(idx_t));
  iperm = (idx_t*) calloc (ACOO.M, sizeof(idx_t));
  sizes = (idx_t*) calloc (2*nparts, sizeof(idx_t));

  csc2metisgraph3 (C.M, C.J[C.N], C.I, C.J, C.V,
		   xadj, adjncy, vwgt, adjwgt);

  free(C.I); free(C.J); free(C.V);

  ncon = 1; nn = (idx_t) ACOO.M;

  start = clock();
  switch (wgtflag)
    {
    case 0:
      metis_partgraphrecursive(&nn, &ncon, xadj, adjncy, NULL, NULL, NULL, &nparts, NULL, NULL, options, &edgecut, part);
      break;
    case 1:
      metis_partgraphrecursive(&nn, &ncon, xadj, adjncy, vwgt, NULL, NULL, &nparts, NULL, NULL, options, &edgecut, part);
      break;
    case 2:
      metis_partgraphrecursive(&nn, &ncon, xadj, adjncy, vwgt, NULL, adjwgt, &nparts, NULL, NULL, options, &edgecut, part);
      break;
    case 3:
      METIS_NodeNDP(nn, xadj, adjncy, NULL, nparts, options, perm, iperm, sizes);
      break;

    }
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("Metis Aggregation: %g seconds\n", cpu_time_used);


  FILE *fp = fopen( "perm_array.txt", "w" );
  FILE *fp2 = fopen( "size_array.txt", "w" );

  for (i = 0; i < ACOO.M; i++)
    fprintf(fp, "%d\n", (int)(perm[i]+1));

  for (i = 0; i < 2*nparts-1; i++)
    fprintf(fp2, "%d\n", (int)sizes[i]);

  fclose(fp);
  free(xadj); free(adjncy); 

  return 0;

}