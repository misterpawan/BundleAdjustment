// function to add two csc matrix

// #include "brahm_defs.h"
#include "print_part.h"

void Brahm_Add_Csc(MATCSC *A, MATCSC *B, MATCSC *C, double alpha, double beta){

  int p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w; // values=0 ;
  double *x, *Bx, *Cx ;
  
  if(A->M != B->M || A->N != B->N ) {
    printf("MAtrices are not conformable \n");
    exit(1);
  }

  m = A->M ; anz = A->nnz ;
  n = B->N ; Bp = B->J ; Bx = B->V ; bnz = B->nnz ;

  w = (int*)calloc(m, sizeof(int));
  x = (double*)calloc(m, sizeof(double));

  C->M   = A->M;
  C->N   = A->N;
  C->nnz = A->nnz + B->nnz;
  C->I   = (int*)malloc(C->nnz*sizeof(int));
  C->J   = (int*)malloc((C->N+1)*sizeof(int));
  C->V   = (double*)malloc(C->nnz*sizeof(double));

  Cp = C->J ; Ci = C->I ; Cx = C->V ;

  for (j = 0 ; j < n ; j++)
    {
      Cp [j] = nz ;                   /* column j of C starts here */
      nz = Brahm_scatter (A, j, alpha, w, x, j+1, C, nz) ;   /* alpha*A(:,j)*/
      nz = Brahm_scatter (B, j, beta, w, x, j+1, C, nz) ;    /* beta*B(:,j) */
      //      if (values) 
      for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
    }
  Cp [n] = nz ;                       /* finalize the last column of C */
  
}
