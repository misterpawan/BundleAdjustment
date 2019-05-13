/* The definitions file for print_part.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "metis.h"
#include "mmio.h"
#include <time.h>

//typedef int idx_t;
//typedef int idxtype idx_t; 

typedef struct {
  int nnz;
  int M;     // row dim
  int N;     // col dim
  int *I;    // row indices
  int *J;    // col indices
  double *V; // values
} MATCOO  ;

typedef struct {
  int nnz;
  int M;     // row dim
  int N;     // col dim
  int *I;    // row indices
  int *J;    // col indices
  double *V; // values
} MATCSC ;

int coocsc3( int ncol, int nnz, double *val, int *rowind, int *colind, double *val_new, int *rowind_new, int *colptr);

void Brahm_Add_Csc(MATCSC *A, MATCSC *B, MATCSC *C, double alpha, double beta);


//int csc2metisgraph3 (int n, int nnz, int *ir, int *jc, double *pr, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt);

int csc2metisgraph3 (int n, int nnz, int *ir, int *jc, double *pr, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, idx_t *adjwgt);

int readMatrixCoo( FILE *f, MATCOO *ACOO );

int Brahm_scatter(MATCSC *A, int j, double beta, 
		  int *w, double *x, int mark,
		  MATCSC *C, int nz);

void Interface_Metis_KWay_Part(int npes, MATCOO *ACOO,
			       idx_t *part );

