/* Two=Grid Solve as a preconditioner for GMRES */


// C++ libraries
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstring>
#include <ctime>
#include <cmath>

using namespace std;    // C++ standard libary namespace

// umfpack libraries
#include "cs.h"
#include "umfpack.h"
#include "util.h"

//mkl libraries
#include "mkl.h"
#include "mkl_solvers_ee.h"

// macros
#define max(a, b) (a) < (b) ? (b): (a)  // for eigenvalue solve routine
#define size 128                       // for ILU 
#define sizeG 441 //hardcoding it for the time being



int main()
{
  

  double* JTe; //rhs
  int i;
  int num_rows;
  int num_cols; // no of columns
  int ncc = 0; // no of non zeros
  int *null = ( int * ) NULL;
  double *solve_null = ( double * ) NULL;
  void *Numeric, *Numeric_D,*Numeric_MSC;
  string prefix = "JTJ49_1";
  double r;
  int status,sym_status,num_status,solve_status;
  void* Symbolic,*Symbolic_D,*Symbolic_MSC;
  int sizeD;
  //int sizeG;
  int sizeB;   // Block Diagonal size
  int j,k;
  int nzD = 0;
  int nzG = 0;
  int nzL = 0; // nzL = nzU
  int nzB = 0; // Non zero in Block Diagonal
  int iterD = 0;
  int iterL = 0;
  int iterG = 0;
  int iterB = 0; 

  //creating structures of cs_diPARSE
  cs_di* A = new cs_di;
  cs_di* MSC = new cs_di;
  cs_di* D = new cs_di;
  cs_di* L = new cs_di;
  cs_di* G = new cs_di;
  cs_di* U = new cs_di;
  cs_di* B = new cs_di;


  //string line;
  string rhs_filename = "JTe49_1.txt";


  // getting the matrix size: n -> no of cols, ncc -> no of non zeros
  cc_header_read ( prefix, ncc, num_cols );
  num_cols = num_cols+1;                          //DOUBT!!!!!!!!
  cout << "\nNo of non zeros = "<< ncc << "\n";
  cout << "\nNo of columns = "<< num_cols << "\n";

  //size of the D matrix
  sizeD = num_cols - sizeG; 
  //size of the B matrix
  sizeB = num_cols;

  A->nzmax = ncc;
  A->m = num_cols;
  A->n = num_cols; // since square matrix
  A->nz = -1;

  //Allocate space for rhs
  JTe = new double[num_cols];

  // reading the rhs
  r8vec_data_read ( rhs_filename, num_cols, JTe);


  //Allocate space for the coefficient matrix
  A->x = new double[ncc];
  A->p = new int[num_cols+1];
  A->i = new int[ncc];

  //read the matrix data
  cc_data_read ( prefix, ncc, num_cols, A->i, A->p, A->x );
  A->p[num_cols] = ncc;

  cout << "\nFile read complete!\n";
  /**************File Read in CSC format Complete***************/


  /**********Converting A from CSC to CSR*****************/
  MKL_INT ivar;
  char cvar;
  ivar = num_cols;
  cvar = 'N'; //no transpose
  MKL_INT job[6] = {1,1,0,0,0,1};    // various parameters for conversion
  double *acsr = new double[ncc];    // values
  MKL_INT *ja = new MKL_INT[ncc];    // column
  MKL_INT *ia = new MKL_INT[ivar+1]; // rows
  MKL_INT info1;                      // status

  // call to routine dcsrcsr
  mkl_dcsrcsc(job,&ivar,acsr,ja,ia,A->x,A->i,A->p,&info1);

  cout << "\n Conversion info A : "<< info1 << "\n";

  /**************Matrix Conversion in CSR format Complete***************/


  /***********Now Call MKL Eigenvalue Solve Routine************/


/*
!   Content: Example for k Max/Min eigenvalue problem based on Intel MKL 
!            Extended Eigensolver (CSR sparse format, double precision)
!
!*******************************************************************************
!
! The following routines are used in the example:
!          MKL_SPARSE_D_EV
!
! Consider the 4x4 matrix A
!
!                 |  6   2   0   0   |
!                 |  2   3   0   0   |
!     A   =       |  0   0   2  -1   |
!                 |  0   0  -1   2   |
!
! stored as sparse matrix.
!
!
!  The test calls mkl_sparse_d_ev routine to find several largest singular 
!  values and corresponding right-singular vectors. Orthogonality of singular  
!  vectors is tested using DGEMM routine
!
!*******************************************************************************/

 /* Matrix A of size N in CSR format */
  MKL_INT N = num_cols;               /* number of rows in matrix A */
  MKL_INT M = num_cols;               /* number of columns in matrix A , since symmetric # cols = # rows*/
  MKL_INT nnz = ncc;             /* number of non-zeros in matrix */

 
//MKL_INT ia[5] = {1,3,5,7,9};                         /* ia array from CSR format */
//MKL_INT ja[8] = {1,2,1,2,3,4,3,4};                   /* ja array from CSR format */
//double   a[8] = {6.0,2.0,2.0,3.0,2.0,-1.0,-1.0,2.0}; /* val array from CSR format */

//double   Eig[4] = {1.0, 2.0, 3.0, 7.0}; /* Exact eigenvalues */

/* mkl_sparse_d_ev input parameters */
  char         which = 'L'; /* Which eigenvalues to calculate. ('L' - largest (algebraic) eigenvalues, 'S' - smallest (algebraic) eigenvalues) */
  MKL_INT      pm[128];     /* This array is used to pass various parameters to Extended Eigensolver Extensions routines. */
  MKL_INT      k0  = 5;     /* Desired number of max/min eigenvalues */   
  MKL_INT      k1  = k0*num_cols;  /* array size for Eigenvalues */

  /* mkl_sparse_d_ev output parameters */        
  MKL_INT      kk;           /* Number of eigenvalues found (might be less than k0). */    
  double       E[k0];        /* Eigenvalues */
  double       X[k1];        /* Eigenvectors */
  double       res[k0];      /* Residual */
 
  /* Local variables */
  MKL_INT      info;               /* Errors */
  MKL_INT      compute_vectors = 1;/* Flag to compute eigenvectors */
  MKL_INT      tol = 7;            /* Tolerance */
  double       Y[k0];               /* Y=(X')*X-I */
  double       sparsity;           /* Sparsity of randomly generated matrix */
  MKL_INT      i1, j1;
  double       smax, t;    
    
  /* Sparse BLAS IE variables */
  sparse_status_t status1;
  sparse_matrix_t AA = NULL; /* Handle containing sparse matrix in internal data structure */
  struct matrix_descr descr; /* Structure specifying sparse matrix properties */

  /* Create handle for matrix A stored in CSR format */
  descr.type = SPARSE_MATRIX_TYPE_GENERAL; /* Full matrix is stored */
  status1 = mkl_sparse_d_create_csr ( &AA, SPARSE_INDEX_BASE_ONE, N, N, ia, ia+1, ja, acsr );

  /* Step 2. Call mkl_sparse_ee_init to define default input values */
  mkl_sparse_ee_init(pm);

  pm[1] = tol; /* Set tolerance */
  pm[6] = compute_vectors; /* Compute Eigenvectors */

  /* Step 3. Solve the standard Ax = ex eigenvalue problem. */
  info = mkl_sparse_d_ev(&which, pm, AA, descr, k0, &kk, E, X, res);   

  printf("mkl_sparse_d_ev output info %d \n",info);
  if ( info != 0 )
  {
      printf("Routine mkl_sparse_d_ev returns code of ERROR: %i", (int)info);
      return 1;
  }

  printf("*************************************************\n");
  printf("************** REPORT ***************************\n");
  printf("*************************************************\n");
  printf("#mode found/subspace %d %d \n", kk, k0);
  printf("Index/Exact Eigenvalues/Residuals\n");
 

  for (i1=0; i1<kk; i1++)
 {
     printf("   %d  %.15e %.15e \n" ,i1, E[i1], res[i1]);
 }
 
  mkl_sparse_destroy(AA);  // free memory for csr matrix handle

/*******Eigen value/vector Computation Complete *********/


/****Coarse Grid Matrix = Vector of Eigenvalues AND Prolongation Matrix = Eigenvectors****/


//   Construct Smoother: D = blkdiag(A(1:nJ1, 1:nJ1), A(nJ+1:n, nJ+1:n))

/**********Domain Decomposition************/


  //Allocating memory for blocks
  //MSC->p = new int[sizeG+1];
  D->p = new int[sizeD+1];
  L->p = new int[sizeD+1];
  G->p = new int[sizeG+1];  
  U->p = new int[sizeG+1];
  B->p = new int[sizeB+1];


  //MSC->nz = -1;MSC->m = sizeG;MSC->n = sizeG;
  D->nz = -1;D->m = sizeD;D->n = sizeD;
  G->nz = -1;G->m = sizeG;G->n = sizeG;
  L->nz = -1;L->m = sizeG;L->n = sizeD;
  U->nz = -1;U->m = sizeD;U->n = sizeG;
  B->nz = -1;B->m = sizeB;B->n = sizeB;  


  //cout << "\nCounting non zeros ...\n";
  for(j=0;j<sizeD;j++)
  {
    for(k=A->p[j]; k < A->p[j+1]; k++)
    {
      if(A->i[k] < sizeD){ 
        ++nzD;
        //++nzB;
      }
      else ++nzL;
    }
  }

  
  nzG = ncc - (nzD + 2*nzL); 
  nzB = nzD + nzG;

  //Allocating memory
  D->i = new int[nzD];
  D->x = new double[nzD];
  L->i = new int[nzL];
  L->x = new double[nzL];
  U->i = new int[nzL];
  U->x = new double[nzL];
  G->i = new int[nzG];
  G->x = new double[nzG];
  B->i = new int[nzB];   // num of cols = nzB
  B->x = new double[nzB];   // nzB = nzD + nzG

  //MSC->i = new int[sizeG*sizeG];
  //MSC->x = new double[sizeG*sizeG];

  //setting values
  D->nzmax = nzD;
  L->nzmax = nzL; 
  U->nzmax = nzL;
  G->nzmax = nzG;  
  B->nzmax = nzB;
  //MSC->nzmax = sizeG*sizeG;

  
  cout << "\nFilling non zeros ...\n";
  for(j=0;j<sizeD;j++)
  {
    D->p[j] = iterD;
    L->p[j] = iterL;
    B->p[j] = iterD;

    for(k=A->p[j]; k < A->p[j+1]; k++)
    {
      if(A->i[k] < sizeD) 
      {
        D->i[iterD] = A->i[k];
        D->x[iterD] = A->x[k];
        B->i[iterD] = A->i[k];     // Filling B block upto D blocks
        B->x[iterD] = A->x[k];
        iterD += 1;
        iterB += 1;
        //++nzD_col;
      }
      else
      {
        L->i[iterL] = A->i[k] - sizeD;
        L->x[iterL] = A->x[k];
        iterL += 1;
      }
    }
  }
  D->p[sizeD] = iterD;
 // B->p[sizeD] = iterB;    // Doubt !!!, iterB = iterD
  L->p[sizeD] = iterL;
  
  status = umfpack_di_transpose(sizeG,sizeD, L->p,L->i,L->x,null,null,U->p,U->i,U->x) ;
  //cout << "\n TRANSPOSE STATUS : "<< status<< "\n";
  

  for(j = sizeD;j<num_cols;j++)
  {
    G->p[j - sizeD] = iterG;
    B->p[j] = iterB;            // as B already filled upto D
    if(j < num_cols-1)
    {
      //cout << "\n In first condition...\n";
      for(k = A->p[j]; k < A->p[j+1];k++)
      {
        if(A->i[k] < sizeD) continue;
        else
        {
          //cout << "\n row : " << A->i[k] << "\n";
          G->i[iterG] = A->i[k]-sizeD;
          G->x[iterG] = A->x[k];
          B->i[iterB] = A->i[k];     // Don't subtract sizeD
          B->x[iterB] = A->x[k];
          iterG += 1;
          iterB += 1;
        }     
      } 
    }
    //for the last column
    else 
    {
      //cout << "\n In second condition..."<<A->nzmax <<"\n";
      for(k = A->p[j]; k < A->nzmax;k++)
      {
        if(A->i[k] < sizeD) continue;
        else
        {
          //cout << "\n row : " << A->i[k] << "\n";
          G->i[iterG] = A->i[k] - sizeD;
          G->x[iterG] = A->x[k];
          B->i[iterB] = A->i[k];     // Don't Subtract sizeD
          B->x[iterB] = A->x[k];
          iterG += 1;
          iterB += 1;
        }
        //cout << "\n k : " << k << "\n";
        //break;
      }
    }
      //break;
  }

  G->p[sizeG] = iterG;
  B->p[sizeB] = iterB;    // iterB = iterD + iterG


  cout << "\nFilling non zeros complete!!\n";

  
  /**********Convert B block to CSR****************/

  MKL_INT ivar1;
  char cvar1;
  ivar1 = sizeB;
  cvar1 = 'N'; //no transpose
  MKL_INT job1[6] = {1,1,0,0,0,1};         // various parameters for conversion
  double *bacsr = new double[B->nzmax];    // values
  MKL_INT *bja = new MKL_INT[B->nzmax];    // column
  MKL_INT *bia = new MKL_INT[sizeB];       // rows
  MKL_INT info2;                           // status

  // call to routine dcsrcsr
  mkl_dcsrcsc(job1,&ivar1,bacsr,bja,bia,B->x,B->i,B->p,&info2);
  //mkl_dcsrcsc(job,&ivar,acsr,ja,ia,A->x,A->i,A->p,&info1);

  cout << "\n Conversion info B : "<< info2 << "\n";

/*************Matrix Conversion of B block in CSR format Complete!!!*************/

/**********Compute LU factorization of Smoother BLK-diagonal Matrix - B*********/

/*

  dcsrilu0 (const MKL_INT *n , const double *a , const MKL_INT *ia , const MKL_INT *ja , double *bilu0 , const MKL_INT *ipar , const double *dpar , MKL_INT *ierr );
*/

  // paramters for ilu routine
  MKL_INT ipar[size];
  double dpar[size];
  MKL_INT ierr = 0;
  double *bilu0 = new double[B->nzmax];

  // Initialize the parameters
  ipar[30] = 1;
  dpar[30] = 1.E-20;  // 1.0e-16
  dpar[31] = 1.E-16;

  // call to ILU routine for B
  dcsrilu0 (&ivar1, bacsr, bia, bja, bilu0, ipar, dpar, &ierr);

  printf ("dcsrilu0 has returned the code %d\n", ierr);

  /*******Approximate LU factorization complete using ILU**********/


  /* Two-grid Solve Function */
  // As of now just do in main build function later 

  /* 1st Step : Pre-smoothing  */

  printf("   Solve t = U\\(L\\x) system   \n");

  // Paramteres for triangular solver routine - mkl_cspblas_scsrtrsv

  double *t_0 = new double[num_cols];
  double *t_1 = new double[num_cols];

  uplo = 'l';        // lower triangular
  nonunit = 'n';     // diagonal is not unit triangular
  transa = 'n';      // not tranpose

  mkl_dcsrtrsv(&uplo,&transa,&nonunit,&ivar1,bilu0,bia,bja,&tmp[ipar[21]-1],t_0);
  /* Computation of L\X complete: L\x = t_0  */

  uplo = 'u';        // upper triangular
  nonunit = 'n';     // diagonal is not unit triangular
  transa = 'n';      // not tranpose

  /* Compute U\(t_0); t_0 = L\x */
  mkl_dcsrtrsv(&uplo,&transa,&nonunit,&ivar1,bilu0,bia,bja,t_0,t_1);
  /*  Solution obtained in t_1  */

  printf(" Pre-smoothing Complete!!!");


  /* 2nd Step: Coarse Grid Solve:  g = P* (Ac \ (P'*x))  */

  //1st perform the innermost matrix vector multiplication i.e P'*x

  double *x_0 = new double[k];  // size = number of eigenvalues

  MKL_INT ldy  = k;     /* Leading dimension for destination array in GEMM */
  MKL_INT ldx  = N;    /* Leading dimension for source arrays in GEMM */
  MKL_INT Cols = 1;
  double one   = 1.0;   /* alpha parameter for GEMM */
  double zero  = 0.0;   /* alpha parameter for GEMM */
  dgemm(
        &DGEMMC,          /* IN: 'T', transposed case*/
        &DGEMMN,          /* IN: 'N', non-transposed case*/
        &k,               /* IN: Number of rows in matrix Y */
        &Cols,            /* IN: Number of columns in matrix Y */
        &N,               /* IN: Number of rows in matrix X */
        &one,             /* IN: alpha = 1.0 */
        X,                /* IN: Source #1 for GEMM, will be transposed */
        &ldx,             /* IN: Leading dimension of Source 1 */
        &tmp[ipar[21]-1], /* IN: Source #2 for GEMM */
        &ldx,             /* IN: Leading dimension of Source 2 */
        &zero,            /* IN: beta = 0.0 */
        x_0,                /* OUT: Destination */
        &ldy              /* IN: Leading dimension of Destination */
        );

  // Performs pointwise element by element division of vector x_0 by vector Ac
  //vsDiv( n, a, b, y );
  double *x_1 = new double[k];
  vsDiv(k, x_0, E, x_1);  // number of elements = k (number of eigenvalues) 

  // Computed x_1 = (Ac \ (P'*x)); x_0 = P'*x

  /*   
    Now compute matrix -vector multiplication :g = P * x_1 ; x1 = (Ac \ (P'*x));
  */

  double *x_2  = new double[num_cols]; 
  MKL_INT ldy1 = N;     /* Leading dimension for destination array in GEMM */
  MKL_INT ldx1 = k;    /* Leading dimension for source arrays in GEMM */

  dgemm(
        &DGEMMN,          /* IN: 'N', non-transposed case*/
        &DGEMMN,          /* IN: 'N', non-transposed case*/
        &N,               /* IN: Number of rows in matrix Y */
        &Cols,            /* IN: Number of columns in matrix Y */
        &N,               /* IN: Number of rows in matrix X */
        &one,             /* IN: alpha = 1.0 */
        X,                /* IN: Source #1 for GEMM, will be transposed */
        &ldx,             /* IN: Leading dimension of Source 1 */
        &x_1, /* IN: Source #2 for GEMM */
        &ldx1,             /* IN: Leading dimension of Source 2 */
        &zero,            /* IN: beta = 0.0 */
        x_2,                /* OUT: Destination */
        &ldy1              /* IN: Leading dimension of Destination */
        );

  // Compute : q(x_3) = A*(x_2); 

  // paramters initialization
  char  matdescra[6];
    
  transa = 'n';           // not-transpose
  matdescra[1] = 's';    // Symmetric Matrix
  matdescra[3] = 'c';    // zero-based indexing
  double  alpha = 1.0, beta = 0.0;

  double *x_3  = new double[num_cols]; 
  mkl_dcsrmv(&transa, num_rows, num_cols, &alpha, matdescra, acsr, ja, ia, ia+1, &x_2, &beta, x_3);

  // Computation : q(x_3) = A*g; where g = (x_2) Complete;

  

  return 0;

}