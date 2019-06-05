/* Two=Grid Solve as a preconditioner for GMRES */

// C++ libraries
#include <cstdlib>
#include <iomanip>         
#include <iostream>
#include <fstream>
#include <cstring>
#include <ctime>
#include <cmath>            
#include <algorithm>

#include "cs.h"                  // header for matrix structure
#include "umfpack.h"             // header for UMFPACK library
#include "util.h"                // header for matrix read from file
#include "mkl.h"                 // header for intel mkl genneral routines 
#include "mkl_solvers_ee.h"      // header for the Eigenvalue solve routine
#include "mkl_types.h"           // header for mkl_cspblas_scsrtrsv
#include "mkl_spblas.h"          // header for calling blas routines 

using namespace std;             // C++ standard libary namespace
#define max(a, b) (a) < (b) ? (b): (a)    // for eigenvalue solve routine
#define n_eig (num_cols*eig_N)
int main()

{

/* reads matrix in CSC format, note that: transpose(CSC) = CSR */

  double* JTe; //rhs
  int i;
  int num_rows;
  int num_cols; // no of columns
  int ncc = 0; // no of non zeros
  int *null = ( int * ) NULL;
  double *solve_null = ( double * ) NULL;
  void *Numeric;
  string prefix = "JTJ49_1.txt";
  double r;
  int status,sym_status,num_status,solve_status;
  void* Symbolic;
  int sizeD;
  int j,k;
  int nzD = 0;
  int nzG = 0;
  int nzL = 0; // nzL = nzU
  int iterD = 0;
  int iterL = 0;
  int iterG = 0;

  //creating structures of cs_diPARSE
  cs_di* A = new cs_di;
  cs_di* MSC = new cs_di;
  cs_di* D = new cs_di;
  cs_di* L = new cs_di;
  cs_di* G = new cs_di;
  cs_di* U = new cs_di;

  //string line;
  string rhs_filename = "JTe49_1.txt";


  // getting the matrix size: n -> no of cols, ncc -> no of non zeros
  cc_header_read ( prefix, ncc, num_cols );
  num_cols = num_cols+1;                          //DOUBT!!!!!!!!
  cout << "\nNo of non zeros = "<< ncc << "\n";
  cout << "\nNo of columns = "<< num_cols << "\n";



  A->nzmax = ncc;
  A->m = num_cols;
  A->n = num_cols; // since square matrix
  A->nz = -1;

  //Allocate space for rhs
  JTe = new double[num_cols];

  // reading the rhs
  r8vec_data_read ( rhs_filename, num_cols, JTe);

  cout << "\nRHS file read complete!!\n";

  //Allocate space for the coefficient matrix
  A->x = new double[ncc];
  A->p = new int[num_cols+1];
  A->i = new int[ncc];

  //read the matrix data
  cc_data_read ( prefix, ncc, num_cols, A->i, A->p, A->x );

  cout << "\nFile read complete!\n";

/*

The above read JTJ matrix is in CSC to convert to CSR: 
- Pass row pointers in place of column & column in place of rows 
- And then pass them to the routines.
This is true because for symmetric matrix JTJ, CSR = transpose(CSC)

----------------------------------------------------------------------------------

 Call the Eigenvalue/Vector Solver for Intel MKL*/
  /*
!   Content: Example for Intel(R) MKL Extended Eigensolvers (sparse format,
!            double precision)
!
!*******************************************************************************
!
! The following routines are used in the example:
!          DGEMM  DFEAST_SCSREV  DFEAST_SCSRGV  FEASTINIT
!
! Consider the matrix A
!                 |  5   2   1   1   0   0   0   0   0   0   0  |
!                 |  2   6   3   1   1   0   0   0   0   0   0  |
!                 |  1   3   6   3   1   1   0   0   0   0   0  |
!                 |  1   1   3   6   3   1   1   0   0   0   0  |
!                 |  0   1   1   3   6   3   1   1   0   0   0  |
!    A    =       |  0   0   1   1   3   6   3   1   1   0   0  |,
!                 |  0   0   0   1   1   3   6   3   1   1   0  |
!                 |  0   0   0   0   1   1   3   6   3   1   1  |
!                 |  0   0   0   0   0   1   1   3   6   3   1  |
!                 |  0   0   0   0   0   0   1   1   3   6   2  |
!                 |  0   0   0   0   0   0   0   1   1   2   5  |
!
! stored as sparse matrix.
!
!  In what follows the symbol ' represents a transpose operation.
!
!  The test performs the following operations:
!
!  Step 1. Calls FEASTINIT to define the default values for the input
!          FEAST parameters.
!
!  Step  2. The code solves the standard eigenvalue problem Ax=ex using
!          DFEAST_SCSREV.
!
! 
!*******************************************************************************/


  char UPLO = 'F'; /* Type of matrix: (F means full matrix, L/U - lower/upper triangular part of matrix) */
    /* Matrix A of size N in CSR format. Size N and 3 arrays are used to store matrix in CSR format */
  const MKL_INT N = num_cols;
  const MKL_INT eig_N = 5;     // Hard-coding it change later 

/*MKL_INT  rows[12] = { 1, 5, 10, 16, 23, 30, 37, 44, 51, 57, 62, 66 };
  MKL_INT  cols[65] = {    1,   2,   3,   4,
                                  1,   2,   3,   4,   5,
                                  1,   2,   3,   4,   5,   6,
                                  1,   2,   3,   4,   5,   6,   7,
                                       2,   3,   4,   5,   6,   7,   8,
                                            3,   4,   5,   6,   7,   8,   9,
                                                 4,   5,   6,   7,   8,   9,  10,
                                                      5,   6,   7,   8,   9,  10,  11,
                                                           6,   7,   8,   9,  10,  11,
                                                                7,   8,   9,  10,  11,
                                                                     8,   9,  10,  11
                            };

  double  val[65] = {   5.0, 2.0, 1.0, 1.0,
                                2.0, 6.0, 3.0, 1.0, 1.0,
                                1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                     1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                          1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                               1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                                    1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                                         1.0, 1.0, 3.0, 6.0, 3.0, 1.0,
                                                              1.0, 1.0, 3.0, 6.0, 2.0,
                                                                   1.0, 1.0, 2.0, 5.0
                            };*/

                                
    /* Declaration of FEAST variables */
  MKL_INT      fpm[128];      /* Array to pass parameters to Intel(R) MKL Extended Eigensolvers */
  double       Emin, Emax;    /* Lower/upper bound of search interval [Emin,Emax] */

  double       epsout;        /* Relative error on the trace */
  MKL_INT      loop;          /* Number of refinement loop */
  MKL_INT      L = 5;
  MKL_INT      M0;            /* Initial guess for subspace dimension to be used */
  MKL_INT      M;             /* Total number of eigenvalues found in the interval */

  double       E[5];          /* # of Eigenvalues = 5*/
  double       X[n_eig];      /* Eigenvectors, n_eig defined in macro*/
  double       res[5];        /* Residual */
  double       Ac[25];        /* Coarse Grid Matrix        
  /* Declaration of local variables */
  MKL_INT      info;          /* Errors */
  double       Eig[5];        /* Eig - array for storing exact eigenvalues */
//double       R[11];         /* R = |E-Eig| */
    
  /* Exact eigenvalues in range (3.0, 7.0) */
    /* Needs to change for JTJ*/
    /*Eig[0] = 3.1715728752538100;
    Eig[1] = 4.0000000000000000;
    Eig[2] = 4.0000000000000000;
    Eig[3] = 4.1292484841890931;
    Eig[4] = 4.4066499006731521;
    Eig[5] = 6.0000000000000000;
    */
  for (i=0; i<eig_N; i++)
  {
      Eig[i] = 0.0;
  }

  printf("\n FEAST DFEAST_SCSREV AND DFEAST_SCSRGV EXAMPLE\n");
  /* Initialize matrix X */
  for (i=0; i<n_eig; i++)
  {
        X[i] = zero;
  }

  printf("Sparse matrix size %i\n", (int)num_cols);

  /* Search interval [Emin,Emax] */
  Emin = 1.0e+11;   // change this for JTJ
  Emax = 1.0e+9;
  printf("Search interval [ %.15e, %.15e  ]  \n", Emin, Emax);

  M0     = L;
  M      = L;
  loop   = 0;
  info   = 0;
  epsout = 0.0;

  /* Step 1. Call  FEASTINIT to define the default values for the input FEAST parameters */
  feastinit(
     fpm /* OUT: Array is used to pass parameters to Intel(R) MKL Extended Eigensolvers */
  );

  fpm[0] =  1; /* Extended Eigensolver routines print runtime status to the screen. */

    

  /* Step 2. Solve the standard Ax = ex eigenvalue problem. */

  printf("Testing dfeast_scsrev routine:\n");

  dfeast_scsrev(
      &UPLO,   /* IN: UPLO = 'F', stores the full matrix */
      &N,      /* IN: Size of the problem */
      A->x,     /* IN: CSR matrix A, values of non-zero elements */
      A->i,    /* IN: CSR matrix A, index of the first non-zero element in row */
      A->p,    /* IN: CSR matrix A, columns indices for each non-zero element */
      fpm,     /* IN/OUT: Array is used to pass parameters to Intel(R) MKL Extended Eigensolvers */
      &epsout, /* OUT: Relative error of on the trace */
      &loop,   /* OUT: Contains the number of refinement loop executed */
      &Emin,   /* IN: Lower bound of search interval */
      &Emax,   /* IN: Upper bound of search interval */
      &M0,     /* IN: The initial guess for subspace dimension to be used. */
      E,       /* OUT: The first M entries of Eigenvalues */
      X,       /* IN/OUT: The first M entries of Eigenvectors */
      &M,      /* OUT: The total number of eigenvalues found in the interval */
      res,     /* OUT: The first M components contain the relative residual vector */
      &info    /* OUT: Error code */
  );

  printf("FEAST OUTPUT INFO %d \n",info);   // info status of dfeasr_scsrev routine
  if ( info != 0 )
  {
      printf("Routine dfeast_scsrev returns code of ERROR: %i", (int)info);
      return 1;
  }

  /* Construction of Coarse -Grid Matrix
     Ac = diagonal matrix of eigenvalues
     hardcoding it for top 5 eigenvalues as of now, change later!
  */

  int idx = 0;
  int j = 0;
  for ( int i = 0; i < 25 ; i++)
  {
    if( i % 6 == 0){
     idx = 6*i;
     Ac[idx] = E[j];    
     j++;
    }
  } // end of for: Coarse Grid Matrix Constructed
  

  // Prolongation Matrix , P = X (Matrix of Eigenvectors)
  // coarse grid matrix , Ac = Diagonal matrix of eigenvalues


  // Construct Smoother: D = blkdiag(A(1:nJ1, 1:nJ1), A(nJ+1:n, nJ+1:n))

  /***Domain Decomposition***/

  sizeD = num_cols;   // take D and G in D only
  D->p = new int[sizeD+1];   // number of columns
  D->nz = -1;D->m = sizeD;D->n = sizeD;   // same for m and n as JTJ' = JTJ

  
  cout << "\nCounting non zeros ...\n";
  
  for(j=0;j<sizeD;j++)
  {
    //j = 1000;
    
    for(k=A->p[j]; k < A->p[j+1]; k++)
    {
      if(A->i[k] < sizeD) ++nzD;
      else ++nzL;
      //cout << "\nk :"<<k<<"\n";
    }
    //cout << "\n nzD = " << nzD << "\n";
    //cout << "\n nzL = " << nzL << "\n";
    //break;
  }
  cout << "\nCounting non zeros complete!!\n";
  //cout << "\n nzD = " << nzD << "\n";
  //cout << "\n nzL = " << nzL << "\n";

  //Allocating memory
  D->i = new int[nzD];
  D->x = new double[nzD];

  //Setting values 
  D->nzmax = nzD;

  cout << "\nFilling non zeros ...\n";
  for(j=0;j<sizeD;j++)
  {
    //j = 1000;
    //nzD_col = 0;
    D->p[j] = iterD;
    for(k=A->p[j]; k < A->p[j+1]; k++)
    {
      if(A->i[k] < sizeD) 
      {
        D->i[iterD] = A->i[k];
        D->x[iterD] = A->x[k];
        iterD += 1;
        //++nzD_col;
      }
      else
      {
        L->i[iterL] = A->i[k] - sizeD;
        L->x[iterL] = A->x[k];
        iterL += 1;
      }
      //cout << "\nk :"<<k<<"\n";
    }
    
    //cout << "\n nzD = " << nzD << "\n";
    //cout << "\n nzL = " << nzL << "\n";
    //break;
  }

  D->p[sizeD] = iterD;
  L->p[sizeD] = iterL;
  
  cout << "\nConstruction of Smoother complete!\n";
  



  /* Two-grid Solve Function */


  float * tg_solve(float *){

  /* 1st Step : Pre-smoothing  */
    
  printf("   Solve t = U\\(L\\x) system   \n");

  // Paramteres for triangular solver routine - mkl_cspblas_scsrtrsv
  
  float t_0[num_cols];
  float t_1[num_cols];  
  uplo = 'l';        // lower triangular
  nonunit = 'n';     // diagonal is not unit triangular
  transa = 'n';      // not tranpose

  /* Compute L\x */
  mkl_cspblas_scsrtrsv(&uplo, &transa, &nonunit, num_cols, D->x, D>p, D->i, x, t_0);
  /* Computation of L\X complete: L\x = t_0  */

  uplo = 'u';        // upper triangular
  nonunit = 'n';     // diagonal is not unit triangular
  transa = 'n';      // not tranpose
  /* Compute U\(t_0); t_0 = L\x */
  mkl_cspblas_scsrtrsv(&uplo, &transa, &nonunit, num_cols, D->x, D>p, D->i, t_0, t_1); 
  /*  Solution obtained in t_1  */

  printf(" Pre-smoothing Complete");


        
  /* 2nd Step: Coarse Grid Solve:  g = P* (Ac \ (P'*x))  

   1st perform the innermost matrix vector multiplication i.e P'*x  

  */

    printf("  INPUT DATA FOR MKL_DCSRMV \n");
         
    float x_0[num_cols];
    float x_1[num_cols];
    char  matdescra[6];
    
    transa = 'T';           // Transpose of the prolongation matrix
    matdescra[1] = 's';    // Symmetric Matrix
    matdescra[3] = 'c';    // zero-based indexing
    float  alpha = 1.0, beta = 0.0;

    for (i = 0; i < num_rows; i++) {
       printf("%7.1f\n", x_0[i]);
    }
    // Constructing pointerE
    MKL_INT  pointerE[M];
    MKL_INT  i, j;   
    for (i = 0; i < num_cols; i++) {
        pointerE[i] = P->i[i+1];
    }
        
    /*
      Parameter initiallization for sparse matrix vector MKL multiplication complete
      Now call sparse matrix vector multiplication routines
    */

    //mkl_dcsrmv(&transa, &m, &m, &alpha, matdescra, values, columns, pointerB, pointerE, sol_vec, &beta, rhs_vec);
    mkl_dcsrmv(&transa, num_rows, num_cols, &alpha, matdescra, P->x, P->p, P->i, pointerE, &x, &beta, x_0);

    /*
      g = P * (Ac \ (P'*x))
      Now : x_0 = P'*x  (done)
      Next : solve (Ac \x)
    */

    printf("   Solve Ac\\ x_0 system   \n");

    uplo = 'u';
    nonunit = 'n';
    transa = 'n'; 

    mkl_cspblas_scsrtrsv(&uplo, &transa, &nonunit, num_cols, Ac->x, Ac->p, Ac->i, x_0, x_1);

    /*
      
      Now compute matrix -vector multiplication :g = P * x_1 ; x1 = (Ac \ (P'*x_0));
    
    */
    transa = 'n';
        
    // x_2 (g) = P*x_1;   

    mkl_dcsrmv(&transa, num_rows, num_cols, &alpha, matdescra, P->x, P->p, P->i, pointerE, &x_1, &beta, x_2);

    // Compute : q(x_3) = A*(x_2); 

    mkl_dcsrmv(&transa, num_rows, num_cols, &alpha, matdescra, A->x, A->p, A->i, pointerE, &x_1, &beta, x_3);

    /*
        y = t + g - U \ (L\q)

        y = t + x_1 - U\(L\x_2)

        
    */

    // Initializing parameters for sparse vector-vector addition 

    double   alpha = 1;
    int incx = 1, incy = 1;

    // t + g; where g = x_2; t= t_1

    cblas_daxpy(n, alpha, t_1, incx, x_2, incy);
    // addition of t+g complete, x_2 = x_2 + t_1

    
    //mkl_cspblas_scsrtrsv(&uplo, &transa, &nonunit, num_cols, D->x, D->p, D->i, x, t_0);
    
    float q_0[num_cols];
    float q_1[num_cols]; 
    
    uplo = 'l';             // lower-triangular
    nonunit = 'n';          // not unit triangular
    transa = 'n';           // not transpose
    /* Compute (L\q)
       q = x_3
    */
    mkl_cspblas_scsrtrsv(&uplo, &transa, &nonunit, num_cols, D->x, D>p, D->i, x_3, q_0);
    /* Computation of L\q complete: L\x = t_0  */
    uplo = 'u';             // upper - triangular
    nonunit = 'n';          // not unit tringular
    transa = 'n';           // not transpose
   /* Compute U\(t_0); t_0 = L\x */
    mkl_cspblas_scsrtrsv(&uplo, &transa, &nonunit, num_cols, D->x, D>p, D->i, q_0, q_1); 
    /*  Solution obtained in t_1  */

    /* Initializing parameters for sparse vector-vector subtraction 
       x_2 - q_1 
       where x_2 = t+g
       and q_1 = U\(L\q)

    */
    double   alpha = -1;      // for subtratcion
    int incx = 1, incy = 1;
    cblas_daxpy(n, alpha, q_1, incx, x_2, incy);      // x_2 - q_1 written in form -q_1 + x_2
    // Subtraction Complete: x_2 = x_2 - q_1

    /*
      Two - Grid Solve Complete 
      Solution obtained in x_2

      Now pass this preconditioned solution to GMRES

    */
    return x_2;  // Solution of preconditioner - tg_solve

 }

/* Free Memory */
  delete [] JTe; 
  delete [] D->p;
  delete [] D->i; 
  delete [] D->x; 
  delete A; 

  return 0;
}







