/* Two=Grid Solve as a preconditioner for GMRES */
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstring>
#include <ctime>
#include <cmath>
#include <algorithm>

using namespace std;

#include "cs.h"
#include "umfpack.h"
#include "util.h"
#include "sort.h"
#include "mkl.h"
#include "mkl_solvers_ee.h"      // Contains the Eigenvalue solve routine

int main()

{

/* reads matrix market format (coordinate) and returns
   csr format */
  double* JTe; //rhs
  int i;
  int num_rows;
  int num_cols; // no of columns
  int ncc = 0; // no of non zeros
  int *null = ( int * ) NULL;
  double *solve_null = ( double * ) NULL;
  void *Numeric;
  string prefix = "JTJ49_1";
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




/* Call the Eigenvalue/Vector Solver for Intel MKL*/
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

#define max(a, b) (a) < (b) ? (b): (a)

int main()
{
    char          UPLO = 'F'; /* Type of matrix: (F means full matrix, L/U - lower/upper triangular part of matrix) */
    /* Matrix A of size N in CSR format. Size N and 3 arrays are used to store matrix in CSR format */
    const MKL_INT N = num_cols;
    const MKL_INT eig_N = 5;     // Hard-coding it change later 

   // MKL_INT       rows[12] = { 1, 5, 10, 16, 23, 30, 37, 44, 51, 57, 62, 66 };
    /*MKL_INT       cols[65] = {    1,   2,   3,   4,
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

    double        val[65] = {   5.0, 2.0, 1.0, 1.0,
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

    double       E[5];         /* Eigenvalues */
    double       X[25];        /* Eigenvectors */
    double       res[5];       /* Residual */

    /* Declaration of local variables */
    MKL_INT      info;          /* Errors */
    double       Eig[5];       /* Eig - array for storing exact eigenvalues */
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
    for (i=0; i<eig_N*eig_N; i++)
    {
        X[i] = zero;
    }

    printf("Sparse matrix size %i\n", (int)N);

    /* Search interval [Emin,Emax] */
    Emin = 3.0;   // change this for JTJ
    Emax = 7.0;
    printf("Search interval [ %.15e, %.15e  ]  \n", Emin, Emax);

    M0   = L;
    M    = L;
    loop = 0;
    info = 0;
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
        A->p,    /* IN: CSR matrix A, index of the first non-zero element in row */
        A->i,    /* IN: CSR matrix A, columns indices for each non-zero element */
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
    printf("FEAST OUTPUT INFO %d \n",info);
    if ( info != 0 )
    {
        printf("Routine dfeast_scsrev returns code of ERROR: %i", (int)info);
        return 1;
    }

// Prolongation Matrix , P = X (matrix of Eigenvectors)



double * tg_solve(double *){

        void* Symbolic;
        void* Numeric;
        double* null = ( double* ) NULL;
        int j,k;
        double* x;
        double* t;
        int sym_status,num_status,solve_status;
                                                   
        
        // 1st Step : Pre-smoothing

        /* t = U\(L\U\x)   */

        //  From the matrix data, create the symbolic factorization information.
        sym_status = umfpack_di_symbolic ( A->n, A->n, A->p, A->i, A->x, &Symbolic, null, null );
        //cout << "\n Symbolic status :" << sym_status << "\n";

        //  From the symbolic factorization information, carry out the numeric factorization.
        num_status = umfpack_di_numeric ( A->p, A->i, A->x, Symbolic, &Numeric, null, null );
        //cout << "\n Numeric status :" << num_status << "\n";

        //  Using the numeric factorization, solve the linear system.
        solve_status = umfpack_di_solve ( UMFPACK_A, A->p, A->i, A->x, t, x, Numeric, null, null );
      //cout << "\n Solve status :" << solve_status << "\n";

        /*  Solution obtained in t  */

        
        /* 2nd Step:   g = P* (Ac \ (P'*x)) */

        // Check for matrix vector multiplication 













        //  Free the symbolic factorization memory.
      umfpack_di_free_symbolic ( &Symbolic );
        //  Free the numeric factorization.
      umfpack_di_free_numeric ( &Numeric );
}



  return 0;
}







