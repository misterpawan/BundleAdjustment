#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstring>
#include <ctime>
#include <cmath>

using namespace std;

#include "cs.h"
#include "umfpack.h"
#include "util.h"


int main()
{

	void* Symbolic;
	void* Numeric;
	double* null = ( double* ) NULL;
	int sym_status=10,num_status = 10,solve_status = 10;

	cs* A = new cs;
	double* b;
	double* x1;
	//string prefix = "test/test_D";
	string prefix = "49/1/JTJ49_1";
	int ncc, num_cols;
	int j;

	//string rhs_filename = "test/test_D_b.txt";
	string rhs_filename = "49/1/JTe49_1.txt";

	// getting the matrix size: n -> no of cols, ncc -> no of non zeros
	cc_header_read ( prefix, ncc, num_cols );
	num_cols = num_cols+1;                          //DOUBT!!!!!!!!
	cout << "\nNo of non zeros = "<< ncc << "\n";
	cout << "\nNo of columns = "<< num_cols << "\n";
	//num_cols = num_cols + 1;

	b = new double[num_cols];

	r8vec_data_read ( rhs_filename, num_cols, b);

	A->nzmax = ncc;
	A->nz = -1;
	A->p = new int[num_cols+1];
	A->i = new int[ncc];
	A->x = new double[ncc];
	A->m = num_cols;
	A->n = num_cols;

	cc_data_read ( prefix, ncc, num_cols, A->i, A->p, A->x );
	A->p[num_cols] = ncc;

	//  From the matrix data, create the symbolic factorization information.
	sym_status = umfpack_di_symbolic ( A->m, A->n, A->p, A->i, A->x, &Symbolic, null, null );
	cout << "\n Symbolic status :" << sym_status << "\n";

	//  From the symbolic factorization information, carry out the numeric factorization.
  	num_status = umfpack_di_numeric ( A->p, A->i, A->x, Symbolic, &Numeric, null, null );
  	cout << "\n Numeric status :" << num_status << "\n";

  	//  Free the symbolic factorization memory.
  	umfpack_di_free_symbolic ( &Symbolic );

  	//  Using the numeric factorization, solve the linear system.
  	x1 = new double[num_cols];
  	solve_status = umfpack_di_solve ( UMFPACK_A, A->p, A->i, A->x, x1, b, Numeric, null, null );
  	cout << "\n Solve status :" << solve_status << "\n";

  	//  Free the numeric factorization.
  	umfpack_di_free_numeric ( &Numeric );

  /*	for(j=0;j<num_cols;j++)
  		cout << "\n"<<x1[j] << "\n";
*/

	delete [] b; delete [] A->p; delete A->i; delete [] A->x;
	delete A;
	return 0;
}