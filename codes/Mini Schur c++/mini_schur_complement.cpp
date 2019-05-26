//this is a test code for mini schur complement

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

#define sizeG 441 //hardcoding it for the time being

/*This function computed the mini schur complement of the input matrix.
PARAMS : 
	n : INPUT, no of columns of the input matrix
	cc_row : INPUT,row indices of the input matrix
	cc_col : INPUT,pointers to start of the columns for the input matrix
	cc_val : INPUT,values of the input matrix 
	row_msc : OUTPUT,row indices of the output matrix
	col_msc : OUTPUT,pointers to start of the columns for the output matrix
	val_msc : OUTPUT,values of the output matrix

The domain decomposition step is not shown as it is assumed that the size of 
block G is available. 
*/
void compute_mini_schur_complement(int num_cols,cs* A,cs* MSC,cs* D,cs* L,cs* G)
{
	


	
	return;
}



int main()
{
	
	double* JTe; //rhs
	int i;
	int num_rows;
	int num_cols; // no of columns
	int ncc = 0; // no of non zeros
	int *null = ( int * ) NULL;
	void *Numeric;
	string prefix = "49/1/JTJ49_1";
	double r;
	int status;
	void *Symbolic;
	int sizeD;
	int j,k;
	int nzD = 0;
	int nzG = 0;
	int nzL = 0; // nzL = nzU
	

	//creating structures of CSPARSE
	cs* A = new cs;
	cs* MSC = new cs;
	cs* D = new cs;
	cs* L = new cs;
	cs* G = new cs;
	cs* U = new cs;

	string line;
	string rhs_filename = "49/1/JTe49_1.txt";

	string col_filename = prefix + "_col.txt";
	ifstream in_file(col_filename,ios::in);

	// getting the matrix size: n -> no of cols, ncc -> no of non zeros
	cc_header_read ( prefix, ncc, num_cols );
	num_cols = num_cols+1;                          //DOUBT!!!!!!!!
	cout << "\nNo of non zeros = "<< ncc << "\n";
	cout << "\nNo of columns = "<< num_cols << "\n";

	//size of the D matrix
	sizeD = num_cols - sizeG; 
	cout << "\nsizeD = "<< sizeD << "\n";

	A->nzmax = ncc;
	A->m = num_cols;
	A->n = num_cols; // since square matrix
	A->nz = -1;

	//Allocate space for rhs
	JTe = new double[num_cols+1];

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

	/***Domain Decomposition***/

	//Allocating memory for blocks
	MSC->p = new int[sizeG+1];
	D->p = new int[sizeD+1];
	L->p = new int[sizeD+1];
	G->p = new int[sizeG+1];	
	U->p = new int[sizeG+1];

	MSC->nz = -1;
	MSC->m = sizeG;
	MSC->n = sizeG;
	
	D->nz = -1;
	D->m = sizeD;
	D->n = sizeD;

	G->nz = -1;
	G->m = sizeG;
	G->n = sizeG;

	L->nz = -1;
	L->m = sizeG;
	L->n = sizeD;

	U->nz = -1;
	U->m = sizeD;
	U->n = sizeG;


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

	nzG = ncc - (nzD + 2*nzL);

	//Allocating memory
	D->i = new int[nzD];
	D->x = new double[nzD];
	L->i = new int[nzL];
	L->x = new double[nzL];
	U->i = new int[nzL];
	U->x = new double[nzL];
	G->i = new int[nzG];
	G->x = new double[nzG];

	//setting values
	D->nzmax = nzD;
	L->nzmax = nzL;
	U->nzmax = nzL;
	G->nzmax = nzG;

	int iterD = 0;
	int iterL = 0;
	cout << "\nFilling non zeros ...\n";
	for(j=0;j<sizeD;j++)
	{
		//j = 1000;
		//nzD_col = 0;
		D->p[j] = iterD;
		L->p[j] = iterL;
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
	cout << "\nFilling non zeros complete!!\n";

	//cout << "\n nzD = " << D->p[sizeD] << "\n";
	//cout << "\n nzL = " << L->p[sizeD] << "\n";
	

	//status = umfpack_di_transpose (n_row, n_col, Ap, Ai, Ax, P, Q, Rp, Ri, Rx) ;
	status = umfpack_di_transpose(sizeG,sizeD, L->p,L->i,L->x,null,null,U->p,U->i,U->x) ;

	cout << "\n TRANSPOSE STATUS : "<< status<< "\n";
	//cout << "\nD[9][9](591.616) = "<< L->x[0]<<"\n";
	//cout << "\nL[23714][9](335.400) = "<< U->x[0]<<"\n";

	
	
	//compute the mini schur complement in the CCS format.
	//compute_mini_schur_complement(num_cols,A,MSC,D,L,G);

	/****GMRES CALL AFTER THIS ****/


	delete [] JTe;
	delete [] MSC->p;
	delete [] D->p;
	delete [] L->p;
	delete [] G->p;
	delete [] U->p;
	delete [] D->i;
	delete [] L->i;
	delete [] U->i;
	delete [] G->i;
	delete [] D->x;
	delete [] L->x;
	delete [] U->x;
	delete [] G->x;
	delete U;
	delete A;
	delete MSC;
	delete D;
	delete L;
	delete G;

	return 0;
}