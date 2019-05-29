/* Mini Schur Complement implementation using the CHOLMOD Library
*/
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstring>
#include <ctime>
#include <cmath>

using namespace std;

#include "cholmod.h" 
#include "util.h"
#include "umfpack.h" 
#include "SuiteSparseQR.hpp"

#define TRUE 1
#define FALSE 0
#define sizeG 441 //hardcoding it for the time being


/* This function densifies the jth column of input matrix A
*/
//cholmod_dense* densify_column(cholmod_sparse* A, int j,cholmod_common* Common)
double* densify_column(cholmod_sparse* A, int j,cholmod_common* Common)
{
	//cholmod_dense* b = cholmod_l_allocate_dense(A->nrow,1,A->nrow,CHOLMOD_REAL,Common);
	double* b = new double[A->nrow];
	int k;
	int l = ((int*)A->p)[j]; //start index for rows of column j

	//b->nzmax = b->nrow;
	for(k = 0; k < A->nrow && l < ((int*)A->p)[j+1]; k++ )
	{
		if(k == ((int*)A->i)[l])
		{
			//((double*)b->x)[k] = ((double*)A->x)[l];
			b[k] = ((double*)A->x)[l];
			l += 1;
		}
		//else ((double*)b->x)[k] = 0;
		else b[k] = 0;
	}

	return b;
}


/* This function does the computation PD\PU (MATLAB '\') by solving 
a linear system with multiple right hand sides. The output is written 
in the AA structure.
*/
void compute_PDinv_times_PU(cholmod_sparse* PD,cholmod_sparse* PU,cholmod_sparse* AA,cholmod_common* Common)
{
	void* Symbolic;
	void* Numeric;
	double* null = ( double* ) NULL;
	int j,k;
	double* b;
	double* x1;
	int sym_status=10,num_status=10,solve_status=10;
	int nnz_b = 0;
	int index;
	int free_status;
	int factor;

	for(j=0;j<PU->nrow;j++)
	//for(j=0;j<1;j++)
	{
		//cholmod_factor* F; //  for storing factorization information
		//cholmod_dense* b;// = cholmod_l_allocate_dense(PU->nrow,1,PU->nrow,CHOLMOD_REAL,Common);
		//cholmod_dense* x1;// = cholmod_l_allocate_dense(PD->ncol,1,PD->ncol,CHOLMOD_REAL,Common);;
		//cholmod_sparse* b; //vector to store the rhs at each iteration
		//cholmod_sparse* x1;//vector to store the solution at each iteration
		double* b;
		double* x1 = new double[PD->ncol];

		// Copying the jth column of PU to b(rhs)
		b = densify_column(PU,j,Common);
		/*if(j < PU->ncol-1)
			nnz_b = ((int*)PU->p)[j+1] - ((int*)PU->p)[j]; 
		else
			nnz_b = PU->nzmax - ((int*)PU->p)[j]; 
		b = cholmod_l_allocate_sparse(PU->nrow,1,nnz_b,TRUE,FALSE,0,CHOLMOD_REAL,Common);
		
		((int*)b->p)[j] = ((int*)PU->p)[j];
		for(k = 0; k<nnz_b; k++)
		{	
			index = ((int*)PU->p)[j]+k;
			((int*)b->i)[k] = ((int*)PU->i)[index];
			((double*)b->x)[k] = ((double*)PU->x)[index];
		}
		*/

		//Allocate factor
		//F = cholmod_l_allocate_factor(PD->nrow,Common);  cout << "\n PD-> stype : "<< PD->stype << "\n";

		//  From the matrix data, create the symbolic factorization information.
		//F = cholmod_l_analyze(PD,Common);
		sym_status = umfpack_di_symbolic ( PD->nrow, PD->ncol, (int*)PD->p, (int*)PD->i, (double*)PD->x, &Symbolic, null, null );
		//cout << "\n Symbolic status :" << sym_status << "\n";

		//  From the symbolic factorization information, carry out the numeric factorization.
		//factor = cholmod_l_factorize(PD,F,Common);  cout << "\n Factor: "<< factor << "\n";
  		num_status = umfpack_di_numeric ( (int*)PD->p, (int*)PD->i, (double*)PD->x, Symbolic, &Numeric, null, null );
  		//cout << "\n Numeric status :" << num_status << "\n";

  		//  Free the symbolic factorization memory.
  		umfpack_di_free_symbolic ( &Symbolic );

  		//  Using the numeric factorization, solve the linear system.
		//x1 = SuiteSparseQR<double>(PD,b,Common);

  		//x1 = cholmod_l_spsolve(CHOLMOD_A,F,b,Common);
  		solve_status = umfpack_di_solve ( UMFPACK_A, (int*)PD->p, (int*)PD->i, (double*)PD->x, x1, b, Numeric, null, null );
  		//cout << "\n Solve status :" << solve_status << "\n";
		//  Free the numeric factorization.
  		umfpack_di_free_numeric ( &Numeric );

  		//for(k = 0; k < PU->m)

		//free_status = cholmod_l_free_dense(&b,Common); //cout << "\nfree status : " << free_status <<"\n";
		//free_status = cholmod_l_free_factor(&F,Common);
		//free_status = cholmod_l_free_dense(&x1,Common);
		delete [] b; delete [] x1;
	}
	cout << "\n inv(PD) PU dome ! \n";
	return;
}


/* This function computes the full Schur complement using the blocks 
given as input. The output id written in MSC.
*/
void compute_full_schur_complement(int r1,int r2,int rD1,int rD2,cholmod_sparse* MSC,cholmod_sparse* PD,cholmod_sparse* PL,
									cholmod_sparse* PU,cholmod_sparse* PG,cholmod_common* Common)
{
	int free_status;
	cholmod_sparse* S;
	cholmod_sparse* AA;

	//Allocating sparse matrices
	S = cholmod_l_allocate_sparse(PG->nrow,PG->ncol,sizeG*sizeG,TRUE,FALSE,0,CHOLMOD_REAL,Common);
	AA = cholmod_l_allocate_sparse(PG->nrow,PG->ncol,sizeG*sizeG,TRUE,FALSE,0,CHOLMOD_REAL,Common);

	compute_PDinv_times_PU(PD,PU,AA,Common);

	//deallocate matrix
	free_status = cholmod_l_free_sparse(&S,Common);
	free_status = cholmod_l_free_sparse(&AA,Common);

	return;
}


/* This function computed the mini schur complement of the input matrix.
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
void compute_mini_schur_complement(cholmod_sparse* A,cholmod_sparse* MSC,cholmod_sparse* D,cholmod_sparse* L,
								   cholmod_sparse* U,cholmod_sparse* G,cholmod_common* Common)
{
	int free_status;
	int nmsc_block = 3 ;  //no of blocks for G
	int r = sizeG % nmsc_block;
	int sz = (sizeG - r)/nmsc_block ; //size of each nmsc block for G
	int r1 = 0,r2 = sz; //C++ convention   
	int nD = nmsc_block ; //no of blocks for D (same as that for G)
	int sizeD = D->nrow;
	int rD = sizeD % nD;
	int szD = (sizeD - rD)/nD; //size of each block for D
	int rD1 = 0,rD2 = szD; //C++ convention

	//cout << "\n r1 :"<<r1 <<"\tr2 :"<<r2<<"\n";
	//cout << "\n rD1 :"<<rD1 <<"\trD2 :"<<rD2<<"\n";

	int ol = 0; //overlap size
	int i=0;
	int nzPD,nzPL,nzPU,nzPG;
	int j;
	int k,l;
	int iterPD,iterPL,iterPG;
	int status;
	int *null = ( int * ) NULL;

	cout<< "\nComputing Mini Schur Complement...\n";
	for(i ; i < nmsc_block -2;i++ )
	{
		//cout << "\n i : " << i << "\n";
		nzPD = 0; nzPL = 0; nzPU = 0; nzPG = 0;
		cholmod_sparse* PD; cholmod_sparse* PL; cholmod_sparse* PU ; cholmod_sparse* PG ;
		
		cout<< "\n Computing nnz..." << "\n";
		
		for(j = rD1;j<rD2+ol;j++)
		{
			for(k = ((int*)D->p)[j]; k < ((int*)D->p)[j+1];k++)
			{
				if(((int*)D->i)[k] >= rD1 && ((int*)D->i)[k] < rD2+ol) ++nzPD;
				else continue;
			}
			for(l = ((int*)L->p)[j];l<((int*)L->p)[j+1];l++)
			{
				if(((int*)L->i)[l] >= r1 && ((int*)L->i)[l] < r2+ol) ++nzPL;
				else continue;
			}
		}

		for(j = r1;j<r2+ol;j++)
		{
			for(k = ((int*)U->p)[j];k<((int*)U->p)[j+1];k++)
			{
				if(((int*)U->i)[k] >= rD1 && ((int*)U->i)[k] < rD2+ol) ++nzPU;
				else continue;
			}
			for(l = ((int*)G->p)[j];l< ((int*)G->p)[j+1];l++)
			{
				//cout << "\n G->i[j] "<< G->i[j] << "\n"; 
				if(((int*)G->i)[l] >= r1 && ((int*)G->i)[l] < r2+ol) ++nzPG;
				else continue;
			}
		}
		cout << "\n nzPD : "<< nzPD<<"\tnzPL :"<<nzPL<<"\tnzPU :"<<nzPU<<"\tnzPG :"<<nzPG<<"\n";

		//Allocating sparse matrices
		PD = cholmod_l_allocate_sparse((rD2+ol)-rD1,(rD2+ol)-rD1,nzPD,TRUE,FALSE,1,CHOLMOD_REAL,Common);
		PL = cholmod_l_allocate_sparse((r2+ol)-r1,(rD2+ol)-rD1,nzPL,TRUE,FALSE,0,CHOLMOD_REAL,Common);
		PG = cholmod_l_allocate_sparse((r2+ol)-r1,(r2+ol)-r1,nzPG,TRUE,FALSE,1,CHOLMOD_REAL,Common);
		PU = cholmod_l_allocate_sparse((rD2+ol)-rD1,(r2+ol)-r1,nzPL,TRUE,FALSE,0,CHOLMOD_REAL,Common);

	
		cout << "\n Allocating PD and PL..."<<"\n";
		iterPD =0; iterPL = 0;
		//filling the blocks PD and PL
		for(j=rD1;j<rD2+ol;j++)
		{
			//cout << "\n j : "<< j<<"\n";
			((int*)PD->p)[j-rD1] = iterPD;
			for(k = ((int*)D->p)[j];k<((int*)D->p)[j+1];k++)
			{
				if(((int*)D->i)[k] >= rD1 && ((int*)D->i)[k] < rD2+ol) 
				{
					((int*)PD->i)[iterPD] = ((int*)D->i)[k]-rD1;
					((double*)PD->x)[iterPD] = ((double*)D->i)[k];
					iterPD += 1;
				}
				else continue;
			}
			
			((int*)PL->p)[j-rD1] = iterPL;
			for(l = ((int*)L->p)[j];l<((int*)L->p)[j+1];l++)
			{
				if(((int*)L->i)[l] >= r1 && ((int*)L->i)[l] < r2+ol)
				{
					((int*)PL->i)[iterPL] = ((int*)L->i)[l] -rD1;
					((double*)PL->x)[iterPL] = ((int*)L->x)[l];
					iterPL += 1;
				}
				else continue;
			}
		}
		((int*)PD->p)[PD->ncol] = iterPD; ((int*)PL->p)[PL->ncol] = iterPL;
		//cout << "\n PL->p[ncol] :" << PL->p[PL->n] << "\n";
		//cout << "\n Diff PD : "<< nzPD - iterPD << "\n";
		//cout << "\n Diff PL : "<< nzPL - iterPL << "\n";

		status = umfpack_di_transpose(PL->nrow,PL->ncol, (int*)PL->p,(int*)PL->i,(double*)PL->x,null,null,
																							(int*)PU->p,(int*)PU->i,(double*)PU->x) ;
		//cout << "\n Transpose status : "<< status << "\n";

		iterPG = 0;
		for(j = r1;j<r2+ol;j++)
		{
			for(l = ((int*)G->p)[j];l< ((int*)G->p)[j+1];l++)
			{
				//cout << "\n G->i[j] "<< G->i[j] << "\n"; 
				((int*)G->p)[j-r1] = iterPG;
				if(((int*)G->i)[l] >= r1 && ((int*)G->i)[l] < r2+ol) 
				{
					((int*)PG->i)[iterPG] = ((int*)G->i)[l] - r1;
					((double*)PG->x)[iterPG] = ((double*)G->x)[l];
					iterPG += 1;
				}
				else continue;
			}
		}
		((int*)PG->p)[PG->ncol] = iterPG;
		//cout << "\n Diff PG : "<< nzPG - iterPG << "\n";

		//compute the full schur complement of the extracted blocks
		compute_full_schur_complement(r1,r2,rD1,rD2,MSC,PD,PL,PU,PG,Common);
		
		//deallocate matrix
	    free_status = cholmod_l_free_sparse(&PD,Common);
	    free_status = cholmod_l_free_sparse(&PU,Common);
	    free_status = cholmod_l_free_sparse(&PL,Common);
	    free_status = cholmod_l_free_sparse(&PG,Common);
	}


	

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
	int status,free_status,transpose_status;
	void *Symbolic;
	int sizeD;
	int j,k;
	int nzD = 0;
	int nzG = 0;
	int nzL = 0; // nzL = nzU
	int iterD = 0;
	int iterL = 0;
	int iterG = 0;

	//creating CHOLMOD structures
	cholmod_common Common;
	cholmod_sparse* A ;
	cholmod_sparse* MSC ;
	cholmod_sparse* D ;
	cholmod_sparse* L ;
	cholmod_sparse* G ;
	cholmod_sparse* U ;



	string line;
	string rhs_filename = "49/1/JTe49_1.txt";

	string col_filename = prefix + "_col.txt";
	//ifstream in_file(col_filename,ios::in);

	// getting the matrix size: n -> no of cols, ncc -> no of non zeros
	cc_header_read ( prefix, ncc, num_cols );
	num_cols = num_cols+1;                          //DOUBT!!!!!!!!
	cout << "\nNo of non zeros = "<< ncc << "\n";
	cout << "\nNo of columns = "<< num_cols << "\n";

	//size of the D matrix
	sizeD = num_cols - sizeG; cout << "\n sizeD : " << sizeD << "\n";

	//start CHOLMOD
	cholmod_l_start(&Common);

	//allocating the sparse matrices ....STYPE CHANGE 0 OR 1!!!!!!!!!!!!!!!!!!
	//params : (rows,cols,nnz,sorted_col,packed,symm(upper),xtype, cholmod_common*)
	A = cholmod_l_allocate_sparse(num_cols,num_cols,ncc,TRUE,FALSE,0,CHOLMOD_REAL,&Common);
	//A->p = (int*)A->p; A->i = (int*)(A->i); A->x = (double*)(A->x); 

	//cout << "\nA -> nzmax : " << A->nzmax <<"\n";

	//read the matrix data
	cc_data_read ( prefix, ncc, num_cols, (int*)A->i, (int*)A->p, (double*)A->x );
	cout << "\nFile read complete!\n";

	//cout << "\n A = " << ((int*)(A->p))[5]<< "\t"<<((int*)(A->i))[5] << "\t"<<((double*)(A->x))[5]<<"\n";

	

	cout << "\nCounting non zeros ...\n";
	for(j=0;j<sizeD;j++)
	{
		//j = 1000;
		
		for(k=((int*)A->p)[j]; k < ((int*)A->p)[j+1]; k++)
		{
			if(((int*)A->i)[k] < sizeD) ++nzD;
			else ++nzL;
			//cout << "\nk :"<<k<<"\n";
		}
	}
	//cout << "\n nzD = " << nzD << "\n";
	//cout << "\n nzL = " << nzL << "\n";

	nzG = ncc - (nzD + 2*nzL); 
	cout << "\nCounting non zeros complete!!\n";

	//Allocating sparse matrices
	MSC = cholmod_l_allocate_sparse(sizeG,sizeG,sizeG*sizeG,TRUE,FALSE,0,CHOLMOD_REAL,&Common);
	D = cholmod_l_allocate_sparse(sizeD,sizeD,nzD,TRUE,FALSE,1,CHOLMOD_REAL,&Common);
	L = cholmod_l_allocate_sparse(sizeG,sizeD,nzL,TRUE,FALSE,0,CHOLMOD_REAL,&Common);
	G = cholmod_l_allocate_sparse(sizeG,sizeG,nzG,TRUE,FALSE,1,CHOLMOD_REAL,&Common);
	U = cholmod_l_allocate_sparse(sizeD,sizeG,nzL,TRUE,FALSE,0,CHOLMOD_REAL,&Common);


	cout << "\nFilling non zeros ...\n";
	for(j=0;j<sizeD;j++)
	{
		//j = 1000;
		//nzD_col = 0;
		((int*)D->p)[j] = iterD;
		((int*)L->p)[j] = iterL;
		for(k = ((int*)A->p)[j]; k < ((int*)A->p)[j+1]; k++)
		{
			if(((int*)A->i)[k] < sizeD) 
			{
				((int*)D->i)[iterD] = ((int*)A->i)[k];
				((double*)D->x)[iterD] = ((double*)A->x)[k];
				iterD += 1;
				//++nzD_col;
			}
			else
			{
				((int*)L->i)[iterL] = ((int*)A->i)[k] - sizeD;
				((double*)L->x)[iterL] = ((double*)A->x)[k];
				iterL += 1;
			}
			//cout << "\nk :"<<k<<"\n";
		}
		
		//cout << "\n nzD = " << nzD << "\n";
		//cout << "\n nzL = " << nzL << "\n";
		//break;
	}
	((int*)D->p)[sizeD] = iterD;
	((int*)L->p)[sizeD] = iterL;

/*
	SuiteSparse_long fset[num_cols-1];

	for(k = 0; k < num_cols-1; k++)
	{
		fset[k] = ((int*)L->p)[k];
	}

	transpose_status = cholmod_l_transpose_unsym(L,1,NULL,fset,num_cols,U,&Common);
*/
	status = umfpack_di_transpose(L->nrow,L->ncol, (int*)L->p,(int*)L->i,(double*)L->x,null,null,(int*)U->p,(int*)U->i,(double*)U->x) ;

	//cout << "\n nzL : " << L->nzmax << "\tnzU : " << U->nzmax << "\n";

	for(j = sizeD;j<num_cols;j++)
	{
		((int*)G->p)[j - sizeD] = iterG;
		if(j < num_cols-1)
		{
			//cout << "\n In first condition...\n";
			for(k = ((int*)A->p)[j]; k < ((int*)A->p)[j+1];k++)
			{
				if(((int*)A->i)[k] < sizeD) continue;
				else
				{
					//cout << "\n row : " << A->i[k] << "\n";
					((int*)G->i)[iterG] = ((int*)A->i)[k]-sizeD;
					((double*)G->x)[iterG] = ((double*)A->x)[k];
					iterG += 1;
				}			
			}	
		}
		//for the last column
		else 
		{
			//cout << "\n In second condition..."<<A->nzmax <<"\n";
			for(k = ((int*)A->p)[j]; k < A->nzmax;k++)
			{
				if(((int*)A->i)[k] < sizeD) continue;
				else
				{
					//cout << "\n row : " << A->i[k] << "\n";
					((int*)G->i)[iterG] = ((int*)A->i)[k] - sizeD;
					((double*)G->x)[iterG] = ((double*)A->x)[k];
					iterG += 1;
				}
				//cout << "\n k : " << k << "\n";
				//break;
			}
		}
	    //break;
	}

	((int*)G->p)[sizeG] = iterG;
	//cout << "\n Diff G : "<<iterG - G->nzmax<<"\n";
	cout << "\nFilling non zeros complete!!\n";

	//compute the mini schur complement in the CCS format.
	compute_mini_schur_complement(A,MSC,D,L,U,G,&Common);


    //deallocate matrix
    free_status = cholmod_l_free_sparse(&A,&Common);
    free_status = cholmod_l_free_sparse(&MSC,&Common);
    free_status = cholmod_l_free_sparse(&D,&Common);
    free_status = cholmod_l_free_sparse(&L,&Common);
    free_status = cholmod_l_free_sparse(&G,&Common);
    free_status = cholmod_l_free_sparse(&U,&Common);

	//finish CHOLMOD
	cholmod_l_finish(&Common);

	return 0;
}