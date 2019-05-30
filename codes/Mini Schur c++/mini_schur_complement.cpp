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

/* This function densifies the jth column of input matrix A
*/
double* densify_column(cs_di* A, int j)
{

	double* b = new double[A->m];
	
	int k;
	int l = A->p[j]; //start index for rows of column j

	for(k = 0; k < A->m && l < A->p[j+1]; k++ )
	{
		if(k == A->i[l])
		{
			b[k] = A->x[l];
			l += 1;
		}
		else b[k] = 0;
		//cout << "\n In densify .. j,k : "<< j<< "\t" << k << "\n";
	}

	return b;
}



/* This function stores the dense array x1 into the jth column of AA as a sparse array
*/
void make_x_sparse(cs_di* AA,int j,double* x1)
{
	int k,l;
	int row_arr_count; //displacement on the row and val array where the next iterm must be stored.
	int nz_count = 0;

	if(j==0) AA->p[j] = 0;
	
	row_arr_count = AA->p[j];

	for(k=0;k<AA->m;k++)
	{
		if(x1[k] != 0)
		{
			//nz_count += 1;
			AA->i[row_arr_count+nz_count] = k;
			AA->x[row_arr_count+nz_count] = x1[k];  
			nz_count += 1;
		}
		else continue;
	}

	if(j < AA->n) AA->p[j+1] = row_arr_count+nz_count;

	//cout << "\nrow_arr_count+nz_count : " << row_arr_count+nz_count << "\n";

	return;
}


/* This function does the computation PD\PU (MATLAB '\') by solving 
a linear system with multiple right hand sides. The output is written 
in the AA structure.
*/
void compute_PDinv_times_PU(cs_di* PD,cs_di* PU,cs_di* AA)
{
	void* Symbolic;
	void* Numeric;
	double* null = ( double* ) NULL;
	int j,k;
	double* b;
	double* x1;
	int sym_status=10,num_status = 10,solve_status = 10;
	//int nz_count;

	for(j=0;j<PU->n;j++)
	//for(j=0;j<1;j++)
	{
		b = new double[PU->m]; //vector to store the densified rhs at each iteration
		x1 = new double[PD->n]; //vector to store the solution at each iteration

		b = densify_column(PU,j); // densifies the jth column of PU 

		//  From the matrix data, create the symbolic factorization information.
		sym_status = umfpack_di_symbolic ( PD->m, PD->n, PD->p, PD->i, PD->x, &Symbolic, null, null );
		//cout << "\n Symbolic status :" << sym_status << "\n";
		//  From the symbolic factorization information, carry out the numeric factorization.
  		num_status = umfpack_di_numeric ( PD->p, PD->i, PD->x, Symbolic, &Numeric, null, null );
  		//cout << "\n Numeric status :" << num_status << "\n";
  		//  Free the symbolic factorization memory.
  		umfpack_di_free_symbolic ( &Symbolic );
  		//  Using the numeric factorization, solve the linear system.
  		solve_status = umfpack_di_solve ( UMFPACK_A, PD->p, PD->i, PD->x, x1, b, Numeric, null, null );
  		//cout << "\n Solve status :" << solve_status << "\n";
		//  Free the numeric factorization.
  		umfpack_di_free_numeric ( &Numeric );
  		
		
  		//Store x1 in the jth column of AA
  		make_x_sparse(AA,j,x1);

		delete []  x1; delete [] b;
	}

	//cout << "\nAA->x[210] : " << AA->x[210] << "\n";

	return;
}



/* This function computes the full Schur complement using the blocks 
given as input. The output id written in MSC.
*/
void compute_full_schur_complement(int r1,int r2,int rD1,int rD2,int ol,cs_di* MSC,cs_di* PD,cs_di* PL,cs_di* PU,cs_di* PG,int* total_nz)
{
	
	cs_di* AA = new cs_di;
	AA->nz = -1;
	AA->m = PD->m; AA->n = PU->n; AA->nzmax = AA->m * AA->n;
	AA->p = new int[AA->n+1]; AA->i = new int[AA->nzmax]; AA->x = new double[AA->nzmax];

	cs_di* LDU;// = new cs_di;
	
	compute_PDinv_times_PU(PD,PU,AA);

	cout << "\nAA->p[AA->n] : " << AA->p[AA->n] << "\n";

	LDU = cs_di_multiply(PL,AA);
	
	delete [] AA->p; delete [] AA->i; delete AA->x;
	delete AA; //delete LDU;

	cs_di* S ;//= new cs_di;
	//S->nz = -1;
	//S->m = PG->m; S->n = PG->n; S->nzmax = S->m * S->n;
	//S->p = new int[S->n+1]; S->i = new int[S->nzmax]; S->x = new double[S->nzmax];

	
	S = cs_di_add(PG,LDU,1.0,-1.0); //subtraction
	
	//filling up the msc block
	int j,k,l;
	for(j=0;j< S->n; j++)
	{
		if(r1 == 0 && j == 0) MSC->p[0] = 0;
		else
			MSC->p[r1+j] = S->p[j]+MSC->p[r1+j-1];
	}	
	
	l = 0;
	for(j = r1; j < r2+ol; j++)
	{
		if(j < sizeG-1)
		{
			for(k = MSC->p[j];k  < MSC->p[j+1] && l < S->nzmax; k++)
			{
				MSC->i[k] = S->i[l] + r1;
				MSC->x[k] = S->x[l];
				l++;
			}
		}
		else
		{
			for(k = MSC->p[j];k  < MSC->nzmax && l < S->nzmax; k++)
			{
				MSC->i[k] = S->i[l] + r1;
				MSC->x[k] = S->x[l];
				l++;
			}
		}
	}
	*total_nz += S->nzmax;  cout << "\nTotal nnz : " << *total_nz << "\n";
	//cout << "\n Total nnz : "<< S->nzmax << "\n";
	//cout << "\n Last element : "<< S->x[S->nzmax-1] << "\n";

	delete [] LDU->p; delete [] LDU->i; delete LDU->x;
	delete [] S->p; delete [] S->i; delete S->x;
	//delete [] AA->p; delete [] AA->i; delete AA->x;
	//delete [] LD->p; delete [] LD->i; delete LD->x;
	delete S; delete LDU;//delete PSC;delete LD; 

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
void compute_mini_schur_complement(cs_di* A,cs_di* MSC,cs_di* D,cs_di* L,cs_di* U,cs_di* G)
{
	int nmsc_block = 3 ;  //no of blocks for G
	int r = sizeG % nmsc_block;
	int sz = (sizeG - r)/nmsc_block ; //size of each nmsc block for G
	int r1 = 0,r2 = sz; //C++ convention   
	int nD = nmsc_block ; //no of blocks for D (same as that for G)
	int sizeD = D->m;
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
	int total_nz = 0;
	int oll = 0; //overlap for 2nd last msc


	cout<< "\nComputing Mini Schur Complement...\n";
	for(i ; i < nmsc_block -2;i++ )
	{
		//cout << "\n i : " << i << "\n";
		nzPD = 0; nzPL = 0; nzPU = 0; nzPG = 0;
		cs_di* PD = new cs_di; cs_di* PL = new cs_di; cs_di* PU = new cs_di; cs_di* PG = new cs_di;
		PD->m = (rD2+ol) -rD1; PD->n = (rD2+ol) -rD1;
		PL->m = (r2+ol) -r1; PL->n = (rD2+ol) -rD1;  
		PU->m = (rD2+ol) -rD1; PU->n = (r2+ol) -r1;
		PG->m = (r2+ol) -r1; PG->n = (r2+ol) -r1;
		PD->p = new int[PD->n+1];PL->p = new int[PL->n+1];PU->p = new int[PU->n+1];PG->p = new int[PG->n+1];
		PD->nz = -1; PL->nz = -1; PU->nz = -1;PG->nz = -1;

		//cout<< "\n Computing nnz..." << "\n";
		//This is the wrong way as we only want the non zeros of the blocks and not the whole column
		//nzPD = D->p[rD2] - D->p[rD1]; nzPL = L->p[rD2] - L->p[rD1];
		//nzPU = U->p[r2] - U->p[r1];          // cout<< "U->p[r2] :"<< U->p[r2] <<"\t U->p[r1] :"<<U->p[r1] << "\n";
		//nzPG = G->p[r2] - G->p[r1];
		for(j = rD1;j<rD2+ol;j++)
		{
			for(k = D->p[j];k<D->p[j+1];k++)
			{
				if(D->i[k] >= rD1 && D->i[k] < rD2+ol) ++nzPD;
				else continue;
			}
			for(l = L->p[j];l<L->p[j+1];l++)
			{
				if(L->i[l] >= r1 && L->i[l] < r2+ol) ++nzPL;
				else continue;
			}
		}

		for(j = r1;j<r2+ol;j++)
		{
			for(k = U->p[j];k<U->p[j+1];k++)
			{
				if(U->i[k] >= rD1 && U->i[k] < rD2+ol) ++nzPU;
				else continue;
			}
			for(l = G->p[j];l< G->p[j+1];l++)
			{
				//cout << "\n G->i[j] "<< G->i[j] << "\n"; 
				if(G->i[l] >= r1 && G->i[l] < r2+ol) ++nzPG;
				else continue;
			}
		}
		//cout << "\n nzPD : "<< nzPD<<"\tnzPL :"<<nzPL<<"\tnzPU :"<<nzPU<<"\tnzPG :"<<nzPG<<"\n";


		//cout << "\n Allocating memory..."<<"\n";		
		PD->nzmax = nzPD; PL->nzmax = nzPL; PU->nzmax = nzPU;PG->nzmax = nzPG;
		PD->i = new int[nzPD];PL->i = new int[nzPL];PU->i = new int[nzPU];PG->i = new int[nzPG];
		PD->x = new double[nzPD];PL->x = new double[nzPL];PU->x = new double[nzPU];PG->x = new double[nzPG];

		//cout << "\n Allocating PD and PL..."<<"\n";
		iterPD =0; iterPL = 0;
		//filling the blocks PD and PL
		for(j=rD1;j<rD2+ol;j++)
		{
			//cout << "\n j : "<< j<<"\n";
			
			PD->p[j-rD1] = iterPD;
			for(k = D->p[j];k<D->p[j+1];k++)
			{
				if(D->i[k] >= rD1 && D->i[k] < rD2+ol) 
				{
					PD->i[iterPD] = D->i[k]-rD1;
					PD->x[iterPD] = D->x[k];
					iterPD += 1;
				}
				else continue;
			}
			
			PL->p[j-rD1] = iterPL;
			for(l = L->p[j];l<L->p[j+1];l++)
			{
				if(L->i[l] >= r1 && L->i[l] < r2+ol)
				{
					PL->i[iterPL] = L->i[l] -r1;
					PL->x[iterPL] = L->x[l];
					iterPL += 1;
				}
				else continue;
			}
		}
		PD->p[PD->n] = iterPD; PL->p[PL->n] = iterPL;
		//cout << "\n PL->i[0] : " << PL->i[0] << "\n";
		//cout << "\n PL->p[ncol] :" << PL->p[PL->n] << "\n";
		//cout << "\n Diff PD : "<< nzPD - iterPD << "\n";
		//cout << "\n Diff PL : "<< nzPL - iterPL << "\n";

		status = umfpack_di_transpose(PL->m,PL->n, PL->p,PL->i,PL->x,null,null,PU->p,PU->i,PU->x) ;
		//cout << "\n Transpose status : "<< status << "\n";

		iterPG = 0;
		for(j = r1;j<r2+ol;j++)
		{
			PG->p[j-r1] = iterPG;
			for(l = G->p[j];l< G->p[j+1];l++)
			{
				//cout << "\n G->i[j] "<< G->i[j] << "\n"; 
				
				if(G->i[l] >= r1 && G->i[l] < r2+ol) 
				{
					PG->i[iterPG] = G->i[l] - r1;
					PG->x[iterPG] = G->x[l];
					iterPG += 1;
				}
				else continue;
			}
		}
		PG->p[PG->n] = iterPG;
		//cout << "\n MSC G->x[0]: "<< G->x[0] << "\n";
		//cout << "\n MSC PG->p[0]: "<< PG->p[0] << "\n";
		//cout << "\n Diff PG : "<< nzPG - iterPG << "\n";

		//compute the full schur complement of the extracted blocks
		compute_full_schur_complement(r1,r2,rD1,rD2,ol,MSC,PD,PL,PU,PG,&total_nz);

		//cout << "\n MSC->p[0] : "<< MSC->p[0] << "\n";
		//cout << "\n MSC->nnz : "<< total_nz << "\n";
		//cout << "\n MSC->x[total_nz-1] : "<< MSC->x[total_nz-1] << "\n";

		//updating the coordinates
		r1 = r2;  r2 = r2+sz;
		rD1 = rD2; rD2 = rD2 + szD;

		delete [] PD->p;delete [] PL->p;delete [] PU->p;delete [] PG->p;
		delete [] PD->i;delete [] PL->i;delete [] PU->i;delete [] PG->i;
		delete [] PD->x;delete [] PL->x;delete [] PU->x;delete [] PG->x;
		delete PD;delete PL;delete PU;delete PG;
	}

	//i = i + 1;
	//cout << "\n i : " << i <<"\n";
	if(i == nmsc_block - 2)
	{
		//cout << "\n i : " << i << "\n";
		nzPD = 0; nzPL = 0; nzPU = 0; nzPG = 0;
		cs_di* PD = new cs_di; cs_di* PL = new cs_di; cs_di* PU = new cs_di; cs_di* PG = new cs_di;
		PD->m = (rD2+oll) -rD1; PD->n = (rD2+oll) -rD1;
		PL->m = (r2+oll) -r1; PL->n = (rD2+oll) -rD1;  
		PU->m = (rD2+oll) -rD1; PU->n = (r2+oll) -r1;
		PG->m = (r2+oll) -r1; PG->n = (r2+oll) -r1;
		PD->p = new int[PD->n+1];PL->p = new int[PL->n+1];PU->p = new int[PU->n+1];PG->p = new int[PG->n+1];
		PD->nz = -1; PL->nz = -1; PU->nz = -1;PG->nz = -1;

		//cout << "\nnnz count\n";
		for(j = rD1;j<rD2+oll;j++)
		{
			for(k = D->p[j];k<D->p[j+1];k++)
			{
				if(D->i[k] >= rD1 && D->i[k] < rD2+oll) ++nzPD;
				else continue;
			}
			for(l = L->p[j];l<L->p[j+1];l++)
			{
				if(L->i[l] >= r1 && L->i[l] < r2+oll) ++nzPL;
				else continue;
			}
		}
		//cout << "\nnnz count U and G\n";
		for(j = r1;j<r2+oll;j++)
		{
			for(k = U->p[j];k<U->p[j+1];k++)
			{
				if(U->i[k] >= rD1 && U->i[k] < rD2+oll) ++nzPU;
				else continue;
			}
			for(l = G->p[j];l< G->p[j+1];l++)
			{
				//cout << "\n G->i[j] "<< G->i[j] << "\n"; 
				if(G->i[l] >= r1 && G->i[l] < r2+oll) ++nzPG;
				else continue;
			}
		}
		//cout << "\n nzPD : "<< nzPD<<"\tnzPL :"<<nzPL<<"\tnzPU :"<<nzPU<<"\tnzPG :"<<nzPG<<"\n";


		//cout << "\n Allocating memory..."<<"\n";		
		PD->nzmax = nzPD; PL->nzmax = nzPL; PU->nzmax = nzPU;PG->nzmax = nzPG;
		PD->i = new int[nzPD];PL->i = new int[nzPL];PU->i = new int[nzPU];PG->i = new int[nzPG];
		PD->x = new double[nzPD];PL->x = new double[nzPL];PU->x = new double[nzPU];PG->x = new double[nzPG];

		//cout << "\n Allocating PD and PL..."<<"\n";
		iterPD =0; iterPL = 0;
		//filling the blocks PD and PL
		for(j=rD1;j<rD2+oll;j++)
		{
			//cout << "\n j : "<< j<<"\n";
			PD->p[j-rD1] = iterPD;
			for(k = D->p[j];k<D->p[j+1];k++)
			{
				if(D->i[k] >= rD1 && D->i[k] < rD2+oll) 
				{
					PD->i[iterPD] = D->i[k]-rD1;
					PD->x[iterPD] = D->x[k];
					iterPD += 1;
				}
				else continue;
			}
			
			PL->p[j-rD1] = iterPL;
			for(l = L->p[j];l<L->p[j+1];l++)
			{

				if(L->i[l] >= r1 && L->i[l] < r2+oll)
				{
					PL->i[iterPL] = L->i[l] -r1;
					PL->x[iterPL] = L->x[l];
					iterPL += 1;
				}
				else continue;
				//cout << "\n r1 : "<<r1 <<"\t"<<"L->i[l] : "<< L->i[l] << "\n";
				//cout << "\n r2 : "<<r2 <<"\t"<<"L->i[l] : "<< L->i[l] << "\n"; 
			}
		}
		PD->p[PD->n] = iterPD; PL->p[PL->n] = iterPL;
		//cout << "\n PL->nnz :" << nzPL << "\n";
		//cout << "\n iterPL :" << iterPL << "\n";
		//cout << "\n Diff PD : "<< nzPD - iterPD << "\n";
		//cout << "\n Diff PL : "<< nzPL - iterPL << "\n";
		//cout << "\n PL->p[0] : " << PL->i[0] << "\n";
		//cout << "\n PL->p[sizeD] : " << PD->p[0] << "\n";

		status = umfpack_di_transpose(PL->m,PL->n, PL->p,PL->i,PL->x,null,null,PU->p,PU->i,PU->x) ;
		//cout << "\n Transpose status : "<< status << "\n";

		iterPG = 0;
		for(j = r1;j<r2+oll;j++)
		{
			PG->p[j-r1] = iterPG;
			for(l = G->p[j];l< G->p[j+1];l++)
			{
				//cout << "\n G->i[j] "<< G->i[j] << "\n"; 
				
				if(G->i[l] >= r1 && G->i[l] < r2+oll) 
				{
					PG->i[iterPG] = G->i[l] - r1;
					PG->x[iterPG] = G->x[l];
					iterPG += 1;
				}
				else continue;
			}
		}
		PG->p[PG->n] = iterPG;
		//cout << "\n MSC G->x[0]: "<< G->x[0] << "\n";
		//cout << "\n MSC PG->p[0]: "<< PG->p[0] << "\n";
		//cout << "\n Diff PG : "<< nzPG - iterPG << "\n";

		//compute the full schur complement of the extracted blocks
		compute_full_schur_complement(r1,r2,rD1,rD2,oll,MSC,PD,PL,PU,PG,&total_nz);

		//cout << "\n MSC->p[0] : "<< MSC->p[0] << "\n";
		//cout << "\n MSC->nnz : "<< total_nz << "\n";
		//cout << "\n MSC->x[total_nz-1] : "<< MSC->x[total_nz-1] << "\n";

		//updating the coordinates
		r1 = r2;  
		rD1 = rD2; 
		//cout << "\n 2nd last block done! \n";
		delete [] PD->p;delete [] PL->p;delete [] PU->p;delete [] PG->p;
		delete [] PD->i;delete [] PL->i;delete [] PU->i;delete [] PG->i;
		delete [] PD->x;delete [] PL->x;delete [] PU->x;delete [] PG->x;
		delete PD;delete PL;delete PU;delete PG;
	}

	i = i+1;  //cout << "\n i :"<< i << "\n";
	r2 = sizeG; rD2 = sizeD;

	if(i == nmsc_block-1)
	{
		//cout << "\n i : " << i << "\n";
		nzPD = 0; nzPL = 0; nzPU = 0; nzPG = 0;
		cs_di* PD = new cs_di; cs_di* PL = new cs_di; cs_di* PU = new cs_di; cs_di* PG = new cs_di;
		PD->m = rD2 -rD1; PD->n = rD2 -rD1;
		PL->m = r2 -r1; PL->n = rD2 -rD1;  
		PU->m = rD2 -rD1; PU->n = r2 -r1;
		PG->m = r2 -r1; PG->n = r2 -r1;
		PD->p = new int[PD->n+1];PL->p = new int[PL->n+1];PU->p = new int[PU->n+1];PG->p = new int[PG->n+1];
		PD->nz = -1; PL->nz = -1; PU->nz = -1;PG->nz = -1;

		for(j = rD1;j<rD2;j++)
		{
			for(k = D->p[j];k<D->p[j+1];k++)
			{
				if(D->i[k] >= rD1 && D->i[k] < rD2) ++nzPD;
				else continue;
			}
			for(l = L->p[j];l<L->p[j+1];l++)
			{
				if(L->i[l] >= r1 && L->i[l] < r2) ++nzPL;
				else continue;
			}
		}

		for(j = r1;j<r2;j++)
		{
			for(k = U->p[j];k<U->p[j+1];k++)
			{
				if(U->i[k] >= rD1 && U->i[k] < rD2) ++nzPU;
				else continue;
			}
			for(l = G->p[j];l< G->p[j+1];l++)
			{
				//cout << "\n G->i[j] "<< G->i[j] << "\n"; 
				if(G->i[l] >= r1 && G->i[l] < r2) ++nzPG;
				else continue;
			}
		}
		//cout << "\n nzPD : "<< nzPD<<"\tnzPL :"<<nzPL<<"\tnzPU :"<<nzPU<<"\tnzPG :"<<nzPG<<"\n";


		//cout << "\n Allocating memory..."<<"\n";		
		PD->nzmax = nzPD; PL->nzmax = nzPL; PU->nzmax = nzPU;PG->nzmax = nzPG;
		PD->i = new int[nzPD];PL->i = new int[nzPL];PU->i = new int[nzPU];PG->i = new int[nzPG];
		PD->x = new double[nzPD];PL->x = new double[nzPL];PU->x = new double[nzPU];PG->x = new double[nzPG];

		//cout << "\n Allocating PD and PL..."<<"\n";
		iterPD =0; iterPL = 0;
		//filling the blocks PD and PL
		for(j=rD1;j<rD2;j++)
		{
			//cout << "\n j : "<< j<<"\n";
			PD->p[j-rD1] = iterPD;
			for(k = D->p[j];k<D->p[j+1];k++)
			{
				if(D->i[k] >= rD1 && D->i[k] < rD2) 
				{
					PD->i[iterPD] = D->i[k]-rD1;
					PD->x[iterPD] = D->x[k];
					iterPD += 1;
				}
				else continue;
			}
			
			PL->p[j-rD1] = iterPL;
			for(l = L->p[j];l<L->p[j+1];l++)
			{
				if(L->i[l] >= r1 && L->i[l] < r2)
				{
					PL->i[iterPL] = L->i[l] -r1;
					PL->x[iterPL] = L->x[l];
					iterPL += 1;
				}
				else continue;
			}
		}
		PD->p[PD->n] = iterPD; PL->p[PL->n] = iterPL;
		//cout << "\n PL->p[ncol] :" << PL->p[PL->n] << "\n";
		//cout << "\n Diff PD : "<< nzPD - iterPD << "\n";
		//cout << "\n Diff PL : "<< nzPL - iterPL << "\n";

		status = umfpack_di_transpose(PL->m,PL->n, PL->p,PL->i,PL->x,null,null,PU->p,PU->i,PU->x) ;
		//cout << "\n Transpose status : "<< status << "\n";

		iterPG = 0;
		for(j = r1;j<r2;j++)
		{
			PG->p[j-r1] = iterPG;
			for(l = G->p[j];l< G->p[j+1];l++)
			{
				//cout << "\n G->i[j] "<< G->i[j] << "\n"; 
				
				if(G->i[l] >= r1 && G->i[l] < r2) 
				{
					PG->i[iterPG] = G->i[l] - r1;
					PG->x[iterPG] = G->x[l];
					iterPG += 1;
				}
				else continue;
			}
		}
		PG->p[PG->n] = iterPG;
		//cout << "\n MSC G->x[0]: "<< G->x[0] << "\n";
		//cout << "\n MSC PG->p[0]: "<< PG->p[0] << "\n";
		//cout << "\n Diff PG : "<< nzPG - iterPG << "\n";

		//compute the full schur complement of the extracted blocks
		compute_full_schur_complement(r1,r2,rD1,rD2,0,MSC,PD,PL,PU,PG,&total_nz);

		MSC->p[sizeG] = total_nz;
		//cout << "\n MSC->p[0] : "<< MSC->p[0] << "\n";
		cout << "\n MSC->nnz : "<< MSC->p[sizeG] << "\n";
		//cout << "\n MSC->x[total_nz-1] : "<< MSC->x[total_nz-1] << "\n";


		delete [] PD->p;delete [] PL->p;delete [] PU->p;delete [] PG->p;
		delete [] PD->i;delete [] PL->i;delete [] PU->i;delete [] PG->i;
		delete [] PD->x;delete [] PL->x;delete [] PU->x;delete [] PG->x;
		delete PD;delete PL;delete PU;delete PG;
	}

	cout << "\n Mini Schur Complement computation done! \n";

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
	string rhs_filename = "49/1/JTe49_1.txt";

	//string col_filename = prefix + "_col.txt";
	//ifstream in_file(col_filename,ios::in);

	// getting the matrix size: n -> no of cols, ncc -> no of non zeros
	cc_header_read ( prefix, ncc, num_cols );
	num_cols = num_cols+1;                          //DOUBT!!!!!!!!
	cout << "\nNo of non zeros = "<< ncc << "\n";
	cout << "\nNo of columns = "<< num_cols << "\n";

	//size of the D matrix
	sizeD = num_cols - sizeG; 
	//cout << "\nsizeD = "<< sizeD << "\n";

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
	MSC->i = new int[sizeG*sizeG];
	MSC->x = new double[sizeG*sizeG];

	//setting values
	D->nzmax = nzD;
	L->nzmax = nzL; //cout << "\n L nnz :"<< L->nzmax << "\n";
	U->nzmax = nzL;
	G->nzmax = nzG;  
	MSC->nzmax = sizeG*sizeG;

	
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
	

	//cout << "\n nzD = " << D->p[sizeD] << "\n";
	//cout << "\n nzL = " << L->p[sizeD] << "\n";
	

	//status = umfpack_di_transpose (n_row, n_col, Ap, Ai, Ax, P, Q, Rp, Ri, Rx) ;
	status = umfpack_di_transpose(sizeG,sizeD, L->p,L->i,L->x,null,null,U->p,U->i,U->x) ;

	//cout << "\n TRANSPOSE STATUS : "<< status<< "\n";
	//cout << "\nD[9][9](591.616) = "<< L->x[0]<<"\n";
	//cout << "\nL[23714][9](335.400) = "<< U->x[0]<<"\n";
	//cout << "\n Diff D: "<<iterD - D->nzmax<<"\n";
	//cout << "\n Diff L: "<<iterL - L->nzmax<<"\n";
	//cout << "\n Diff U: "<<iterU - U->nzmax<<"\n";

	for(j = sizeD;j<num_cols;j++)
	{
		G->p[j - sizeD] = iterG;
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
					iterG += 1;
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
					iterG += 1;
				}
				//cout << "\n k : " << k << "\n";
				//break;
			}
		}
	    //break;
	}

	G->p[sizeG] = iterG;
	//cout << "\n Col Count : " << col_count << "\n";
	//cout << "\n Diff G : "<<iterG - G->nzmax<<"\n";
	//cout << "\n U Last element : "<<L->x[54975]<<"\n";
	cout << "\nFilling non zeros complete!!\n";
	
	
	//compute the mini schur complement in the Ccs_di format.
	compute_mini_schur_complement(A,MSC,D,L,U,G);

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