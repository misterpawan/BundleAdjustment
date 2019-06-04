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
#include "sort.h"
#include "mkl.h"

#define sizeG 441 //hardcoding it for the time being
#define size_MKL_IPAR 128

/* This function densifies the jth column of input matrix A
*/
double* densify_column(cs_di* A, int j)
{
	//cout << "\n In densify....! \n";
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
	int k;
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
	int sym_status,num_status,solve_status;
	//int nz_count;

	//  From the matrix data, create the symbolic factorization information.
	sym_status = umfpack_di_symbolic ( PD->m, PD->n, PD->p, PD->i, PD->x, &Symbolic, null, null );
	//cout << "\n Symbolic status :" << sym_status << "\n";

	//  From the symbolic factorization information, carry out the numeric factorization.
	num_status = umfpack_di_numeric ( PD->p, PD->i, PD->x, Symbolic, &Numeric, null, null );
	//cout << "\n Numeric status :" << num_status << "\n";

	for(j=0;j<PU->n;j++)
	//for(j=0;j<1;j++)
	{
		b = new double[PU->m]; //vector to store the densified rhs at each iteration
		x1 = new double[PD->n]; //vector to store the solution at each iteration

		b = densify_column(PU,j); // densifies the jth column of PU 

  		
  		//  Using the numeric factorization, solve the linear system.
  		solve_status = umfpack_di_solve ( UMFPACK_A, PD->p, PD->i, PD->x, x1, b, Numeric, null, null );
  		//cout << "\n Solve status :" << solve_status << "\n";
	
		
  		//Store x1 in the jth column of AA
  		make_x_sparse(AA,j,x1);

		delete []  x1; delete [] b;
	}

	//  Free the symbolic factorization memory.
  	umfpack_di_free_symbolic ( &Symbolic );
  	//  Free the numeric factorization.
  	umfpack_di_free_numeric ( &Numeric );

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

	//cout << "\nAA->p[AA->n] : " << AA->p[AA->n] << "\n";


	LDU = cs_di_multiply(PL,AA);
	

	delete[] AA->p; delete[] AA->i; delete[] AA->x;
	delete AA; 

	
	cs_di* S ;//= new cs_di;
	//S->nz = -1;
	//S->m = PG->m; S->n = PG->n; S->nzmax = S->m * S->n;
	//S->p = new int[S->n+1]; S->i = new int[S->nzmax]; S->x = new double[S->nzmax];


	
	S = cs_di_add(PG,LDU,1.0,-1.0); //subtraction
	//S = cs_add(PG,PG,1.0,-1.0); //subtraction
	//cout << "S->nzmax : " << S->nzmax << "\n";
	//cout << "S->nzmax : " << S->p[10] << "\n";

	int j,k,l;

	//filling up the msc block
	
	for(j=0;j< S->n; j++)
	{
		if(r1 == 0 && j == 0) MSC->p[0] = 0;
		else if(j==0) continue;
		else
			MSC->p[r1+j] = (S->p[j]-S->p[j-1])+MSC->p[r1+j-1];

		if(j == S->n-1) 
		{	
			MSC->p[r1+j+1] = (S->p[S->n] - S->p[S->n-1]) + MSC->p[r1+j];
			//cout << "\nMSC->p["<<r1+j+1<<"] : " <<  MSC->p[r1+j+1] << "\n";
		}	
	}	
	
	l = 0;
	for(j = r1; j < r2+ol; j++)
	{
		if(j < sizeG-1)
		{
			//cout <<  "\n MSC->p[j+1] : "<< MSC->p[j+1] << "\n";
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
			//for(k = MSC->p[j];k  < MSC->p[j+1] && l < S->nzmax; k++)
			{
				MSC->i[k] = S->i[l] + r1;
				MSC->x[k] = S->x[l];
				l++;
			}
		}
	}
	*total_nz += S->nzmax; //cout << "\nTotal nnz : " << *total_nz << "\n";
	//cout << "\n Total nnz : "<< S->nzmax << "\n";
	//cout << "\n Last element : "<< S->x[S->nzmax-1] << "\n";

	for(j=r1;j< r2+ol; j++)
	{
		if(j < sizeG-1)
		{
			sort_rows_val(&(MSC->i[MSC->p[j]]),&(MSC->x[MSC->p[j]]),(MSC->p[j+1]-MSC->p[j]));	
		}
		else
		{
			sort_rows_val(&(MSC->i[MSC->p[j]]),&(MSC->x[MSC->p[j]]),((MSC->p[r1]+S->nzmax)-MSC->p[j]));
		}
	}
//	cout << "\nSorting done!\n";
/*
	for(k = MSC->p[0]; k < MSC->p[1]; k++)
	{
		//if(MSC->i[k] > 420 &&  MSC->i[k] < 440)
			cout << "MSC->i["<<k<<"]: "<< MSC->i[k] << "\n";
	}
*/
	delete [] LDU->p; delete [] LDU->i; delete [] LDU->x;
	delete [] S->p; delete [] S->i; delete [] S->x;
	delete S;
	delete LDU; 

	return;
}


/* This function computes the mini schur complement of the input matrix.
PARAMS : 
	A,D,L,U : INPUT, the blocks obtained after domain decomposition
	MSC : OUTPUT,  the mini schur complement block computed.

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
		//cout << "\n Block 1 : r1 : "<<r1 << " "<< "r2 : "<< r2 << "\n";
		//compute the full schur complement of the extracted blocks
		compute_full_schur_complement(r1,r2,rD1,rD2,ol,MSC,PD,PL,PU,PG,&total_nz);

		//cout << "\n MSC->p[0] : "<< MSC->p[0] << "\n";
		//cout << "\n MSC->nnz : "<< total_nz << "\n";
		//cout << "\n MSC->x[total_nz-1] : "<< MSC->x[total_nz-1] << "\n";

		//updating the coordinates
		r1 = r2;  r2 = r2+sz;
		rD1 = rD2; rD2 = rD2 + szD;
		//cout << "\n1st k blocks done!\n";
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
		//cout << "\n Block 2 : r1 : "<<r1 << " "<< "r2 : "<< r2 << "\n";
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
		//cout << "\n Block 3 : r1 : "<<r1 << " "<< "r2 : "<< r2 << "\n";
		//compute the full schur complement of the extracted blocks
		compute_full_schur_complement(r1,r2,rD1,rD2,0,MSC,PD,PL,PU,PG,&total_nz);

		MSC->p[sizeG] = total_nz;
		//cout << "\n MSC->p[0] : "<< MSC->p[0] << "\n";
		//cout << "\n MSC->nnz : "<< MSC->p[sizeG] << "\n";
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
	double *solve_null = ( double * ) NULL;
	void *Numeric, *Numeric_D,*Numeric_MSC;
	string prefix = "49/1/JTJ49_1";
	double r;
	int status,sym_status,num_status,solve_status;
	void* Symbolic,*Symbolic_D,*Symbolic_MSC;
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

	//cout << "\nRHS file read complete!!\n";

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

	MSC->nz = -1;MSC->m = sizeG;MSC->n = sizeG;
	D->nz = -1;D->m = sizeD;D->n = sizeD;
	G->nz = -1;G->m = sizeG;G->n = sizeG;
	L->nz = -1;L->m = sizeG;L->n = sizeD;
	U->nz = -1;U->m = sizeD;U->n = sizeG;


	//cout << "\nCounting non zeros ...\n";
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
	//cout << "\nCounting non zeros complete!!\n";
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


	// Since L is not used anymore
	delete [] L->p;delete [] L->i;delete [] L->x;delete L;

	int ok = cs_di_sprealloc(MSC,MSC->p[sizeG]);
	//cout << "\n ok : "<< ok << "\n";

	//cout << "\nMSC ->nzmax : " << MSC->p[sizeG+1] << "\n";
	//cout << "\nMSC ->i[MSC->p[sizeG]] : " << MSC->i[MSC->p[401]-1] << "\n"; //rows are being stored in reverse

	//std::sort(&MSC->i[MSC->p[0]],&MSC->i[MSC->p[sizeG-1]]);

	/************LU Factorization of D and MSC******************************/

	sym_status = umfpack_di_symbolic ( D->m, D->n, D->p, D->i, D->x, &Symbolic_D, solve_null, solve_null );
	//cout << "\n Symbolic status for D :" << sym_status << "\n";

	num_status = umfpack_di_numeric ( D->p, D->i, D->x, Symbolic_D, &Numeric_D, solve_null, solve_null );
	//cout << "\n Numeric status for D:" << num_status << "\n";

	//  Free the symbolic factorization memory.
  	umfpack_di_free_symbolic ( &Symbolic_D );

	sym_status = umfpack_di_symbolic ( MSC->m, MSC->n, MSC->p, MSC->i, MSC->x, &Symbolic_MSC, solve_null, solve_null );
	//cout << "\n Symbolic status for MSC :" << sym_status << "\n";

	num_status = umfpack_di_numeric ( MSC->p, MSC->i, MSC->x, Symbolic_MSC, &Numeric_MSC, solve_null, solve_null );
	//cout << "\n Numeric status for MSC:" << num_status << "\n";
	umfpack_di_free_symbolic ( &Symbolic_MSC );
  	

	/***TESTING THE CORRECTNESS OF COMPUTED MSC WITH MATLAB***/
/*	double* MSC_b = new double[sizeG];
	string test_filename = "49/1/MSC_b.txt";
	void* MSC_Numeric;
	void* MSC_Symbolic;
	double* x_msc = new double[sizeG];


	r8vec_data_read ( test_filename, sizeG, MSC_b);

	//  From the matrix data, create the symbolic factorization information.
	sym_status = umfpack_di_symbolic ( MSC->m, MSC->n, MSC->p, MSC->i, MSC->x, &MSC_Symbolic, solve_null, solve_null );
	cout << "\n MSC Symbolic status :" << sym_status << "\n";
	//  From the symbolic factorization information, carry out the numeric factorization.
	num_status = umfpack_di_numeric ( MSC->p, MSC->i, MSC->x, MSC_Symbolic, &MSC_Numeric, solve_null, solve_null );
	cout << "\n MSC Numeric status :" << num_status << "\n";
	//  Using the numeric factorization, solve the linear system.
  	solve_status = umfpack_di_solve ( UMFPACK_A, MSC->p, MSC->i, MSC->x, x_msc, MSC_b, MSC_Numeric, solve_null, solve_null );
  	cout << "\n MSC Solve status :" << solve_status << "\n";

  	ofstream outfile("x_msc.txt");

  	if(outfile.is_open())
  	{
  		for(k = 0; k < sizeG; k++)
  			outfile << x_msc[k] << "\n";
  	}
  	outfile.close();

  	// The difference between x from matlab and c++ code here is 6.2413e-06.
*/


	/**********************GMRES CALL******************************/

	//initializing variables and data structures for DFGMRES call
	int restart = 20;  //DFGMRES restarts
	MKL_INT* ipar = new MKL_INT[size_MKL_IPAR];
	ipar[14] = restart;

	//cout << "\n tmp size : "<< num_cols*(2*ipar[14]+1)+ipar[14]*((ipar[14]+9)/2+1) << "\n";

	double* dpar = new double[size_MKL_IPAR]; 
	//double tmp[num_cols*(2*ipar[14]+1)+ipar[14]*((ipar[14]+9)/2+1)];
	double* tmp = new double[num_cols*(2*ipar[14]+1)+ipar[14]*((ipar[14]+9)/2+1)];
	//double expected_solution[num_cols];
	double* rhs = new double[num_cols];
	double* computed_solution = new double[num_cols];
	double* residual = new double[num_cols];   
	double nrm2;

	MKL_INT itercount,ierr=0;
	MKL_INT RCI_request, RCI_count, ivar;
	double dvar;
	char cvar,cvar1,cvar2;

	cout << "\nMKL var init done !\n";



	ivar = num_cols;
	cvar = 'N';  //no transpose

	/*---------------------------------------------------------------------------
	/* Save the right-hand side in vector rhs for future use
	/*---------------------------------------------------------------------------*/
	RCI_count=1;
	// JTe vector is not altered
	//rhs is used for residual calculations
	dcopy(&ivar, JTe, &RCI_count, rhs, &RCI_count);   

	/*---------------------------------------------------------------------------
	/* Initialize the initial guess
	/*---------------------------------------------------------------------------*/
	for(RCI_count=0; RCI_count<num_cols; RCI_count++)
	{
		computed_solution[RCI_count]=0.0;
	}
	computed_solution[0]=100.0;

	/*---------------------------------------------------------------------------
	/* Initialize the solver
	/*---------------------------------------------------------------------------*/
	dfgmres_init(&ivar, computed_solution, JTe, &RCI_request, ipar, dpar, tmp); 
	if (RCI_request!=0) goto FAILED;

	//cout << "\n RCI_request : " << RCI_request << "\n";
	//cout << "\n ipar[3]] : " << ipar[3] << "\n";
	//cout << "\n ipar[14]] : " << ipar[14] << "\n";
	ipar[7] = 0;
	ipar[4] = 20;  // Max Iterations
	ipar[10] = 1;  //Preconditioner used

	//cout << "\n dpar[0] : "<<dpar[0] << "\n";
	//cout << "\n dpar[1] : "<<dpar[1] << "\n";

	dpar[0] = 1.0e-04; //Relative Tolerance

	/*---------------------------------------------------------------------------
	/* Check the correctness and consistency of the newly set parameters
	/*---------------------------------------------------------------------------*/
	dfgmres_check(&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp); 
	if (RCI_request!=0) goto FAILED;

	//cout << "\n RCI_request : " << RCI_request << "\n";
	/*---------------------------------------------------------------------------
	/* Compute the solution by RCI (P)FGMRES solver with preconditioning
	/* Reverse Communication starts here
	/*---------------------------------------------------------------------------*/
	ONE:  cout << "\n in gmres... \n";  
	      cout << "\n computed_solution[0] : " << computed_solution[0] << "\n";
	      cout << "\n JTe[1] : " << JTe[1] << "\n";
	      cout << "\n ipar[0] : " << ipar[0] << "\n";
	      cout << "\n dpar[0] : " << dpar[0] << "\n";
	dfgmres(&ivar, computed_solution, JTe, &RCI_request, ipar, dpar, tmp);
	cout << "\n dfgmres RCI_request : "<<RCI_request << "\n";
	/*---------------------------------------------------------------------------
	/* If RCI_request=1, then compute the vector A*tmp[ipar[21]-1]
	/* and put the result in vector tmp[ipar[22]-1]	
	/*---------------------------------------------------------------------------
	/* NOTE that ipar[21] and ipar[22] contain FORTRAN style addresses,
	/* therefore, in C code it is required to subtract 1 from them to get C style
	/* addresses
	/*---------------------------------------------------------------------------*/
	/*---------------------------------------------------------------------------
	/* Since the input matrix is symmetric , so the CSC format for coefficient matrix
	   will be the same as CSR format with row and column arrays reversed. In the 
	   mat-vec function call below, the pointers for row and col arrays 
	   have been exchanged.
	/*-----------------------------------------------------------------------------*/
	/*------------------DEPRACATED ROUTINE (FIND ANOTHER )-------------------------*/
	if (RCI_request==1)
	{	
		//cout << "\n In mat-vec\n";
		//cout << "\n A->x[ncc-1]: " << A->x[ncc-1] << "\n";
		mkl_dcsrgemv(&cvar, &ivar, A->x, A->p, A->i, &tmp[ipar[21]-1], &tmp[ipar[22]-1]);
		//cout << "\n Mat-Vec done!!\n";
		//cout << "\n tmp[(ipar[22]-1)] : " << tmp[ipar[22]-1]<< "\n";
		goto ONE;
	}

	/*---------------------------------------------------------------------------
	/* If RCI_request=2, then do the user-defined stopping test
	/* The residual stopping test for the computed solution is performed here
	/*---------------------------------------------------------------------------
	/* NOTE: from this point vector b[N] is no longer containing the right-hand
	/* side of the problem! It contains the current FGMRES approximation to the
	/* solution. If you need to keep the right-hand side, save it in some other
	/* vector before the call to dfgmres routine. Here we saved it in vector
	/* rhs[N]. The vector b is used instead of rhs to preserve the
	/* original right-hand side of the problem and guarantee the proper
	/* restart of FGMRES method. Vector b will be altered when computing the
	/* residual stopping criterion!
	/*---------------------------------------------------------------------------*/
	if (RCI_request==2)
	{
		/* Request to the dfgmres_get routine to put the solution into b[N] via ipar[12]
		/*---------------------------------------------------------------------------
		/* WARNING: beware that the call to dfgmres_get routine with ipar[12]=0 at this stage may
		/* destroy the convergence of the FGMRES method, therefore, only advanced users should
		/* exploit this option with care */
		ipar[12]=1;
		/* Get the current FGMRES solution in the vector b[N] */
		dfgmres_get(&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp, &itercount);
		/* Compute the current true residual via MKL (Sparse) BLAS routines */
		mkl_dcsrgemv(&cvar, &ivar, A->x, A->p, A->i, rhs, residual); // A*x for new solution x
		dvar=-1.0E0;
		RCI_count=1;
		daxpy(&ivar, &dvar, JTe, &RCI_count, residual, &RCI_count);  // Ax - A*x_correct
		dvar=dnrm2(&ivar,residual,&i);
		if (dvar<1.0E-4) goto COMPLETE;   //taking tolerance as 1e-04

		else goto ONE;
	}

	/*---------------------------------------------------------------------------
	/* If RCI_request=3, then apply the preconditioner on the vector
	/* tmp[ipar[21]-1] and put the result in vector tmp[ipar[22]-1]
	/*---------------------------------------------------------------------------
	/* NOTE that ipar[21] and ipar[22] contain FORTRAN style addresses,
	/* therefore, in C code it is required to subtract 1 from them to get C style
	/* addresses
	/* Here is the recommended usage of the result produced by ILU0 routine
    /* via standard MKL Sparse Blas solver routine mkl_dcsrtrsv.
    /*---------------------------------------------------------------------------*/
	if (RCI_request==3)
	{
		cout << "\n Prec solve ..." << "\n";
		double y1[sizeD], y2[sizeG];
		double z1[sizeD], z2[sizeG],Lz1[sizeG];
		int prec_solve_status;
		double* null = (double*)NULL;
		int zvar = sizeG ; //no of rows of L
		int kk;

		for(kk = 0; kk<ivar; kk++)
		{
			if(kk < sizeD) y1[kk] = tmp[(ipar[21]-1) + kk];
			else y2[kk-sizeD] = tmp[(ipar[21]-1)+kk];  //splitting vector into y1,y2
		}

		prec_solve_status = umfpack_di_solve ( UMFPACK_A, D->p, D->i, D->x, z1, y1, Numeric_D, null, null );
  		cout << "\n  Prec solve status D :" << prec_solve_status << "\n";

  		

  		// performing L*z1 ... since U and L are transposes of each other, so
  		// using csc of U as CSR of L
  		mkl_dcsrgemv(&cvar, &zvar, U->x, U->p, U->i, z1, Lz1);
  		RCI_count = 1;
  		dvar = -1.0E0;
  		//y2 = y2 - L*z1
  		daxpy(&zvar, &dvar, Lz1, &RCI_count, y2, &RCI_count);  // Ax - A*x_correct

  		prec_solve_status = umfpack_di_solve ( UMFPACK_A, MSC->p, MSC->i, MSC->x, z2, y2, Numeric_MSC, null, null );
  		cout << "\n  Prec solve status MSC :" << prec_solve_status << "\n";

  		for(kk = 0; kk < ivar; kk++)
  		{
  			if(kk < sizeD) tmp[(ipar[22]-1)+kk] = z1[kk];
  			else tmp[(ipar[22]-1)+kk] = z2[kk - sizeD];
  		}

		goto ONE;
	}

	/*---------------------------------------------------------------------------
	/* If RCI_request=4, then check if the norm of the next generated vector is
	/* not zero up to rounding and computational errors. The norm is contained
	/* in dpar[6] parameter
	/*---------------------------------------------------------------------------*/
	if (RCI_request==4)
	{
		if (dpar[6]<1.0E-12) goto COMPLETE;
		else goto ONE;
	}
	/*---------------------------------------------------------------------------
	/* If RCI_request=anything else, then dfgmres subroutine failed
	/* to compute the solution vector: computed_solution[N]
	/*---------------------------------------------------------------------------*/
	else
	{
		goto FAILED;
	}

	/*---------------------------------------------------------------------------
	/* Reverse Communication ends here
	/* Get the current iteration number and the FGMRES solution (DO NOT FORGET to
	/* call dfgmres_get routine as computed_solution is still containing
	/* the initial guess!). Request to dfgmres_get to put the solution
	/* into vector computed_solution[N] via ipar[12]
	/*---------------------------------------------------------------------------*/
	COMPLETE:   ipar[12]=0;
	dfgmres_get(&ivar, computed_solution, JTe, &RCI_request, ipar, dpar, tmp, &itercount);
	cout << "The system has been solved \n";
//	cout << "\n RCI_request : "<< RCI_request << "\n";

	FAILED: cout << "The solver has returned the ERROR code " << RCI_request << "\n";


	//  Free the numeric factorization.
  	umfpack_di_free_numeric ( &Numeric_D );
  	umfpack_di_free_numeric ( &Numeric_MSC );

	delete [] JTe; 
	delete [] MSC->p; delete [] D->p;  delete [] G->p; delete [] U->p;
	delete [] D->i;  delete [] G->i; delete [] U->i;
	delete [] D->x; delete [] U->x; delete [] G->x;
	delete A; delete MSC; delete D;  delete G; delete U; 
	delete [] tmp; delete [] ipar; delete [] dpar;
	delete [] rhs; delete [] computed_solution; delete [] residual;

	return 0;
}