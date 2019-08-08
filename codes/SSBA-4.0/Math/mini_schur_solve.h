#ifndef V3D_MINI_SCHUR_SOLVE_H
#define V3D_MINI_SCHUR_SOLVE_H

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstring>
#include <ctime>
#include <cmath>
#include "Base/v3d_timer.h"

using namespace std;

#include "cs.h"
#include "umfpack.h"
#include "sort.h"
#include "mkl.h"

using namespace V3D;



//namespace V3D
//{
	#define sizeG 27612 
	#define size_MKL_IPAR 128
	//#define NUM_MSC_BLOCKS 50

	//This function densifies the jth column of input matrix PU
	void densify_column(cs_di* PU, int j,double* b)
	{
		int k;
		int  l = PU->p[j];
		int u = PU->p[j+1];
		for(int k=l; k < u; k++)
			b[PU->i[k]] = PU->x[k];

		return ;
	}

	//This function stores the dense array x1 into the jth column of AA as a sparse array
	void make_x_sparse(cs_di* AA,int j,double* x1)
	{
		int k;
		int row_arr_count; 
		int nz_count = 0;

		if(j==0) AA->p[j] = 0;
		
		row_arr_count = AA->p[j];

		for(k=0;k<AA->m;k++)
		{
			if(abs(x1[k]) >= 1e-16)
			{
				AA->i[row_arr_count+nz_count] = k;
				AA->x[row_arr_count+nz_count] = x1[k];  
				nz_count += 1;
			}
			else continue;
		}

		AA->p[j+1] = row_arr_count+nz_count;
		return;
	}

	/* This function does the computation PD\PU (MATLAB inverse operator) by solving 
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

		sym_status = umfpack_di_symbolic ( PD->m, PD->n, PD->p, PD->i, PD->x, &Symbolic, null, null );
		//cout << "\n Symbolic status :" << sym_status << "\n";

		num_status = umfpack_di_numeric ( PD->p, PD->i, PD->x, Symbolic, &Numeric, null, null );
		//cout << "\n Numeric status :" << num_status << "\n";

	  	umfpack_di_free_symbolic ( &Symbolic );

		for(j=0;j<PU->n;j++)
		{
			b = new double[PU->m](); //vector to store the densified rhs at each iteration
			x1 = new double[PD->n](); //vector to store the solution at each iteration

			//for(int k = 0; k < PU->m; k++)
            //	b[k] = 0.0;

			densify_column(PU,j,b); // densifies the jth column of PU and stores it in b

	  		
	  		//  Using the numeric factorization, solve the linear system.
	  		//  No need to solve if the jth column has no zeros. Simply set the solution vector x1=0 
	  		if((PU->p[j+1] - PU->p[j]) > 0)
	  			solve_status = umfpack_di_solve ( UMFPACK_A, PD->p, PD->i, PD->x, x1, b, Numeric, null, null );
	  		/*else
	  		{
	  			for(int k = 0; k < PD->n; k++)
	  				x1[k] = 0.0;
	  		}*/
	  		//cout << "\n Solve status :" << solve_status << "\n";
		
			
	  		//Store x1 in the jth column of AA
	  		make_x_sparse(AA,j,x1);

			delete []  x1; delete [] b;
		}

		
	  	umfpack_di_free_numeric ( &Numeric );

		return;
	}

	// This function computes the full Schur complement using the blocks 
	//given as input. The output is written in MSC.
	void compute_full_schur_complement(int r1,int r2,int rD1,int rD2,int ol,cs_di* MSC,cs_di* PD,cs_di* PL,cs_di* PU,cs_di* PG,int* total_nz)
	{
		int j,k,l,ok;
		cs_di* LDU;
		cs_di* S_tol ;
		cs_di* AA = new cs_di;
		
		AA->nz = -1;
		AA->m = PD->m; AA->n = PU->n; AA->nzmax = AA->m * AA->n;
		AA->p = new int[AA->n+1](); AA->i = new int[AA->nzmax](); AA->x = new double[AA->nzmax]();

		
		
		compute_PDinv_times_PU(PD,PU,AA);
		ok = cs_di_sprealloc(AA,AA->p[AA->n]);  //reallocate unused space
		// cout << "\nAA->nzmax : " << AA->nzmax << endl;
		
		LDU = cs_di_multiply(PL,AA);

		delete[] AA->p; delete[] AA->i; delete[] AA->x;
		delete AA; 

		

		S_tol = cs_di_add(PG,LDU,1.0,-1.0); // sparse matrix - matrix subtraction
		

		for(j=0;j< S_tol->n; j++)
		{
			if(j < (S_tol->n)-1)
			{
				sort_rows_val(&(S_tol->i[S_tol->p[j]]),&(S_tol->x[S_tol->p[j]]),(S_tol->p[j+1]-S_tol->p[j]));	
			}
			else
			{
				sort_rows_val(&(S_tol->i[S_tol->p[j]]),&(S_tol->x[S_tol->p[j]]),((S_tol->nzmax)-(S_tol->p[j])));
			}
		}
		/*
	   // keeping only values greater than 1e-16 in S 
		cs_di* S_tol = new cs_di;
		S_tol->nz = -1;
		S_tol->m = S->m; S_tol->n = S->n;
		S_tol->p = new int[S_tol->n+1];
		int nzStol = 0;
		S_tol->p[0] = 0;

		for(l=0; l < S->n;l++)
		{
			if(l < ((S->n)-1))
			{
				for(j = S->p[l]; j<S->p[l+1]; j++)
				{
					if(abs(S->x[j]) >= 1e-16)
					{	
						nzStol += 1;
					}
				}
				S_tol->p[l+1] = nzStol;
			}
			else
			{
				for(j = S->p[l]; j < S->nzmax; j++)
				{
					if(abs(S->x[j]) >= 1e-16)
					{
						nzStol += 1;
					}
				}
				S_tol->nzmax = nzStol;
			}
		}
		S_tol->p[S_tol->n] = S_tol->nzmax;
		S_tol->i = new int[S_tol->nzmax]; S_tol->x = new double[S_tol->nzmax];

		k =0;
		for(l=0; l < S->n; l++)
		{
			if(l < ((S->n)-1))
			{
				for(j = S->p[l]; j<S->p[l+1]; j++)
				{
					if(abs(S->x[j]) >= 1e-16)
					{
						S_tol->i[k] = S->i[j];
						S_tol->x[k] = S->x[j];
						k++;
					}
				}
			}
			else
			{
				for(j = S->p[l]; j < S->nzmax; j++)
				{
					if(abs(S->x[j]) >= 1e-16)
					{
						S_tol->i[k] = S->i[j];
						S_tol->x[k] = S->x[j];
						k++;
					}
				}
			}
		}
	
		delete [] S->p; delete [] S->i; delete [] S->x; delete S;
		*/

		//filling up the msc block
		for(j=0;j< S_tol->n; j++)
		{
			if(r1 == 0 && j == 0) MSC->p[0] = 0;
			else if(j==0) continue;
			else
				MSC->p[r1+j] = (S_tol->p[j]-S_tol->p[j-1])+MSC->p[r1+j-1];
			if(j == S_tol->n-1) 
			{	
				MSC->p[r1+j+1] = (S_tol->p[S_tol->n] - S_tol->p[S_tol->n-1]) + MSC->p[r1+j];
			}	
		}	


		l = 0;
		for(j = r1; j < r2+ol; j++)
		{
			if(j < sizeG-1)
			{
				for(k = MSC->p[j];k  < MSC->p[j+1] && l < S_tol->nzmax; k++)
				{
					MSC->i[k] = S_tol->i[l] + r1;
					MSC->x[k] = S_tol->x[l];
					l++;
				}
			}
			else
			{
				for(k = MSC->p[j];k  < MSC->nzmax && l < S_tol->nzmax; k++)
				{
					{
						MSC->i[k] = S_tol->i[l] + r1;
						MSC->x[k] = S_tol->x[l];
						l++;
					}
				}
			}
		}
		*total_nz += S_tol->nzmax; 
		//cout << "\nTotal nnz : " << *total_nz << "\n";
	

		delete [] LDU->p; delete [] LDU->i; delete [] LDU->x;
		delete [] S_tol->p; delete [] S_tol->i; delete [] S_tol->x;
		delete LDU;delete S_tol;


		return;
	}

	/* This function computes the mini schur complement of the input matrix.
		PARAMS : 
			A,D,L,U : INPUT, the blocks obtained after domain decomposition
			MSC : OUTPUT,  the mini schur complement block computed.

		The domain decomposition step is not shown as it is assumed that the size of 
		block G is available. 
	*/
	void compute_mini_schur_complement(cs_di* A,cs_di* MSC,cs_di* D,cs_di* L,cs_di* U,cs_di* G,int msc_block)
	{
		int nmsc_block = msc_block ;  //no of blocks for G
		int r = sizeG % nmsc_block;
		int sz = (sizeG - r)/nmsc_block ; //size of each nmsc block for G
		int r1 = 0,r2 = sz;   
		int nD = nmsc_block ; //no of blocks for D (same as that for G)
		int sizeD = D->m;
		int rD = sizeD % nD;
		int szD = (sizeD - rD)/nD; //size of each block for D
		int rD1 = 0,rD2 = szD; 


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


		//cout<< "\nComputing Mini Schur Complement...\n";
		for(i ; i < nmsc_block - 2;i++ )
		{
			nzPD = 0; nzPL = 0; nzPU = 0; nzPG = 0;
			cs_di* PD = new cs_di; cs_di* PL = new cs_di; cs_di* PU = new cs_di; cs_di* PG = new cs_di;
			PD->m = (rD2+ol) -rD1; PD->n = (rD2+ol) -rD1;
			PL->m = (r2+ol) -r1; PL->n = (rD2+ol) -rD1;  
			PU->m = (rD2+ol) -rD1; PU->n = (r2+ol) -r1;
			PG->m = (r2+ol) -r1; PG->n = (r2+ol) -r1;
			PD->p = new int[PD->n+1]();PL->p = new int[PL->n+1]();PU->p = new int[PU->n+1]();PG->p = new int[PG->n+1]();
			PD->nz = -1; PL->nz = -1; PU->nz = -1;PG->nz = -1;

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
					if(G->i[l] >= r1 && G->i[l] < r2+ol) ++nzPG;
					else continue;
				}
			}


			//cout << "\n Allocating memory..."<<"\n";		
			PD->nzmax = nzPD; PL->nzmax = nzPL; PU->nzmax = nzPU;PG->nzmax = nzPG;
			PD->i = new int[nzPD]();PL->i = new int[nzPL]();PU->i = new int[nzPU]();PG->i = new int[nzPG]();
			PD->x = new double[nzPD]();PL->x = new double[nzPL]();PU->x = new double[nzPU]();PG->x = new double[nzPG]();

			//cout << "\n Allocating PD and PL..."<<"\n";
			iterPD =0; iterPL = 0;
			
			//filling the blocks PD and PL
			for(j=rD1;j<rD2+ol;j++)
			{
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

			status = umfpack_di_transpose(PL->m,PL->n, PL->p,PL->i,PL->x,null,null,PU->p,PU->i,PU->x) ;
			//cout << "\n Transpose status : "<< status << "\n";

			iterPG = 0;
			for(j = r1;j<r2+ol;j++)
			{
				PG->p[j-r1] = iterPG;
				for(l = G->p[j];l< G->p[j+1];l++)
				{
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

			//compute the full schur complement of the extracted blocks
			compute_full_schur_complement(r1,r2,rD1,rD2,ol,MSC,PD,PL,PU,PG,&total_nz);


			//updating the coordinates
			r1 = r2;  r2 = r2+sz;
			rD1 = rD2; rD2 = rD2 + szD;

			delete [] PD->p;delete [] PL->p;delete [] PU->p;delete [] PG->p;
			delete [] PD->i;delete [] PL->i;delete [] PU->i;delete [] PG->i;
			delete [] PD->x;delete [] PL->x;delete [] PU->x;delete [] PG->x;
			delete PD;delete PL;delete PU;delete PG;
		}

		if(i == nmsc_block - 2)
		{
			nzPD = 0; nzPL = 0; nzPU = 0; nzPG = 0;
			cs_di* PD = new cs_di; cs_di* PL = new cs_di; cs_di* PU = new cs_di; cs_di* PG = new cs_di;
			PD->m = (rD2+oll) -rD1; PD->n = (rD2+oll) -rD1;
			PL->m = (r2+oll) -r1; PL->n = (rD2+oll) -rD1;  
			PU->m = (rD2+oll) -rD1; PU->n = (r2+oll) -r1;
			PG->m = (r2+oll) -r1; PG->n = (r2+oll) -r1;
			PD->p = new int[PD->n+1]();PL->p = new int[PL->n+1]();PU->p = new int[PU->n+1]();PG->p = new int[PG->n+1]();
			PD->nz = -1; PL->nz = -1; PU->nz = -1;PG->nz = -1;

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

			for(j = r1;j<r2+oll;j++)
			{
				for(k = U->p[j];k<U->p[j+1];k++)
				{
					if(U->i[k] >= rD1 && U->i[k] < rD2+oll) ++nzPU;
					else continue;
				}
				for(l = G->p[j];l< G->p[j+1];l++)
				{
					if(G->i[l] >= r1 && G->i[l] < r2+oll) ++nzPG;
					else continue;
				}
			}


			//cout << "\n Allocating memory..."<<"\n";		
			PD->nzmax = nzPD; PL->nzmax = nzPL; PU->nzmax = nzPU;PG->nzmax = nzPG;
			PD->i = new int[nzPD]();PL->i = new int[nzPL]();PU->i = new int[nzPU]();PG->i = new int[nzPG]();
			PD->x = new double[nzPD]();PL->x = new double[nzPL]();PU->x = new double[nzPU]();PG->x = new double[nzPG]();

			//cout << "\n Allocating PD and PL..."<<"\n";
			iterPD =0; iterPL = 0;

			//filling the blocks PD and PL
			for(j=rD1;j<rD2+oll;j++)
			{
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
				}
			}
			PD->p[PD->n] = iterPD; PL->p[PL->n] = iterPL;

			status = umfpack_di_transpose(PL->m,PL->n, PL->p,PL->i,PL->x,null,null,PU->p,PU->i,PU->x) ;
			//cout << "\n Transpose status : "<< status << "\n";

			iterPG = 0;
			for(j = r1;j<r2+oll;j++)
			{
				PG->p[j-r1] = iterPG;
				for(l = G->p[j];l< G->p[j+1];l++)
				{
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

			//compute the full schur complement of the extracted blocks
			compute_full_schur_complement(r1,r2,rD1,rD2,oll,MSC,PD,PL,PU,PG,&total_nz);

			

			//updating the coordinates
			r1 = r2;  
			rD1 = rD2; 

			delete [] PD->p;delete [] PL->p;delete [] PU->p;delete [] PG->p;
			delete [] PD->i;delete [] PL->i;delete [] PU->i;delete [] PG->i;
			delete [] PD->x;delete [] PL->x;delete [] PU->x;delete [] PG->x;
			delete PD;delete PL;delete PU;delete PG;
		}

		i = i+1;  
		r2 = sizeG; rD2 = sizeD;

		if(i == nmsc_block-1)
		{
			nzPD = 0; nzPL = 0; nzPU = 0; nzPG = 0;
			cs_di* PD = new cs_di; cs_di* PL = new cs_di; cs_di* PU = new cs_di; cs_di* PG = new cs_di;
			PD->m = rD2 -rD1; PD->n = rD2 -rD1;
			PL->m = r2 -r1; PL->n = rD2 -rD1;  
			PU->m = rD2 -rD1; PU->n = r2 -r1;
			PG->m = r2 -r1; PG->n = r2 -r1;
			PD->p = new int[PD->n+1]();PL->p = new int[PL->n+1]();PU->p = new int[PU->n+1]();PG->p = new int[PG->n+1]();
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
					if(G->i[l] >= r1 && G->i[l] < r2) ++nzPG;
					else continue;
				}
			}
			
			//cout << "\n Allocating memory..."<<"\n";		
			PD->nzmax = nzPD; PL->nzmax = nzPL; PU->nzmax = nzPU;PG->nzmax = nzPG;
			PD->i = new int[nzPD]();PL->i = new int[nzPL]();PU->i = new int[nzPU]();PG->i = new int[nzPG]();
			PD->x = new double[nzPD]();PL->x = new double[nzPL]();PU->x = new double[nzPU]();PG->x = new double[nzPG]();

			//cout << "\n Allocating PD and PL..."<<"\n";
			iterPD =0; iterPL = 0;
			
			//filling the blocks PD and PL
			for(j=rD1;j<rD2;j++)
			{
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

			status = umfpack_di_transpose(PL->m,PL->n, PL->p,PL->i,PL->x,null,null,PU->p,PU->i,PU->x) ;
			//cout << "\n Transpose status : "<< status << "\n";

			iterPG = 0;
			for(j = r1;j<r2;j++)
			{
				PG->p[j-r1] = iterPG;
				for(l = G->p[j];l< G->p[j+1];l++)
				{		
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
			
			//compute the full schur complement of the extracted blocks
			compute_full_schur_complement(r1,r2,rD1,rD2,0,MSC,PD,PL,PU,PG,&total_nz);

			MSC->p[sizeG] = total_nz;


			delete [] PD->p;delete [] PL->p;delete [] PU->p;delete [] PG->p;
			delete [] PD->i;delete [] PL->i;delete [] PU->i;delete [] PG->i;
			delete [] PD->x;delete [] PL->x;delete [] PU->x;delete [] PG->x;
			delete PD;delete PL;delete PU;delete PG;
		}

		return;
	}

	/* This function computes the preconditioner solve for the input array y_in
		and writes the output in z_out
	*/
	void prec_solve(cs_di *A,cs_di *D,cs_di *MSC,void *Numeric_D,void *Numeric_MSC,double *lcsr,int *il,int *jl,double *y_in,double *z_out)
	{
		double* y1 = new double[D->n]();
		double* y2 = new double[MSC->n]();
		double* z1 = new double[D->n]();
		double* z2 = new double[MSC->n]();
		double* Lz1 = new double[MSC->n]();
		int prec_solve_status;
		double* null = (double*)NULL;
		int zvar = MSC->n ; //no of rows of L
		int kk;
		char cvar = 'N';

		int ivar = A->n;
		int p = 1;
		double dvar = -1.0E0;


		for(kk = 0; kk<ivar; kk++)
		{
			if(kk < (D->n)) y1[kk] = y_in[kk];
			else y2[kk-(D->n)] = y_in[kk];  //splitting vector into y1,y2
		}

		prec_solve_status = umfpack_di_solve ( UMFPACK_A, D->p, D->i, D->x, z1, y1, Numeric_D, null, null );
		//cout << "\n Prec Solve status : " << prec_solve_status << "\n";
		
		mkl_dcsrgemv(&cvar, &zvar, lcsr, il, jl, z1, Lz1);
		
		//y2 = y2 - L*z1
		daxpy(&zvar, &dvar, Lz1, &p, y2, &p);  // Ax - A*x_correct

		prec_solve_status = umfpack_di_solve ( UMFPACK_A, MSC->p, MSC->i, MSC->x, z2, y2, Numeric_MSC, null, null );
		//cout << "\n  Prec solve status MSC :" << prec_solve_status << "\n";

		for(kk = 0; kk < ivar; kk++)
		{
			if(kk < (D->n)) z_out[kk] = z1[kk];
			else z_out[kk] = z2[kk - (D->n)];
		}


		delete [] y1; delete [] y2; delete [] z1; delete [] z2; delete [] Lz1;
		
		return;
	}

	 void mini_schur_solve(int num_cols,int ncc,int *colStarts,int *rowIdxs,double *values,double *Jt_e,
	 					   double *delta,double *prev_sol,int msc_block, double *MSC_time)
	 {	
		int i;
		int *null = ( int * ) NULL;
		double *solve_null = ( double * ) NULL;
		void *Numeric, *Numeric_D,*Numeric_MSC;
		double r;
		int status,sym_status,num_status,solve_status;
		void* Symbolic,*Symbolic_D,*Symbolic_MSC;
		int sizeD;
		int j,k;
		int nzD = 0 ,nzG = 0 ,nzL = 0; 
		int iterD = 0, iterL = 0 ,iterG = 0;
		int count = 1;


		//creating structures of cs_diPARSE
		cs_di* A = new cs_di;cs_di* MSC = new cs_di;cs_di* D = new cs_di;
		cs_di* L = new cs_di;cs_di* G = new cs_di;cs_di* U = new cs_di;


		//size of the D matrix
		sizeD = num_cols - sizeG; 

		A->nzmax = ncc; A->m = num_cols;A->n = num_cols;A->nz = -1;


		//Allocate space for the coefficient matrix
		A->x = new double[ncc](); A->p = new int[num_cols+1](); A->i = new int[ncc]();


		for(int q = 0; q < num_cols; q++) A->p[q] = colStarts[q];
		for(int q = 0; q < ncc; q++) A->i[q] = rowIdxs[q];
		for(int q = 0; q < ncc; q++) A->x[q] = values[q];

		A->p[num_cols] = ncc;

		//Allocating memory for blocks
		MSC->p = new int[sizeG+1]();
		D->p = new int[sizeD+1]();
		L->p = new int[sizeD+1]();
		G->p = new int[sizeG+1]();	
		U->p = new int[sizeG+1]();

		MSC->nz = -1;MSC->m = sizeG;MSC->n = sizeG;
		D->nz = -1;D->m = sizeD;D->n = sizeD;
		G->nz = -1;G->m = sizeG;G->n = sizeG;
		L->nz = -1;L->m = sizeG;L->n = sizeD;
		U->nz = -1;U->m = sizeD;U->n = sizeG;


		//cout << "\nCounting non zeros ...\n";
		for(j=0;j<sizeD;j++)
		{
			for(k=A->p[j]; k < A->p[j+1]; k++)
			{
				if(A->i[k] < sizeD) ++nzD;
				else ++nzL;
			}
		}

		
		nzG = ncc - (nzD + 2*nzL); 

		//Allocating memory
		D->i = new int[nzD](); D->x = new double[nzD]();
		L->i = new int[nzL](); L->x = new double[nzL]();
		U->i = new int[nzL](); U->x = new double[nzL]();
		G->i = new int[nzG](); G->x = new double[nzG]();
		MSC->i = new int[sizeG*sizeG](); MSC->x = new double[sizeG*sizeG]();

		//setting values
		D->nzmax = nzD; L->nzmax = nzL; U->nzmax = nzL; G->nzmax = nzG;MSC->nzmax = sizeG*sizeG;

		
		//cout << "\nFilling non zeros ...\n";
		for(j=0;j<sizeD;j++)
		{
			D->p[j] = iterD;
			L->p[j] = iterL;
			for(k=A->p[j]; k < A->p[j+1]; k++)
			{
				if(A->i[k] < sizeD) 
				{
					D->i[iterD] = A->i[k];
					D->x[iterD] = A->x[k];
					iterD += 1;
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
		L->p[sizeD] = iterL;
		
		status = umfpack_di_transpose(sizeG,sizeD, L->p,L->i,L->x,null,null,U->p,U->i,U->x) ;
		//cout << "\n TRANSPOSE STATUS : "<< status<< "\n";
		

		for(j = sizeD;j<num_cols;j++)
		{
			G->p[j - sizeD] = iterG;
			if(j < num_cols-1)
			{
				for(k = A->p[j]; k < A->p[j+1];k++)
				{
					if(A->i[k] < sizeD) continue;
					else
					{
						G->i[iterG] = A->i[k]-sizeD;
						G->x[iterG] = A->x[k];
						iterG += 1;
					}			
				}	
			}
			else 
			{
				for(k = A->p[j]; k < A->nzmax;k++)
				{
					if(A->i[k] < sizeD) continue;
					else
					{
						G->i[iterG] = A->i[k] - sizeD;
						G->x[iterG] = A->x[k];
						iterG += 1;
					}
				}
			}
		}

		G->p[sizeG] = iterG;

		//cout << "\nFilling non zeros complete!!\n";
		
		Timer t("MSC");
		t.start();
		//compute the mini schur complement in the cs_di format.
		compute_mini_schur_complement(A,MSC,D,L,U,G,msc_block);
		//cout << "\n Mini Schur Complement computation done! \n";
		t.stop();
		//cout << "\n Time for MSC construction : " << t.getTime() << endl;
		*MSC_time = t.getTime();

		// Since L is not used anymore
		delete [] U->p;delete [] U->i;delete [] U->x;delete U;

		int ok = cs_di_sprealloc(MSC,MSC->p[sizeG]);
		
		/*
		for(int k = 98; k < 99; k++)
		{
			//cout << "\n In loop .. " << endl;
			for(int l = MSC->p[k]; l<MSC->p[k+1]; l++)
				printf("\nMSC->i[%d] = %d\t\tMSC->x[%d] = %lf",l,MSC->i[l],l,MSC->x[l]);
				//cout << "\n l : " << l << endl;
		}
		*/

		/************LU Factorization of D and MSC******************************/

		sym_status = umfpack_di_symbolic ( D->m, D->n, D->p, D->i, D->x, &Symbolic_D, solve_null, solve_null );
		//cout << "\n Symbolic status for D :" << sym_status << "\n";

		num_status = umfpack_di_numeric ( D->p, D->i, D->x, Symbolic_D, &Numeric_D, solve_null, solve_null );
		//cout << "\n Numeric status for D:" << num_status << "\n";

	  	umfpack_di_free_symbolic ( &Symbolic_D );

		sym_status = umfpack_di_symbolic ( MSC->m, MSC->n, MSC->p, MSC->i, MSC->x, &Symbolic_MSC, solve_null, solve_null );
		//cout << "\n Symbolic status for MSC :" << sym_status << "\n";

		num_status = umfpack_di_numeric ( MSC->p, MSC->i, MSC->x, Symbolic_MSC, &Numeric_MSC, solve_null, solve_null );
		//cout << "\n Numeric status for MSC:" << num_status << "\n";
		umfpack_di_free_symbolic ( &Symbolic_MSC );
	  	
	  	delete [] G->p; delete [] G->i; delete [] G->x; delete G;

	  	//cout <<  "\n Starting MKL routines ... " << endl;
		/**********************GMRES CALL******************************/

		//initializing variables and data structures for DFGMRES call
		MKL_INT* ipar = new MKL_INT[size_MKL_IPAR]();

		double* dpar = new double[size_MKL_IPAR](); 
		
		double* tmp = new double[num_cols*(2*40+1)+(40*(40+9))/2+1]();
		double* rhs = new double[num_cols]();
		double* computed_solution = new double[num_cols]();
		double* residual = new double[num_cols]();   
		double nrm2,rhs_nrm,relres_nrm,dvar,relres_prev,prec_rhs_nrm,prec_relres_nrm;
		double *prec_rhs = new double[num_cols]();
		double tol = 1.0E-02;
		

		MKL_INT itercount,ierr=0;
		MKL_INT RCI_request, RCI_count, ivar;
		char cvar;

		//cout << "\nMKL var init done !\n";


		ivar = num_cols;
		cvar = 'N';  //no transpose
		
		/**********Converting A & L from CSC to CSR*****************/
	  	MKL_INT job[6] = {1,1,0,0,0,1};
	    double *acsr =  new double[ncc]();
	    double *lcsr =  new double[L->nzmax]();
	    MKL_INT *ja = new MKL_INT[ncc]();
	    MKL_INT *jl = new MKL_INT[L->nzmax]();
	    MKL_INT *ia = new MKL_INT[ivar+1]();
	    MKL_INT *il = new MKL_INT[sizeG+1]();
	    MKL_INT info;
	    MKL_INT lvar = sizeG;

	      //converting COO to CSR
	    mkl_dcsrcsc(job,&ivar,acsr,ja,ia,A->x,A->i,A->p,&info);
	    //cout << "\n Conversion info A : "<< info << "\n";


	    mkl_dcsrcsc(job,&lvar,lcsr,jl,il,L->x,L->i,L->p,&info);
	    //cout << "\n Conversion info L : "<< info << "\n";

		
		/*---------------------------------------------------------------------------
		/* Save the right-hand side in vector rhs for future use
		/*---------------------------------------------------------------------------*/
		RCI_count=1;
		/** extracting the norm of the rhs for computing rel res norm**/
		rhs_nrm = dnrm2(&ivar,Jt_e,&RCI_count); 
		//cout << "\n rhs_nrm : " << rhs_nrm << "\n";
		// Jt_e vector is not altered
		//rhs is used for residual calculations
		for(int q = 0; q < num_cols; q++) rhs[q] = Jt_e[q];
		// PRECONDITIONED RHS
		prec_solve(A,D,MSC,Numeric_D,Numeric_MSC,lcsr,il,jl,Jt_e,prec_rhs);
		//norm of the preconditioned rhs
		prec_rhs_nrm = dnrm2(&ivar,prec_rhs,&RCI_count);
		delete [] prec_rhs; 

		
		/*---------------------------------------------------------------------------
		/* Initialize the solver
		/*---------------------------------------------------------------------------*/
		dfgmres_init(&ivar, computed_solution, Jt_e, &RCI_request, ipar, dpar, tmp); 
		if (RCI_request!=0) goto FAILED;

		
		ipar[7] = 1;
		ipar[4] = 100;  // Max Iterations
		ipar[10] = 1;  //  Preconditioner used
		ipar[14] = 40; //  Internal iterations
		
		dpar[0] = tol; //Relative Tolerance

		/*---------------------------------------------------------------------------
		/* Initialize the initial guess
		/*---------------------------------------------------------------------------*/
		for(RCI_count=0; RCI_count<num_cols; RCI_count++)
		{
			computed_solution[RCI_count]=prev_sol[RCI_count];
		}
		
		/*---------------------------------------------------------------------------
		/* Check the correctness and consistency of the newly set parameters
		/*---------------------------------------------------------------------------*/
		dfgmres_check(&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp); 
		if (RCI_request!=0) goto FAILED;

		/*---------------------------------------------------------------------------
		/* Compute the solution by RCI (P)FGMRES solver with preconditioning
		/* Reverse Communication starts here
		/*---------------------------------------------------------------------------*/
		ONE:  dfgmres(&ivar, computed_solution, Jt_e, &RCI_request, ipar, dpar, tmp);

		
		if(RCI_request==0) goto COMPLETE;

		/*---------------------------------------------------------------------------
		/* If RCI_request=1, then compute the vector A*tmp[ipar[21]-1]
		/* and put the result in vector tmp[ipar[22]-1]	
		/*------------------DEPRECATED ROUTINE (FIND ANOTHER )-------------------------*/
		if (RCI_request==1)
		{
			mkl_dcsrgemv(&cvar, &ivar, acsr, ia, ja, &tmp[ipar[21]-1], &tmp[ipar[22]-1]);

			goto ONE;
		}

		/*---------------------------------------------------------------------------
		/* If RCI_request=2, then do the user-defined stopping test
		/* The residual stopping test for the computed solution is performed here*/
		if (RCI_request==2)
		{
			/* Request to the dfgmres_get routine to put the solution into b[N] via ipar[12]*/
			ipar[12]=1;
			
			/* Get the current FGMRES solution in the vector rhs[N] */
			dfgmres_get(&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp, &itercount);

			

			/* Compute the current true residual via MKL (Sparse) BLAS routines */
			mkl_dcsrgemv(&cvar, &ivar, acsr, ia, ja, rhs, residual); // A*x for new solution x
			dvar=-1.0E0;
			RCI_count=1;
			daxpy(&ivar, &dvar, Jt_e, &RCI_count, residual, &RCI_count);  // Ax - A*x_solution
			

			if(ipar[10] == 0)   // non preconditioned system
			{
				dvar=dnrm2(&ivar,residual,&RCI_count);
				relres_nrm = dvar/rhs_nrm;
				
			}
			else if(ipar[10] == 1)  //preconditioned system
			{
				double *prec_relres = new double[num_cols]();

				prec_solve(A,D,MSC,Numeric_D,Numeric_MSC,lcsr,il,jl,residual,prec_relres);
				prec_relres_nrm = dnrm2(&ivar,prec_relres,&RCI_count); 
				delete [] prec_relres;
				relres_nrm = prec_relres_nrm/prec_rhs_nrm;

			}

			if (relres_nrm<=tol) goto COMPLETE; 

			else goto ONE;
			
		}

		/*---------------------------------------------------------------------------
		/* If RCI_request=3, then apply the preconditioner on the vector
		/* tmp[ipar[21]-1] and put the result in vector tmp[ipar[22]-1]*/
		
		if (RCI_request==3)
		{
			prec_solve(A,D,MSC,Numeric_D,Numeric_MSC,lcsr,il,jl,&tmp[(ipar[21]-1)],&tmp[(ipar[22]-1)]);
		
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
		dfgmres_get(&ivar, computed_solution, Jt_e, &RCI_request, ipar, dpar, tmp, &itercount);
		//cout << "The system has been solved  in " << itercount << " iterations!\n";


		//store the solution into delta 
		RCI_count = 1;
		for(int q = 0; q < num_cols; q++) delta[q] = computed_solution[q];
	
		MKL_Free_Buffers();


		//  Free the numeric factorization.
	  	umfpack_di_free_numeric ( &Numeric_D );
	  	umfpack_di_free_numeric ( &Numeric_MSC );

		
		
		delete [] MSC->p; delete [] MSC->i;delete [] MSC->x;
		delete [] D->p;  delete [] D->i;delete [] D->x; 
		delete [] A->p; delete [] A->i; delete [] A->x; 
		delete [] L->p; delete [] L->i; delete [] L->x; 
		delete L; delete A;delete MSC; delete D;  

		delete [] tmp; delete [] ipar; delete [] dpar;
		delete [] rhs; delete [] computed_solution; delete [] residual;
		delete [] acsr; delete [] ia; delete [] ja;
		delete [] lcsr; delete [] il; delete [] jl;

		return ;

		FAILED: cout << "The solver has returned the ERROR code " << RCI_request << "\n";

		MKL_Free_Buffers();

		return ;

	    NOT_CONVERGE: cout << "The relative residual did not change for successive iterations !\n";

	    return ;
	 }
//}
#endif