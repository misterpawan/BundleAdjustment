#ifndef V3D_BLOCK_JACOBI_SOLVE_H
#define V3D_BLOCK_JACOBI_SOLVE_H

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstring>
#include <ctime>
#include <cmath>


#include "cs.h"
#include "umfpack.h"
#include "mkl.h"


using namespace std;
using namespace V3D;



namespace V3D
{
	//#define sizeG 3399
	#define size_MKL_IPAR 128
	//#define MAX_ITERS 500
	//#define TOL 1E-1
	//#define RESTARTS 50

	/* This function computes the preconditioner solve for the input array y_in
		and writes the output in z_out
	*/
	void prec_solve(cs_di *A,cs_di *D,cs_di *G,void *Numeric_D,void *Numeric_G,double *y_in,double *z_out)
	{
		double* y1 = new double[D->n]();
		double* y2 = new double[G->n]();
		double* z1 = new double[D->n]();
		double* z2 = new double[G->n]();
		int prec_solve_status;
		double* null = (double*)NULL;
		int zvar = G->n ; //no of rows of L
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

		// computing D_inv(y1)
		prec_solve_status = umfpack_di_solve ( UMFPACK_A, D->p, D->i, D->x, z1, y1, Numeric_D, null, null );
		//cout << "\n Prec Solve status : " << prec_solve_status << "\n";
		
		//computing G_inv(y2)
		prec_solve_status = umfpack_di_solve ( UMFPACK_A, G->p, G->i, G->x, z2, y2, Numeric_G, null, null );
		//cout << "\n  Prec solve status MSC :" << prec_solve_status << "\n";

		for(kk = 0; kk < ivar; kk++)
		{
			if(kk < (D->n)) z_out[kk] = z1[kk];
			else z_out[kk] = z2[kk - (D->n)];
		}


		delete [] y1; delete [] y2; delete [] z1; delete [] z2;
		
		return;
	}

	 void block_jacobi_solve(int num_cols,int ncc,int *colStarts,int *rowIdxs,
	 							double *values,double *Jt_e,double *delta,int *total_iters,double *LU_time,
	 							int max_gmres_iterations,int gmres_restarts,double tolerance,int sizeG)
	 {
		int i;
		int *null = ( int * ) NULL;
		double *solve_null = ( double * ) NULL;
		void *Numeric, *Numeric_D,*Numeric_G;
		double r;
		int status,sym_status,num_status,solve_status;
		void* Symbolic,*Symbolic_D,*Symbolic_G;
		int sizeD;
		int j,k;
		int nzD = 0 ,nzG = 0, nzL = 0 ; 
		int iterD = 0, iterG = 0;
		int count = 1;

		//creating structures of cs_diPARSE
		cs_di* A = new cs_di;cs_di* D = new cs_di;cs_di* G = new cs_di;


		//size of the D matrix
		sizeD = num_cols - sizeG; 

		A->nzmax = ncc; A->m = num_cols;A->n = num_cols;A->nz = -1;


		//Allocate space for the coefficient matrix
		A->x = new double[ncc](); A->p = new int[num_cols+1](); A->i = new int[ncc]();


		for(int q = 0; q < num_cols; q++) A->p[q] = colStarts[q];
		for(int q = 0; q < ncc; q++) A->i[q] = rowIdxs[q];
		for(int q = 0; q < ncc; q++) A->x[q] = values[q];
		//dcopy(&num_cols, values, &count, A->x, &count); 
		A->p[num_cols] = ncc;

		//Allocating memory for blocks
		D->p = new int[sizeD+1]();
		G->p = new int[sizeG+1]();	

		D->nz = -1;D->m = sizeD;D->n = sizeD;
		G->nz = -1;G->m = sizeG;G->n = sizeG;
		
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
		G->i = new int[nzG](); G->x = new double[nzG]();

		//setting values
		D->nzmax = nzD; G->nzmax = nzG;

		
		//cout << "\nFilling non zeros ...\n";
		for(j=0;j<sizeD;j++)
		{
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
			}
		}
		D->p[sizeD] = iterD;
		

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
						G->i[iterG] = A->i[k]-sizeD;
						G->x[iterG] = A->x[k];
						iterG += 1;
					}			
				}	
			}
			//for the last column
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


		/************LU Factorization of D and MSC******************************/

		sym_status = umfpack_di_symbolic ( D->m, D->n, D->p, D->i, D->x, &Symbolic_D, solve_null, solve_null );
		//cout << "\n Symbolic status for D :" << sym_status << "\n";

		num_status = umfpack_di_numeric ( D->p, D->i, D->x, Symbolic_D, &Numeric_D, solve_null, solve_null );
		//cout << "\n Numeric status for D:" << num_status << "\n";

		//  Free the symbolic factorization memory.
	  	umfpack_di_free_symbolic ( &Symbolic_D );

		sym_status = umfpack_di_symbolic ( G->m, G->n, G->p, G->i, G->x, &Symbolic_G, solve_null, solve_null );
		//cout << "\n Symbolic status for G :" << sym_status << "\n";

		num_status = umfpack_di_numeric ( G->p, G->i, G->x, Symbolic_G, &Numeric_G, solve_null, solve_null );
		//cout << "\n Numeric status for G:" << num_status << "\n";
		umfpack_di_free_symbolic ( &Symbolic_G );

		/**********************GMRES CALL******************************/

		//initializing variables and data structures for DFGMRES call
		MKL_INT* ipar = new MKL_INT[size_MKL_IPAR]();

		double* dpar = new double[size_MKL_IPAR](); 
		
		double* tmp = new double[num_cols*4]();
		//double* rhs = new double[num_cols]();
		double* computed_solution = new double[num_cols]();
		double* residual = new double[num_cols]();   
		double nrm2,rhs_nrm,relres_nrm,dvar,relres_prev,prec_rhs_nrm,prec_relres_nrm;
		double *prec_rhs = new double[num_cols]();
		double tol = tolerance;
		

		MKL_INT itercount,ierr=0;
		MKL_INT RCI_request, RCI_count, ivar;
		char cvar;

		//cout << "\nMKL var init done !\n";


		ivar = num_cols;
		cvar = 'N';  //no transpose
		
		/**********Converting A & L from CSC to CSR*****************/
	  	MKL_INT job[6] = {1,1,0,0,0,1};
	    double *acsr =  new double[ncc]();
	    MKL_INT *ja = new MKL_INT[ncc]();
	    MKL_INT *ia = new MKL_INT[ivar+1]();
	    MKL_INT info;
	    MKL_INT lvar = sizeG;

	    //converting CSC to CSR
	    mkl_dcsrcsc(job,&ivar,acsr,ja,ia,A->x,A->i,A->p,&info);
	    //cout << "\n Conversion info A : "<< info << "\n";
	    

		/*---------------------------------------------------------------------------
		/* Save the right-hand side in vector rhs for future use
		/*---------------------------------------------------------------------------*/
		RCI_count=1;
		/** extracting the norm of the rhs for computing rel res norm**/
		rhs_nrm = dnrm2(&ivar,Jt_e,&RCI_count); 
		//cout << "\n rhs_nrm : " << rhs_nrm << "\n";
		// Jt_e vector is not altered
		//rhs is used for residual calculations
		//for(int q = 0; q < num_cols; q++) rhs[q] = Jt_e[q];
		// PRECONDITIONED RHS
		prec_solve(A,D,G,Numeric_D,Numeric_G,Jt_e,prec_rhs);
		//norm of the preconditioned rhs
		prec_rhs_nrm = dnrm2(&ivar,prec_rhs,&RCI_count);
		delete [] prec_rhs; 

		//test_prec_solve(A, D,MSC,Numeric_D,Numeric_MSC,lcsr,il,jl);
		
		/*---------------------------------------------------------------------------
		/* Initialize the solver
		/*---------------------------------------------------------------------------*/
		//dfgmres_init(&ivar, computed_solution, Jt_e, &RCI_request, ipar, dpar, tmp); 
		dcg_init(&ivar, computed_solution, Jt_e, &RCI_request, ipar, dpar, tmp);
		if (RCI_request!=0) goto FAILED;

		
		ipar[7] = 1;
		ipar[4] = MAX_ITERS;  // Max Iterations
		ipar[10] = 1;  //  Preconditioner used
		//ipar[14] = RESTARTS; //  Internal iterations
		
		dpar[0] = tol; //Relative Tolerance

		/*---------------------------------------------------------------------------
		/* Initialize the initial guess
		/*---------------------------------------------------------------------------*/
		for(RCI_count=0; RCI_count<num_cols; RCI_count++)
		{
			computed_solution[RCI_count] = 0.0;
		}
		
		/*---------------------------------------------------------------------------
		/* Check the correctness and consistency of the newly set parameters
		/*---------------------------------------------------------------------------*/
		//dfgmres_check(&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp); 
		dcg_check(&ivar,computed_solution,Jt_e,&RCI_request,ipar,dpar,tmp);
		if (RCI_request!=0) goto FAILED;

		/*---------------------------------------------------------------------------
		/* Compute the solution by RCI (P)FGMRES solver with preconditioning
		/* Reverse Communication starts here
		/*---------------------------------------------------------------------------*/
		//ONE:  dfgmres(&ivar, computed_solution, Jt_e, &RCI_request, ipar, dpar, tmp);
		ONE:  dcg(&ivar,computed_solution,Jt_e,&RCI_request,ipar,dpar,tmp);

		
		if(RCI_request==0) goto COMPLETE;

		/*---------------------------------------------------------------------------
		/* If RCI_request=1, then compute the vector A*tmp[ipar[21]-1]
		/* and put the result in vector tmp[ipar[22]-1]	
		/*-----------------------DEPRECATED ROUTINE --------------------------------*/
		if (RCI_request==1)
		{
			//mkl_dcsrsymv(&tr, &n, a, ia, ja, tmp, &tmp[n]);
			mkl_dcsrgemv(&cvar, &ivar, acsr, ia, ja, tmp, &tmp[ivar]);
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
			//dfgmres_get(&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp, &itercount);

			/* Compute the current true residual via MKL (Sparse) BLAS routines */
			mkl_dcsrgemv(&cvar, &ivar, acsr, ia, ja, Jt_e, residual); // A*x for new solution x
			dvar=-1.0E0;
			RCI_count=1;
			//daxpy(&ivar, &dvar, Jt_e, &RCI_count, residual, &RCI_count);  // Ax - A*x_solution
			for(j = 0 ; j<ivar ; j++)
				residual[j] = Jt_e[j] - residual[j];
			

			if(ipar[10] == 0)   // non preconditioned system
			{
				dvar=dnrm2(&ivar,residual,&RCI_count);
				relres_nrm = dvar/rhs_nrm;
				
			}
			else if(ipar[10] == 1)  //preconditioned system
			{
				double *prec_relres = new double[num_cols]();

				//prec_solve(A,D,MSC,Numeric_D,Numeric_MSC,lcsr,il,jl,residual,prec_relres);
				prec_solve(A,D,G,Numeric_D,Numeric_G,residual,prec_relres);
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
			//prec_solve(A,D,MSC,Numeric_D,Numeric_MSC,lcsr,il,jl,&tmp[2*ivar],&tmp[3*ivar]);
			prec_solve(A,D,G,Numeric_D,Numeric_G,&tmp[2*ivar],&tmp[3*ivar]);
			goto ONE;
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
		//COMPLETE:  dfgmres_get(&ivar, computed_solution, Jt_e, &RCI_request, ipar, dpar, tmp, &itercount);
		COMPLETE : dcg_get(&ivar,computed_solution,Jt_e,&RCI_request,ipar,dpar,tmp,&itercount);
		//cout << "The system has been solved  in " << itercount << " iterations!\n";
		*total_iters = itercount;

		//store the solution into delta 
		RCI_count = 1;
		for(int q = 0; q < num_cols; q++) delta[q] = computed_solution[q];
	
		MKL_Free_Buffers();


		//  Free the numeric factorization.
	  	umfpack_di_free_numeric ( &Numeric_D );
	  	umfpack_di_free_numeric ( &Numeric_G );

		
		
		delete [] D->p;  delete [] D->i;delete [] D->x; 
		delete [] A->p; delete [] A->i; delete [] A->x; 
		delete [] G->p; delete [] G->i; delete [] G->x; 
		delete A; delete D; delete G; 

		delete [] tmp; delete [] ipar; delete [] dpar; 
		delete [] computed_solution; delete [] residual;
		delete [] acsr; delete [] ia; delete [] ja;

		return ;

		FAILED: cout << "The solver has returned the ERROR code " << RCI_request << "\n";

		MKL_Free_Buffers();

		return ;

	    NOT_CONVERGE: cout << "The relative residual did not change for successive iterations !\n";

	    return ;
	 }
}
#endif