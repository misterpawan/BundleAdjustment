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
	#define sizeG 2025
	#define size_MKL_IPAR 128
	#define MAX_ITERS 100
	#define RESTARTS 40

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
	 							double *values,double *Jt_e,double *delta,int *total_iters)
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

		//int const num_cols = H.num_cols();
		//int const ncc = H.getNonzeroCount();

		//creating structures of cs_diPARSE
		cs_di* A = new cs_di;cs_di* D = new cs_di;cs_di* G = new cs_di;

		//num_cols = num_cols+1;         
		//cout << "\nNo of non zeros = "<< ncc << "\n";
		//cout << "\nNo of columns = "<< num_cols << "\n";

		//size of the D matrix
		sizeD = num_cols - sizeG; 

		A->nzmax = ncc; A->m = num_cols;A->n = num_cols;A->nz = -1;


		//Allocate space for the coefficient matrix
		A->x = new double[ncc](); A->p = new int[num_cols+1](); A->i = new int[ncc]();

		//A->p = colStarts;
		//A->i = rowIdxs;
		//A->x = values;

		for(int q = 0; q < num_cols; q++) A->p[q] = colStarts[q];
		for(int q = 0; q < ncc; q++) A->i[q] = rowIdxs[q];
		for(int q = 0; q < ncc; q++) A->x[q] = values[q];
		//dcopy(&num_cols, values, &count, A->x, &count); 
		A->p[num_cols] = ncc;

		//cout << "\n values[ncc-1] = "<< values[ncc-1] << "\t\t A->x[ncc-1] = "<<A->x[ncc-1]<<"\n";


		/***Domain Decomposition***/

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

		//cout << "\nFilling non zeros complete!!\n";

		//cout << "\n ok : "<< ok << "\n";
	/*
		for(int k = 0; k < 1; k++)
		{
			for(int l = MSC->p[k]; l<20; l++)
				printf("\nMSC->i[%d] = %d\t\tMSC->x[%d] = %10.9f",l,MSC->i[l],l,MSC->x[l]);
		}
	*/
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
		//int restart = 20;  //DFGMRES restarts
		MKL_INT* ipar = new MKL_INT[size_MKL_IPAR]();
		//ipar[14] = 150;  //non restarted iterations

		//cout << "\n tmp size : "<< num_cols*(2*ipar[14]+1)+ipar[14]*((ipar[14]+9)/2+1) << "\n";

		double* dpar = new double[size_MKL_IPAR](); 
		
		double* tmp = new double[num_cols*(2*RESTARTS+1)+(RESTARTS*(RESTARTS+9))/2+1]();
		//double expected_solution[num_cols];
		double* rhs = new double[num_cols]();
		double* computed_solution = new double[num_cols]();
		double* residual = new double[num_cols]();   
		double nrm2,rhs_nrm,relres_nrm,dvar,relres_prev,prec_rhs_nrm,prec_relres_nrm;
		double *prec_rhs = new double[num_cols]();
		double tol = 1.0e-02;
		

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

	      //converting COO to CSR
	    mkl_dcsrcsc(job,&ivar,acsr,ja,ia,A->x,A->i,A->p,&info);
	    //cout << "\n Conversion info A : "<< info << "\n";


	    // Testing A solve with LU and comparing with MATLAB
		//test_A_solve(A);
		//test_D_solve(D);
		//test_MSC_solve(MSC);
		// test_matvec_multiply(A);

	    /***Testing the preconditioner solve with MATLAB***/
	/*    double *matvec = new double[A->n];
	    double *randvec = new double[A->n];
	    string test_filename = "test/rand_vec.txt";

	    r8vec_data_read ( test_filename, A->n, randvec);

		mkl_dcsrgemv(&cvar, &ivar, acsr, ia, ja, randvec, matvec);

		ofstream outfile("test/rand_matvec_cpp.txt");

	  	if(outfile.is_open())
	  	{
	  		for(int k = 0; k < A->n; k++)
	  			outfile << matvec[k] << "\n";
	  	}
	  	outfile.close();

	    test_prec_solve(A,D,MSC,Numeric_D,Numeric_MSC,lcsr,il,jl);
	    delete [] matvec; delete [] randvec;
	*/
		/*---------------------------------------------------------------------------
		/* Save the right-hand side in vector rhs for future use
		/*---------------------------------------------------------------------------*/
		RCI_count=1;
		/** extracting the norm of the rhs for computing rel res norm**/
		rhs_nrm = dnrm2(&ivar,Jt_e,&RCI_count); 
		//cout << "\n rhs_nrm : " << rhs_nrm << "\n";
		// Jt_e vector is not altered
		//rhs is used for residual calculations
		//dcopy(&ivar, Jt_e, &RCI_count, rhs, &RCI_count);   
		for(int q = 0; q < num_cols; q++) rhs[q] = Jt_e[q];

		// PRECONDITIONED RHS
		prec_solve(A,D,G,Numeric_D,Numeric_G,Jt_e,prec_rhs);
		//norm of the preconditioned rhs
		prec_rhs_nrm = dnrm2(&ivar,prec_rhs,&RCI_count);
		delete [] prec_rhs; 
		//cout << "\n prec rhs norm : "<< prec_rhs_nrm << "\n";

		
		/*---------------------------------------------------------------------------
		/* Initialize the solver
		/*---------------------------------------------------------------------------*/
		dfgmres_init(&ivar, computed_solution, Jt_e, &RCI_request, ipar, dpar, tmp); 
		if (RCI_request!=0) goto FAILED;

		
		ipar[7] = 1;
		ipar[4] = MAX_ITERS;  // Max Iterations
		ipar[10] = 1;  //Preconditioner used
		ipar[14] = RESTARTS; //internal iterations
		
		dpar[0] = tol; //Relative Tolerance

		/*---------------------------------------------------------------------------
		/* Initialize the initial guess
		/*---------------------------------------------------------------------------*/
		/*for(RCI_count=0; RCI_count<num_cols; RCI_count++)
		{
			computed_solution[RCI_count]=0.0;
		}*/
		//if(ipar[10] == 1) computed_solution[0]=1000.0;

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
		//cout << "\n dfgmres RCI_request : "<<RCI_request << "\n";

		
		if(RCI_request==0) goto COMPLETE;

		/*---------------------------------------------------------------------------
		/* If RCI_request=1, then compute the vector A*tmp[ipar[21]-1]
		/* and put the result in vector tmp[ipar[22]-1]	
		/*------------------DEPRECATED ROUTINE (FIND ANOTHER )-------------------------*/
		if (RCI_request==1)
		{
			mkl_dcsrgemv(&cvar, &ivar, acsr, ia, ja, &tmp[ipar[21]-1], &tmp[ipar[22]-1]);
		/*	
			ofstream outfile("RCI_1_matvec_cpp.txt");

		  	if(outfile.is_open())
		  	{
		  		for(int k = 0; k < A->n; k++)
		  			outfile << tmp[(ipar[21]-1)+k] << "\n";
		  	}
		  	outfile.close();

		  	char comp_done ;
		  	printf("\nMATLAB  matvec computation done ? :");
		  	scanf(" %c[^\n]",&comp_done);
		  	cout << "\n";
		  	string matvec_filename = "matvec.txt"; 
		  	
		  	if(comp_done == 'y' || comp_done == 'Y')
		  		r8vec_data_read ( matvec_filename, A->n, &tmp[ipar[22]-1]);
		*/
			goto ONE;
		}

		/*---------------------------------------------------------------------------
		/* If RCI_request=2, then do the user-defined stopping test
		/* The residual stopping test for the computed solution is performed here*/
		if (RCI_request==2)
		{
			/* Request to the dfgmres_get routine to put the solution into b[N] via ipar[12]*/
			ipar[12]=1;
			
			//for(int kl = 0; kl < 10; kl++)
			//	printf("\n comp_sol[%d] = %10.9f",kl, rhs[kl]);

			/* Get the current FGMRES solution in the vector rhs[N] */
			dfgmres_get(&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp, &itercount);

			

			/* Compute the current true residual via MKL (Sparse) BLAS routines */
			mkl_dcsrgemv(&cvar, &ivar, acsr, ia, ja, rhs, residual); // A*x for new solution x
			dvar=-1.0E0;
			RCI_count=1;
			daxpy(&ivar, &dvar, Jt_e, &RCI_count, residual, &RCI_count);  // Ax - A*x_solution
			

			//cout << "\n Iteration : " << itercount;
			//for(int k = 0; k < 10; k++)
			//		printf("\nresildual[%d] = %10.9f\n",k,residual[k]);

			if(ipar[10] == 0)   // non preconditioned system
			{
				dvar=dnrm2(&ivar,residual,&RCI_count);
				relres_nrm = dvar/rhs_nrm;
				//relres_nrm = dvar/prec_rhs_nrm;
				//printf("\nresidual norm non prec = %10.9f\n",dvar);
				
			}
			else if(ipar[10] == 1)  //preconditioned system
			{
				double *prec_relres = new double[num_cols]();
				//dvar=dnrm2(&ivar,residual,&RCI_count);
				//printf("\nresidual norm with prec = %10.9f\n",dvar);

				prec_solve(A,D,G,Numeric_D,Numeric_G,residual,prec_relres);
				prec_relres_nrm = dnrm2(&ivar,prec_relres,&RCI_count); 
				delete [] prec_relres;
				//cout << "\n prec relres norm : " << prec_relres_nrm << "\n";
				//printf("\nPrec relres norm : %10.9f",prec_relres_nrm);
				relres_nrm = prec_relres_nrm/prec_rhs_nrm;

			}

			//cout << "\n relres_nrm : " << relres_nrm << "\n";
			//printf("\nRelres norm = %10.9f\n",relres_nrm);

			if (relres_nrm<=tol) goto COMPLETE;   //taking tolerance as 1e-04

			else goto ONE;
			
		}

		/*---------------------------------------------------------------------------
		/* If RCI_request=3, then apply the preconditioner on the vector
		/* tmp[ipar[21]-1] and put the result in vector tmp[ipar[22]-1]*/
		
		if (RCI_request==3)
		{
			//cout << "\n Prec solve ..." << "\n";
			prec_solve(A,D,G,Numeric_D,Numeric_G,&tmp[(ipar[21]-1)],&tmp[(ipar[22]-1)]);
		/*	
			ofstream outfile("RCI_3_prec_solve_cpp.txt");

		  	if(outfile.is_open())
		  	{
		  		for(int k = 0; k < A->n; k++)
		  			outfile << tmp[(ipar[21]-1)+k] << "\n";
		  	}
		  	outfile.close();

		  	char comp_done ;
		  	printf("\nMATLAB prec solve computation done ? :");
		  	scanf(" %c[^\n]",&comp_done);
		  	cout << "\n";
		  	string precsolve_filename = "prec_solve.txt"; 
		  	
		  	if(comp_done == 'y' || comp_done == 'Y')
		  		r8vec_data_read ( precsolve_filename, A->n, &tmp[ipar[22]-1]);
		*/
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
		*total_iters = itercount;
	//	cout << "\n RCI_request : "<< RCI_request << "\n";
	/*
		printf("\nThe following solution has been obtained: \n");
		for (RCI_count=0;RCI_count<10;RCI_count++)                //PRINTING ONLY THE FIRST 10 MEMBERS
		{
			printf("computed_solution[%d]=%f\n",RCI_count,computed_solution[RCI_count]);
		}
	*/

		//store the solution into delta 
		RCI_count = 1;
		//dcopy(&ivar, computed_solution, &RCI_count, delta, &RCI_count); 
		for(int q = 0; q < num_cols; q++) delta[q] = computed_solution[q];
	/*	for (RCI_count=0;RCI_count<10;RCI_count++)                //PRINTING ONLY THE FIRST 10 MEMBERS
		{
			printf("delta[%d]=%f\n",RCI_count,delta[RCI_count]);
		}
	*/
		MKL_Free_Buffers();


		//  Free the numeric factorization.
	  	umfpack_di_free_numeric ( &Numeric_D );
	  	umfpack_di_free_numeric ( &Numeric_G );

		//delete [] Jt_e; 
		delete [] A->p; delete [] A->i; delete [] A->x; 
		delete [] D->p;  delete [] D->i;delete [] D->x; 
		delete D; delete A; 

		delete [] tmp; delete [] ipar; delete [] dpar;
		delete [] rhs; delete [] computed_solution; delete [] residual;
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