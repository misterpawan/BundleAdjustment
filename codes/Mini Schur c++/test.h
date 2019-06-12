/* 
 This file contains test that are to be computed with MATLAB for accuracy
*/

void test_A_solve(cs_di_sparse* A)
{
	double* b = new double[A->n];
	string test_filename = "test/rand_rhs.txt";
	string solve_filename = "test/rand_sol.txt";
	void* Numeric;
	void* Symbolic;
	double* rhs = new double[A->n];
	int status;
	double *null = (double *) NULL; 
	double *rhs_mat = new double[A->n];
	double rhs_nrm, rhs_mat_nrm;
	int incx = 1;
	double diff;
	int zvar = A->n;
	double dvar = -1.0E0;
	int p = 1;

	r8vec_data_read ( test_filename, A->n, b);
	r8vec_data_read ( solve_filename, A->n, rhs_mat); // MATLAB solution

	status = umfpack_di_symbolic ( A->m, A->n, A->p, A->i, A->x, &Symbolic, null, null );
	//cout << "\n Symbolic status :" << status << "\n";

	status = umfpack_di_numeric ( A->p, A->i, A->x, Symbolic, &Numeric, null, null );
	//cout << "\n Numeric status :" << status << "\n";

	umfpack_di_free_symbolic ( &Symbolic );

	status = umfpack_di_solve ( UMFPACK_A, A->p, A->i, A->x, rhs, b, Numeric, null, null );
  	//cout << "\n Solve status :" << status << "\n";

  	umfpack_di_free_numeric ( &Numeric );

  	for(int i = 0; i < 20; i++)
		printf("\nrhs[%d] = %10.9f \t\t rhs_mat[%d] = %10.9f",i,rhs[i],i,rhs_mat[i]);

/*
  	ofstream outfile("rand_sol_cpp.txt");

  	if(outfile.is_open())
  	{
  		for(int k = 0; k < A->n; k++)
  			outfile << rhs[k] << "\n";
  	}
  	outfile.close();
*/
  	daxpy(&zvar, &dvar, rhs, &p, rhs_mat, &p);

  	//rhs_nrm = dnrm2(&(A->n),rhs,&incx); // cout << "\n norm rhs : " << rhs_nrm << "\n";
  	//printf("\nnorm_rhs = %10.9f",rhs_nrm);
  	//rhs_mat_nrm = dnrm2(&(A->n),rhs_mat,&incx);// cout<< "\n norm_rhs_mat : "<< rhs_mat_nrm << "\n";
  	//printf("\nnorm_rhs_mat = %10.9f",rhs_mat_nrm);

  	//diff = abs(rhs_mat_nrm - rhs_nrm);
  	diff = dnrm2(&zvar,rhs_mat,&incx);

  	cout << "\n The difference in norms is " << diff << "\n";

  	delete [] b; delete [] rhs; delete [] rhs_mat;

	return;
}


void test_D_solve(cs_di_sparse* D)
{
	double* b = new double[D->n];
	string test_filename = "test/rand_rhs_D.txt";
	string solve_filename = "test/rand_sol_D.txt";
	void* Numeric;
	void* Symbolic;
	double* rhs = new double[D->n];
	int status;
	double *null = (double *) NULL; 
	double *rhs_mat = new double[D->n];
	double rhs_nrm, rhs_mat_nrm;
	int incx = 1;
	double diff;
	int zvar = D->n;
	double dvar = -1.0E0;
	int p = 1;

	r8vec_data_read ( test_filename, D->n, b);
	r8vec_data_read ( solve_filename, D->n, rhs_mat); // MATLAB solution

	status = umfpack_di_symbolic ( D->m, D->n, D->p, D->i, D->x, &Symbolic, null, null );
	//cout << "\n Symbolic status :" << status << "\n";

	status = umfpack_di_numeric ( D->p, D->i, D->x, Symbolic, &Numeric, null, null );
	//cout << "\n Numeric status :" << status << "\n";

	umfpack_di_free_symbolic ( &Symbolic );

	status = umfpack_di_solve ( UMFPACK_A, D->p, D->i, D->x, rhs, b, Numeric, null, null );
  	//cout << "\n Solve status :" << status << "\n";

  	umfpack_di_free_numeric ( &Numeric );

  	for(int i = 0; i < 20; i++)
		printf("\nrhs[%d] = %10.9f \t\t rhs_mat[%d] = %10.9f",i,rhs[i],i,rhs_mat[i]);
/*
  	ofstream outfile("rand_sol_cpp.txt");

  	if(outfile.is_open())
  	{
  		for(int k = 0; k < D->n; k++)
  			outfile << rhs[k] << "\n";
  	}
  	outfile.close();
*/
  	daxpy(&zvar, &dvar, rhs, &p, rhs_mat, &p);

  	//rhs_nrm = dnrm2(&(A->n),rhs,&incx); // cout << "\n norm rhs : " << rhs_nrm << "\n";
  	//printf("\nnorm_rhs = %10.9f",rhs_nrm);
  	//rhs_mat_nrm = dnrm2(&(A->n),rhs_mat,&incx);// cout<< "\n norm_rhs_mat : "<< rhs_mat_nrm << "\n";
  	//printf("\nnorm_rhs_mat = %10.9f",rhs_mat_nrm);

  	//diff = abs(rhs_mat_nrm - rhs_nrm);
  	diff = dnrm2(&zvar,rhs_mat,&incx);

  	cout << "\n The difference in norms for D is " << diff << "\n";

  	delete [] b; delete [] rhs; delete [] rhs_mat;

	return;
}

void test_MSC_solve(cs_di_sparse* MSC)
{
	double* b = new double[MSC->n];
	string test_filename = "test/rand_rhs_MSC_20.txt";
	string solve_filename = "test/rand_sol_MSC_20.txt";
	void* Numeric;
	void* Symbolic;
	double* rhs = new double[MSC->n];
	int status;
	double *null = (double *) NULL; 
	double *rhs_mat = new double[MSC->n];
	double rhs_nrm, rhs_mat_nrm;
	int incx = 1;
	double diff;
	int zvar = MSC->n;
	double dvar = -1.0E0;
	int p = 1;

	r8vec_data_read ( test_filename, MSC->n, b);
	r8vec_data_read ( solve_filename, MSC->n, rhs_mat); // MATLAB solution

	status = umfpack_di_symbolic ( MSC->m, MSC->n, MSC->p, MSC->i, MSC->x, &Symbolic, null, null );
	//cout << "\n Symbolic status :" << status << "\n";

	status = umfpack_di_numeric ( MSC->p, MSC->i, MSC->x, Symbolic, &Numeric, null, null );
	//cout << "\n Numeric status :" << status << "\n";

	umfpack_di_free_symbolic ( &Symbolic );

	status = umfpack_di_solve ( UMFPACK_A, MSC->p, MSC->i, MSC->x, rhs, b, Numeric, null, null );
  	//cout << "\n Solve status :" << status << "\n";

  	umfpack_di_free_numeric ( &Numeric );
/*
  	ofstream outfile("rand_sol_cpp.txt");

  	if(outfile.is_open())
  	{
  		for(int k = 0; k < MSC->n; k++)
  			outfile << rhs[k] << "\n";
  	}
  	outfile.close();
*/
  	for(int i = 0; i < 20; i++)
		printf("\nprec_MSC_sol[%d] = %10.9f \t\t prec_mat_MSC_sol[%d] = %10.9f",i,rhs[i],i,rhs_mat[i]);

  	daxpy(&zvar, &dvar, rhs, &p, rhs_mat, &p);

  	//rhs_nrm = dnrm2(&(A->n),rhs,&incx); // cout << "\n norm rhs : " << rhs_nrm << "\n";
  	//printf("\nnorm_rhs = %10.9f",rhs_nrm);
  	//rhs_mat_nrm = dnrm2(&(A->n),rhs_mat,&incx);// cout<< "\n norm_rhs_mat : "<< rhs_mat_nrm << "\n";
  	//printf("\nnorm_rhs_mat = %10.9f",rhs_mat_nrm);

  	//diff = abs(rhs_mat_nrm - rhs_nrm);
  	diff = dnrm2(&zvar,rhs_mat,&incx);

  	//cout << "\n The difference in norms for MSC is " << diff << "\n";
  	printf("\n The difference in norms for MSC is %g\n",diff);

  	delete [] b; delete [] rhs; delete [] rhs_mat;

	return;
}


void test_matvec_multiply(cs_di_sparse* A)
{
	double* b = new double[A->n];
	string test_filename = "test/rand_vec.txt";
	string solve_filename = "test/rand_product.txt";
	double* vec = new double[A->n];
	int status;
	double *vec_mat = new double[A->n];
	double vec_nrm, vec_mat_nrm;
	double diff;
	int incx = 1;
	char cvar = 'N';
	int ivar = A->n;
	int nnz = A->nzmax;
	double dvar = -1.0E0;
	int p = 1;

	MKL_INT job[6] = {1,1,0,0,0,1};
    double *acsr = (double*) malloc (nnz* sizeof(double));
    MKL_INT *ja = (MKL_INT*) malloc (nnz*sizeof(MKL_INT));
    MKL_INT *ia = (MKL_INT*) malloc ((ivar+1)*sizeof(MKL_INT));
    MKL_INT info;

    //converting COO to CSR
    mkl_dcsrcsc(job,&ivar,acsr,ja,ia,A->x,A->i,A->p,&info);
    //printf("\n Conversion info : %d\n",info);

	r8vec_data_read ( test_filename, A->n, b);
	r8vec_data_read ( solve_filename, A->n, vec_mat); // MATLAB solution
	
	//A->p[ivar] = A->p[ivar]+1;
	//int b_nrm = dnrm2(&ivar,b,&incx); cout << "\n norm b : " << b_nrm << "\n";
	

	mkl_dcsrgemv(&cvar, &ivar, acsr, ia, ja, b, vec);

	for(int i = 0; i < 20; i++)
		printf("\nvec[%d] = %10.9f \t\t vec_mat[%d] = %10.9f",i,vec[i],i,vec_mat[i]);

  	//vec_nrm = dnrm2(&ivar,vec,&incx); cout << "\n norm rhs : " << vec_nrm << "\n";
  	//vec_mat_nrm = dnrm2(&ivar,vec_mat,&incx); cout<< "\n norm_vec_mat : "<< vec_mat_nrm << "\n";

	daxpy(&ivar, &dvar, vec_mat, &p, vec, &p);
  	//diff = abs(vec_nrm - vec_mat_nrm);
	diff = dnrm2(&ivar,vec,&p);

  	cout << "\n The difference in norms for matvec multiply  is " << diff << "\n";

 	delete [] acsr; delete [] ia; delete [] ja;
  	delete [] b; delete [] vec; delete [] vec_mat;


	return;
}


void test_prec_solve(cs_di_sparse* A, cs_di_sparse* D,cs_di_sparse* MSC,void *Numeric_D,void *Numeric_MSC,double *lcsr,int *il,int *jl)
{
	double* y1 = new double[D->n];
	double* y2 = new double[MSC->n];
	double* z1 = new double[D->n];
	double* z2 = new double[MSC->n];
	double* Lz1 = new double[MSC->n];
	int prec_solve_status;
	double* null = (double*)NULL;
	int zvar = MSC->n ; //no of rows of L
	int kk;
	char cvar = 'N';

	double *rhs_prec = new double[A->n];
	double *prec_sol = new double[A->n];
	double *prec_mat_sol = new double[A->n];
	string test_filename = "test/rand_rhs_prec.txt";
	string solve_filename = "test/rand_sol_prec.txt";
	int ivar = A->n;
	int p = 1;
	double diff;
	double dvar = -1.0E0;

	r8vec_data_read ( test_filename, A->n, rhs_prec);
	r8vec_data_read ( solve_filename, A->n, prec_mat_sol); // MATLAB solution

	for(kk = 0; kk<ivar; kk++)
	{
		if(kk < (D->n)) y1[kk] = rhs_prec[kk];
		else y2[kk-(D->n)] = rhs_prec[kk];  //splitting vector into y1,y2
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
		if(kk < (D->n)) prec_sol[kk] = z1[kk];
		else prec_sol[kk] = z2[kk - (D->n)];
	}

	for(int i = 0; i < 20; i++)
		printf("\nprec_sol[%d] = %10.9f \t\t prec_mat_sol[%d] = %10.9f",i,prec_sol[i],i,prec_mat_sol[i]);

	daxpy(&ivar, &dvar, prec_sol, &p, prec_mat_sol, &p);
	diff = dnrm2(&ivar,prec_mat_sol,&p);

  	cout << "\n The difference in norms for prec solve  is " << diff << "\n";

	delete [] y1; delete [] y2; delete [] z1; delete [] z2; delete [] Lz1;
	delete [] rhs_prec; delete [] prec_mat_sol; delete [] prec_sol;

	return;
}
