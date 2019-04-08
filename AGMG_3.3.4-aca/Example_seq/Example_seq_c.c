/*
! COPYRIGHT (c) 2012-2018 Yvan Notay - ULB
!
! This file is part of AGMG software package
! Release 3.3.4-aca built on "Oct 22 2018" by Yvan Notay
!
! ALL USAGE OF AGMG IS SUBJECT TO LICENSE. PLEASE REFER TO THE FILE "LICENSE".
! IF YOU OBTAINED A COPY OF THIS SOFTWARE WITHOUT THIS FILE,
! PLEASE CONTACT info@agmg.eu
!
! In particular, if you have a free academic license:
!
! (1) You must be a member of an educational, academic or research institution.
!     The license agreement automatically terminates once you no longer fulfill
!     this requirement.
!
! (2) You are obliged to cite AGMG in any publication or report as:
!     "Yvan Notay, AGMG software and documentation;
!      see http://agmg.eu".
!
! (3) You may not make available to others the software in any form, either
!     as source or as a precompiled object.
!
! (4) You may not use AGMG for the benefit of any third party or for any
!     commercial purposes. Note that this excludes the use within the
!     framework of a contract with an industrial partner.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! See the accompanying userguide for more details on how to use the software,
! and the README file for installation instructions.
!
! See the Web pages <http://agmg.eu> for
! release information and possible upgrade.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DICLAIMER:
!    AGMG is provided on an "AS IS" basis, without any explicit or implied
!    WARRANTY; see the file "LICENSE" for more details.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   If you use AGMG for research, please observe that your work benefits
!   our past research efforts that allowed the development of AGMG.
!   Hence, even if you do not see it directly, the results obtained thanks
!   to the use of AGMG depend on the results in publications [1-3] below,
!   where the main algorithms used in AGMG are presented and justified.
!   It is then a normal duty to cite these publications (besides citing
!   AGMG itself) in any scientific work depending on the usage of AGMG,
!   as you would do with any former research result you are using.
!
! [1] Y. Notay, An aggregation-based algebraic multigrid method,
!    Electronic Transactions on Numerical Analysis, vol. 37, pp. 123-146, 2010
!
! [2] A. Napov and Y. Notay, An algebraic multigrid method with guaranteed
!    convergence rate, SIAM J. Sci. Comput., vol. 34, pp. A1079-A1109, 2012.
!
! [3] Y. Notay, Aggregation-based algebraic multigrid for convection-diffusion
!    equations, SIAM J. Sci. Comput., vol. 34, pp. A2288-A2316, 2012.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
      #include <stdlib.h>
      #include <stdio.h>
      void dagmg_(int*,double*,int*,int*,double*,double*,int*,int*,int*,int*,double*);
      void uni2d(int,double*,double*,int*,int*);

      int main(void) {
/*
   Solves the discrete Laplacian on the unit square by simple call to agmg.
   The right-hand-side is such that the exact solution is the vector of all 1.
*/
        double *a,*f,*x;
        int n,*ja,*ia;
        int zero=0,one=1; 

/*      set inverse of the mesh size (feel free to change) */
        int nhinv=500;
/*      maximal number of iterations */
        int iter=50;
/*      tolerance on relative residual norm */
        double tol=1.e-6;
/*      unit number for output messages: 6 => standard output */
        int iprint=6;

/*      generate the matrix in required format (CSR) */

/*        first allocate the vectors with correct size */
            n=(nhinv-1)*(nhinv-1);
            ia=malloc((n + 1) * sizeof(int));
            ja=malloc(5*n * sizeof(int));
            a=malloc(5*n * sizeof(double));
            f=malloc(n * sizeof(double));
            x=malloc(n * sizeof(double));
/*        next call subroutine to set entries */
            uni2d(nhinv-1,f,a,ja,ia);

/*      call agmg
         argument 5 (ijob)  is 0 because we want a complete solve
         argument 7 (nrest) is 1 because we want to use flexible CG
                             (the matrix is symmetric positive definite) */
        dagmg_(&n,a,ja,ia,f,x,&zero,&iprint,&one,&iter,&tol);

/*      uncomment the following lines to write solution on disk for checking */
/*
        int i;
        FILE *F; 
        F = fopen("sol_c.out","w");
        for (i=0; i<n; i++)
           fprintf(F,"%021.14E\n", f[i]);
        fclose(F); 
*/
      }
/*----------------------------------------------------------------------*/
      void uni2d(int m, double *f, double *a, int *ja, int *ia){
/*
  Fill a matrix in CSR format corresponding to a constant coefficient
  five-point stencil on a square grid */

      int k,km,l,i,j;
      double const zero=0.0,cx=-1.0,cy=-1.0, cd=4.0;
      
      f--; /* so that f[k] (ex f[1]) means f[k-1] (ex f[0]) */
      k=0;
      l=0;
      ia[0]=1;
      for(i=1;i<=m;i++)
        for(j=1;j<=m;j++){
          k++;
          a[l]=cd;
          ja[l]=k;
          l++;
          f[k]=zero;
          if(j < m){
             a[l]=cx;
             ja[l]=k+1;
             l++;
          } else
             f[k]=f[k]-cx;
          if(i < m){
             a[l]=cy;
             ja[l]=k+m;
             l++;
          } else
             f[k]=f[k]-cy;
          if(j > 1){
             a[l]=cx;
             ja[l]=k-1;
             l++;
          } else
             f[k]=f[k]-cx;
          if(i > 1){
             a[l]=cy;
             ja[l]=k-m;
             l++;
          } else
             f[k]=f[k]-cy;
          ia[k]=l+1;
        }
      }
