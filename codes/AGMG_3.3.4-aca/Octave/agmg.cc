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
#include <octave/oct.h>
#include <octave/f77-fcn.h>

/* The two following lines may be commented with most recent Ocatve releases
   (avoids warning messages about deprecated features) */
#define isempty is_empty
#define issparse is_sparse_type

#define PROCOUT(ff)				\
                retval(0) = xm;                  \
                if (iter < 0) {                   \
                  retval(1) = 1;                   \
                  iter = -iter;                     \
                }                                    \
                else {                                \
                  retval(1) = 0;                       \
                }                                       \
                retval(2) = tol;                         \
                retval(3) = iter;                         \
                if (nargout > 4) {                         \
                  nt = iter+1;                              \
                  if (nt > n) nt=n;                          \
                  Matrix resvec (nt,1);                       \
                  rv = resvec.fortran_vec ();                 \
                  for (i=0 ; i < iter+1; i++) {rv[i]=ff;}     \
                  retval(4) = resvec;                         \
                }                                             \
                return retval                                

extern "C"
{
  F77_RET_T
  F77_FUNC (dagmg, FORTSUB)
    (const int&, double *, int *, int *, double *, double *,
     const int&, const int&, const int&, const int&, const double&);
}
extern "C"
{
  F77_RET_T
  F77_FUNC (zagmg, FORTSUB)
    (const int&, Complex *, int *, int *, Complex *, Complex *,
     const int&, const int&, const int&, const int&, const double&);
}

DEFUN_DLD (agmg, args, nargout,
	     "-*- texinfo -*-\n\
@deftypefn {Function} {@var{x} = } agmg(@var{A}, @var{b})\n\
@deftypefnx {Function} {@var{x} = } agmg(@var{A}, @var{b}, @var{restart}, @var{tol}, @var{maxit}, @var{verbose}, @var{x0}, @var{ijob})\n\
@deftypefnx {Function} {[@var{x}, @var{flag}, @var{relres}, @var{iter}, @var{resvec}] = } agmg(@var{A}, @var{b}), ...)\n\
\n\
AGMG: solve the linear system of equations '@var{A} * @var{X}= @var{B}'\n\
by means of an Aggregation-based algebraic Multigrid iterative Method\n\
accelerated either by the Conjugate Gradient method (CG) or by\n\
the Generalized Conjugate Residual method (GCR, a variant of GMRES).\n\
\n\
The N-by-N coefficient matrix @var{A} must be square and sparse,\n\
and the right hand side column vector @var{B} must have length N.\n\
The matrix @var{A} may be real or complex; for complex matrices, however,\n\
AGMG is tentative only. If @var{A} is real, @var{A} has to be real too.\n\
\n\
The input arguments are\n\
\n\
@itemize\n\
@item @var{A}: square sparse (real or complex) matrix.@*\n\
\n\
@item @var{b}: right hand side vector.\n\
\n\
@item @var{restart}:\n\
@* If @var{restart}==1, use CG and perform some simplifications based on\n\
            the assumption that the coefficient matrix A is symmetric;\n\
            should be used only when A is symmetric and positive definite.\n\
@* If @var{restart}>=2, use GCR restarted each @var{restart} iterations.\n\
@* If @var{restart}==0 or [], then AGMG uses the default, which is\n\
                   GCR restarted each 10 iterations.@*\n\
\n\
@item @var{tol}:specifies the tolerance of the method.\n\
@*      If @var{tol} is [], then AGMG uses the default, 1e-6.@*\n\
\n\
@item @var{maxit}:specifies the maximum number of iterations.\n\
@*      If @var{maxit} is [], then AGMG uses the default, 100.\n\
\n\
@item @var{verbose}:\n\
@*  If @var{verbose}==1, information is displayed on the solution process.\n\
@*  If @var{verbose}==0 or [], only error messages are displayed (default).\n\
@*  If @var{verbose}==-1, warning messages about insufficient convergence\n\
                       (in the allowed maximum number of iterations)\n\
                       are further suppressed; i.e., AGMG works silently.\n\
@*     See the user guide provided with AGMG for a description of the\n\
       verbose output.\n\
\n\
@item @var{x0}: specifies an initial guess.\n\
   @*    If X0 is [], then AGMG uses the default, an all zero vector.\n\
\n\
@item @var{ijob}:\n\
@*    If @var{ijob}==1, perform the setup only (preprocessing: \n\
      prepare all parameters for subsequent solves).\n\
      Then, only @var{A}, @var{restart} and @var{verbose} input parameters are\n\
      significant and the calling may be:\n\
      agmg(@var{A},[],@var{restart},[],[],[],[],1).\n\
      The returned @var{x} is empty and other output parameters are meaningless.\n\
@* @* If @var{ijob}==2, solve only, based on previous setup.\n\
      Then, @var{A} may differ from the matrix supplied for set up (former\n\
      call with @var{ijob}==1), but it means using a preconditioner\n\
      computed for a matrix to solve a system with another matrix,\n\
      which is not recommended in general.\n\
@* @* If @var{ijob}==202, solve only, based on previous setup.\n\
      The system matrix is assumed unchanged and argument @var{A} is not significant.\n\
@* @* If @var{ijob}==3, the returned @var{x} is not the solution\n\
      of the linear system, but the result of the action of the multigrid \n\
      preconditioner on the right hand side @var{B}.\n\
      Then, @var{A}, @var{tol}, @var{maxit},  @var{restart},\n\
      and @var{X0} are not significant;\n\
      the calling may be: agmg([],@var{b},[],[],[],[],[],3)\n\
      Further output parameters (besides @var{x}) are meaningless. \n\
@* @* If @var{ijob}==-1, erases the setup and releases internal memory.\n\
      Other input parameters are not significant and the calling may be\n\
      agmg([],[],[],[],[],[],[],-1).\n\
      The returned @var{x} is empty and other output parameters are meaningless.\n\
@* @* DEFAULT: If @var{ijob}==0 or [], solve the linear system (performing setup and solve),\n\
      and release memory\n\
      (other input & output parameters have their usual meaning).\n\
\n\
       @var{ijob} == 100,101,102: same as, respectively, @var{ijob}==0,1,2,\n\
       but use the TRANSPOSE of the input matrix in @var{A}.\n\
       Hence, agmg(@var{A}, @var{b}, @var{restart}, @var{tol}, @var{maxit},@var{verbose}, @var{x0}, @var{ijob}+100)\n\
       is equivalent to\n\
       agmg(@var{A}', @var{b}, @var{restart}, @var{tol}, @var{maxit},@var{verbose}, @var{x0}, @var{ijob})\n\
       but less memory consuming and often faster.\n\
\n\
      @var{ijob}==2,3,102,202 require that one has previously called AGMG\n\
       with @var{ijob}==1 or @var{ijob}==101\n\
\n\
@end itemize\n\
\n\
The output arguments are\n\
@itemize\n\
\n\
@item @var{x}: solution of the linear system.\n\
@* (Result of preconditioner application in case @var{ijob}==3.)\n\
\n\
@item @var{flag}: convergence flag:\n\
@*    0 if AGMG converged to the desired tolerance;\n\
@*    1 if AGMG did not reached the desired tolerance.\n\
\n\
@item @var{relres}:\n\
     relative residual norm (norm(@var{b}-@var{A}*@var{x0})/norm(@var{b})).\n\
@* (If @var{flag} is 0, then @var{relres} <= @var{tol}.)\n\
@* It may happen that \n\
   @var{relres} > norm(@var{b}-@var{A}*@var{x0})/norm(@var{b}), but\n\
   this may occur only when @var{relres} is\n\
   below the accuracy that can be attained with a backward stable solver;\n\
   then, norm(@var{b}-@var{A}*@var{x0})/norm(@var{b})\n\
   is comparable to what can be obtained with a direct solver.\n\
\n\
@item @var{iter}: number of iterations actually performed.\n\
\n\
@item @var{resvec}: vector of the residual norms at each iteration,\n\
  including norm(@var{b}-@var{A}*@var{x0}).\n\
@end itemize\n\
\n\
Example:\n\
\n\
@* % generate sample matrix:\n\
@* % 5-point discrete Poisson operator on a 50x50 grid\n\
@* m=50;\n\
@* n=m^2;\n\
@* v=-[ones(m-1,1);0]*ones(1,m);\n\
@* fx=v(:);\n\
@* fy=v';\n\
@* fy=fy(:);\n\
@* A=spdiags([fx,fy],[-1 -m],n,n);\n\
@* A=A+A';\n\
@* D=4*ones(n,1);\n\
@* A=(1/m^2)*spdiags(D,0,A);\n\
\n\
   % generate sample right hand side vector:\n\
@* b=ones(n,1);\n\
\n\
   % solve with agmg:\n\
@* x=agmg(A,b,1);         % A is symmetric positive definite\n\
@* x=agmg(A,b);           % Consider A as a general matrix\n\
\n\
   x=agmg(A,b,1,[],[],1); % Verbose output\n\
\n\
   % Example of tolerance below attainable accuracy:\n\
@* [x,flag,relres,iter,resvec]=agmg(A,b,1,1e-20);\n\
@* disp('AGMG relative residual = '),disp(relres)\n\
@* disp('True relative residual = '),disp(norm(b-A*x)/norm(b))\n\
@* % The true residual is above the one reported by AGMG,\n\
@* % but the attained accuracy is similar to the one obtained with\n\
@* % the built-in direct solver ('\' or mldivide):\n\
@* y=mldivide(A,b);\n\
@* disp('Relative residual with direct solver = '),\n\
@* disp(norm(A*y-b)/norm(b))\n\
\n\
\n\
COPYRIGHT (c) 2012-2018 Yvan Notay - ULB\n\
@* This function is part of AGMG software package, release VERSION\n\
@* ALL USAGE OF AGMG IS SUBJECT TO LICENSE\n\
@* Enter agmg() for detailed conditions of use.\n\
@* See the Web page <http://agmg.eu> for\n\
release information, a detailed userguide and possible upgrades.\n\
@end deftypefn")
{
  int i, ijob, ijb, nrest, iter, iprint, nz, irhsreal, ix0;
  double tol;
  int *ja, *ia;
  double *a, *f, *x, *rv;
  double *restartp, *tlp, *maxitp, *verbosep, *ijbbp;
  Complex *ac, *fc, *xc;
  static int preprocessed=0, n, notcpl;
  octave_value_list retval;
  octave_idx_type nr, nc, nn, *iia, *jja, nrrhs, ncrhs, nt;
  Matrix empt (0,1);
  retval(0) = empt;
  retval(1) = 1;
  retval(2) = empt;
  retval(3) = 0;
  retval(4) = empt;
  int nargin = args.length ();
  if (nargin < 2){
  octave_stdout
  << " COPYRIGHT (c) 2012-2018 Yvan Notay - ULB \n"
  << " The function agmg is part of AGMG software package \n"
  << " \n"
  << " ALL USAGE OF AGMG IS SUBJECT TO LICENSE. PLEASE REFER TO THE FILE 'LICENSE'. \n"
  << " IF YOU OBTAINED A COPY OF THIS SOFTWARE WITHOUT THIS FILE, \n"
  << " PLEASE CONTACT info@agmg.eu \n"
  << " \n"
  << " In particular, if you have a free academic license: \n"
  << " \n"
  << " (1) You must be a member of an educational, academic or research institution. \n"
  << "     The license agreement automatically terminates once you no longer fulfill \n"
  << "     this requirement. \n"
  << " \n"
  << " (2) You are obliged to cite AGMG in any publication or report as: \n"
  << "     'Yvan Notay, AGMG software and documentation; \n"
  << "      see http://agmg.eu'. \n"
  << " \n"
  << " (3) You may not make available to others the software in any form,  \n"
  << "     either as source or as a precompiled object. \n"
  << " \n"
  << " (4) You may not use AGMG for the benefit of any third party or for any \n"
  << "     commercial purposes. Note that this excludes the use within the \n"
  << "     framework of a contract with an industrial partner. \n"
  << " \n"
  << " See the Web pages <http://agmg.eu> for \n"
  << "    release information, a detailed userguide and possible upgrades. \n"
  << " \n"
  << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n"
  << " DICLAIMER: \n"
  << "    AGMG is provided on an 'AS IS' basis, without any explicit or implied \n"
  << "    WARRANTY; see the see the file 'LICENSE' for more details. \n"
  << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n"
  << "   If you use AGMG for research, please observe that your work benefits \n"
  << "   our past research efforts that allowed the development of AGMG. \n"
  << "   Hence, even if you do not see it directly, the results obtained thanks \n"
  << "   to the use of AGMG depend on the results in publications [1-3] below,\n"
  << "   where the main algorithms used in AGMG are presented and justified.   \n"
  << "   It is then a normal duty to cite these publications (besides citing \n"
  << "   AGMG itself) in any scientific work depending on the usage of AGMG, \n"
  << "   as you would do with any former research result you are using. \n"
  << " \n"
  << " [1] Y. Notay, An aggregation-based algebraic multigrid method, \n"
  << "    Electronic Transactions on Numerical Analysis, vol. 37, pp. 123-146, 2010 \n"
  << "\n"
  << " [2] A. Napov and Y. Notay, An algebraic multigrid method with guaranteed \n"
  << "    convergence rate, SIAM J. Sci. Comput., vol. 34, pp. A1079-A1109, 2012. \n"
  << " \n"
  << " [3] Y. Notay, Aggregation-based algebraic multigrid for convection-diffusion \n"
  << "    equations, SIAM J. Sci. Comput., vol. 34, pp. A2288-A2316, 2012. \n"
  << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n"
  << " \n"
  << " Error using agmg\n"
  << " agmg requires at least two input arguments\n"
  << " enter 'help agmg' for more information\n";
    return retval;
  }
  nrest=0;
  if (nargin > 2){
    if (! args(2).isempty () ) {
      if ( args(2).is_real_scalar () ) {
        Matrix restart = args(2).matrix_value ();
        restartp = restart.fortran_vec ();
        nrest=restartp[0];
      }
      else
        octave_stdout <<
          " Warning: using default value 0 for argument RESTART\n";
    }
  }
  tol=1e-6;
  if (nargin > 3){
    if (! args(3).isempty () ) {
      if ( args(3).is_real_scalar () ) {
        Matrix tl = args(3).matrix_value ();
        tlp = tl.fortran_vec ();
        tol=tlp[0];
       }
      else
        octave_stdout <<
          " Warning: using default value 1e-6 for argument TOL\n";
   }
  }
  iter=100;
  if (nargin > 4){
    if (! args(4).isempty () ) {
      if ( args(4).is_real_scalar () ) {
        Matrix maxit = args(4).matrix_value ();
        maxitp = maxit.fortran_vec ();
        iter=maxitp[0];
       }
      else
        octave_stdout <<
          " Warning: using default value 100 for argument MAXIT\n";
    }
  }
  iprint=0;
  if (nargin > 5){
    if (! args(5).isempty () ) {
      if ( args(5).is_real_scalar () ) {
        Matrix verbose = args(5).matrix_value ();
        verbosep = verbose.fortran_vec ();
        if (verbosep[0] < 0)
          iprint=-1;
        else if (verbosep[0] > 0)
          iprint=6;
      }
      else
        octave_stdout <<
          " Warning: using default value 0 for argument VERBOSE\n";
    }
  }
  ijob=0;
  if (nargin > 7){
    if (! args(7).isempty () ) {
      if ( args(7).is_real_scalar () ) {
        Matrix ijbb = args(7).matrix_value ();
        ijbbp = ijbb.fortran_vec ();
        ijob=ijbbp[0];
      }
      else
        octave_stdout <<
          " Warning: using default value 0 for argument IJOB\n";
    }
  }
  if (ijob>99 && ijob<103 && ijob!=202) {
    ijb=ijob-100;
    ijob=ijb;
  }
  else {
    ijb=ijob;
    if (ijob>=0 && ijob<=2) ijob=ijob+100;
  }
  if (ijb!=-1 && ijb!=0 && ijb!=1 && ijb!=2 && ijb!=3 && ijob!=202) {
    octave_stdout << " Error: IJOB should be equal to -1, 0, 1, 2, 3, 100, 101, 102 or 202\n";
    return retval;
  }
  else if (ijb>1 && preprocessed==0) {
    octave_stdout << " Error: setup not done: IJOB should be equal to 0, 1, 100 or 101\n";
    return retval;
  }
  else if (ijb == -1) {
    if (preprocessed==0) {
      octave_stdout << " Warning: setup not done: nothing to do for IJOB == -1\n";
    }
    else {
      if (notcpl == 1) {
            F77_XFCN (dagmg, FORTSUB,
              (n, a, ja, ia, f, x, ijob, iprint, nrest, iter, tol) );
      }
      else {
      }
    }
    preprocessed = 0;
    retval(1) = 0;
    return retval ;
  }
  if (ijb==0 || ijb==1 || ijb==2) {
    if ( args(0).isempty ()  || ! args(0).issparse () ) {
      octave_stdout << " Error: the input matrix A must be a sparse nonempty"
                    << " array of double real or double complex\n";
    return retval;
    }
    else if ( args(0).is_real_matrix () ) {
      if (ijb==2 && notcpl==0) {
        octave_stdout << " Error: when IJOB==2 or IJOB==102, the input matrix A"
                      << " must have same type (real or complex) as on"
                      << " previous call with IJOB==1 or IJOB==101\n";
        return retval;
      }
      notcpl = 1;
      nr = args(0).rows ();
      nc = args(0).columns ();
      if (nr != nc) {
        octave_stdout << " Error: the input Matrix A must be square\n";
        return retval;
      }
      if (ijb==2 && nr!=n) {
        octave_stdout << " Error: when IJOB==2 or IJOB==102, the input matrix A"
                      << " must have same dimensions as on previous call with"
                      << " IJOB==1 or IJOB==101\n";
        return retval;
      }
      n = nr;
    }
    else if ( args(0).is_complex_matrix () ) {
      if (ijb==2 && notcpl==1) {
        octave_stdout << " Error: When IJOB==2 or IJOB==102, the input matrix A"
                      << " must have same type (real or complex) as on previous"
                      << " call with IJOB==1 or IJOB==101\n";
        return retval;
      }
      notcpl = 0;
      nr = args(0).rows ();
      nc = args(0).columns ();
      if (nr != nc) {
        octave_stdout << " Error: the input matrix A must be square\n";
        return retval;
      }
      if (ijb==2 && nr!=n) {
        octave_stdout << " Error: When IJOB==2 or IJOB==102, the input matrix A"
                      << " must have same dimensions as on previous call with"
                      << " IJOB==1 or IJOB==101\n";
        return retval;
      }
      n = nr;
    }
    else {
      octave_stdout << " Error: the input matrix A must be a sparse nonempty"
                    << " array of double real or double complex\n";
    return retval;
    }
    if (ijb == 0) preprocessed = 0;
  }
  if (ijb != 1) {
    if (args(1).isempty () ) {
      octave_stdout << " Error: the input right hand side vector B must be nonempty\n";
      return retval ;
    }
    if (args(1).issparse () ) {
      octave_stdout << " Error: the input right hand side B cannot be a sparse vector\n"
                    << "        Please convert it [b=full(b);] before calling AGMG\n";
      return retval ;
    }
    if ( args(1).is_real_matrix () ) {
      irhsreal = 1;
    }
    else if ( args(1).is_complex_matrix () ) {
      if (notcpl == 1) {
        octave_stdout << " Error: for real input matrix A, the input right hand"
                      << " side vector B must be real\n";
        return retval;
      }
      irhsreal = 0;
    }
    else {
      octave_stdout << " Error: the input right hand side B must be a vector of"
                    << " double real or (in case of complex input matrix A) double complex\n";
      return retval;
    }
    nrrhs = args(1).rows ();
    ncrhs = args(1).columns ();
    if (nrrhs != n || ncrhs != 1) {
      octave_stdout << " Error: the input right hand side B must be a column"
                    << " vector with the same number of rows as the input matrix A\n";
      return retval;
    }
    ix0=0;
    if (nargin > 6){
      if (! args(6).isempty () ) {
        if (ijb == 3) {
          octave_stdout << " Warning: the input argument X0 is ignored when IJOB==3\n";
        }
        else if (args(6).issparse () ) {
          octave_stdout << " Error: the input initial approximation X0 cannot"
                        << " be a sparse vector\n";
          return retval ;
        }
        else if ( args(6).is_real_matrix () ) {
          ix0 = 1;
        }
        else if ( args(6).is_complex_matrix () ) {
          if (notcpl == 1) {
            octave_stdout << " Error: for real input matrix A, the input"
                          << " initial approximation X0 must be real\n";
           return retval ;
          }
          ix0=-1;
        }
        else {
          octave_stdout << " Error: the input initial approximation X0 must be"
                        << " a vector of double real or (in case of complex"
                        << " input matrix A) double complex\n";
          return retval;
        }
	if (ix0 != 0) {
          nr = args(6).rows ();
          nc = args(6).columns ();
          if (nr != n || nc != 1) {
             octave_stdout << " Error: the input initial approximation X0 must"
                           << " be a column vector with the same number of rows"
                           << " as the input matrix A\n";
             return retval;
	  }
	  ijob=ijob+10;
        }
      }
    }
  }
  if (notcpl == 1) {
    if (ijb==1) {
      SparseMatrix ai = args(0).sparse_matrix_value ();
      jja = ai.ridx();
      iia = ai.cidx();
      a = ai.data();
      nz = iia[n];;
      OCTAVE_LOCAL_BUFFER (int, ia, n+1);
      OCTAVE_LOCAL_BUFFER (int, ja, nz);
      for (i=0 ; i<=n ; i++) {
         ia[i] = iia[i]+1;}
      for (i=0 ; i<nz ; i++) {
         ja[i] = jja[i]+1;}
      F77_XFCN (dagmg, FORTSUB,
               (n, a, ja, ia, f, x, ijob, iprint, nrest, iter, tol) );
      preprocessed = 1;
      retval(1) = 0;
      return retval;
    }
    nn=n;
    if (ijb==3 || ijb==202) {
      Matrix fi = args(1).matrix_value ();
      f = fi.fortran_vec ();
      if (ix0 == 0) {
        Matrix xm(nn,1);
        x = xm.fortran_vec ();
        F77_XFCN (dagmg, FORTSUB,
                 (n, a, ja, ia, f, x, ijob, iprint, nrest, iter, tol) );
	if (ijb==3) {
          retval(0) = xm;
          retval(1) = 0;
	  return retval;
	}
        else {
	  PROCOUT(f[i]);
	}
      }
      else {
	Matrix xm = args(6).matrix_value ();
        x = xm.fortran_vec ();
        F77_XFCN (dagmg, FORTSUB,
                 (n, a, ja, ia, f, x, ijob, iprint, nrest, iter, tol) );
	PROCOUT(f[i]);
      }
    }
    SparseMatrix ai = args(0).sparse_matrix_value ();
    jja = ai.ridx();
    iia = ai.cidx();
    a = ai.data();
    nz = iia[n];
    OCTAVE_LOCAL_BUFFER (int, ia, n+1);
    OCTAVE_LOCAL_BUFFER (int, ja, nz);
    for (i=0 ; i<=n ; i++) {
         ia[i] = iia[i]+1;}
    for (i=0 ; i<nz ; i++) {
	 ja[i] = jja[i]+1;}
    Matrix fi = args(1).matrix_value ();
    f = fi.fortran_vec ();
    if (ix0 == 0) {
      Matrix xm(nn,1);
      x = xm.fortran_vec ();
      F77_XFCN (dagmg, FORTSUB,
               (n, a, ja, ia, f, x, ijob, iprint, nrest, iter, tol) );
      PROCOUT(f[i]);
    }
    else {
      Matrix xm = args(6).matrix_value ();
      x = xm.fortran_vec ();
      F77_XFCN (dagmg, FORTSUB,
               (n, a, ja, ia, f, x, ijob, iprint, nrest, iter, tol) );
      PROCOUT(f[i]);
    }
  }
  else {
    if (ijb==1) {
      SparseComplexMatrix ai = args(0).sparse_complex_matrix_value ();
      jja = ai.ridx();
      iia = ai.cidx();
      ac = ai.data();
      nz = iia[n];;
      OCTAVE_LOCAL_BUFFER (int, ia, n+1);
      OCTAVE_LOCAL_BUFFER (int, ja, nz);
      for (i=0 ; i<=n ; i++) {
         ia[i] = iia[i]+1;}
      for (i=0 ; i<nz ; i++) {
         ja[i] = jja[i]+1;}
      F77_XFCN (zagmg, FORTSUB,
               (n, ac, ja, ia, fc, xc, ijob, iprint, nrest, iter, tol) );
      preprocessed = 1;
      retval(1) = 0;
      return retval;
    }
    nn=n;
    if (ijb==3 || ijb==202) {
      ComplexMatrix fi = args(1).complex_matrix_value ();
      fc = fi.fortran_vec ();
      if (ix0 == 0) {
        ComplexMatrix xm(nn,1);
        xc = xm.fortran_vec ();
        F77_XFCN (zagmg, FORTSUB,
                 (n, ac, ja, ia, fc, xc, ijob, iprint, nrest, iter, tol) );
	if (ijb==3) {
          retval(0) = xm;
          retval(1) = 0;
	  return retval;
	}
        else {
	  PROCOUT(fc[i].real());
	}
      }
      else {
	ComplexMatrix xm = args(6).complex_matrix_value ();
        xc = xm.fortran_vec ();
        F77_XFCN (zagmg, FORTSUB,
                 (n, ac, ja, ia, fc, xc, ijob, iprint, nrest, iter, tol) );
	PROCOUT(fc[i].real());
      }
    }
    SparseComplexMatrix ai = args(0).sparse_complex_matrix_value ();
    jja = ai.ridx();
    iia = ai.cidx();
    ac = ai.data();
    nz = iia[n];
    OCTAVE_LOCAL_BUFFER (int, ia, n+1);
    OCTAVE_LOCAL_BUFFER (int, ja, nz);
    for (i=0 ; i<=n ; i++) {
         ia[i] = iia[i]+1;}
    for (i=0 ; i<nz ; i++) {
	 ja[i] = jja[i]+1;}
    ComplexMatrix fi = args(1).complex_matrix_value ();
    fc = fi.fortran_vec ();
    if (ix0 == 0) {
      ComplexMatrix xm(nn,1);
      xc = xm.fortran_vec ();
      F77_XFCN (zagmg, FORTSUB,
               (n, ac, ja, ia, fc, xc, ijob, iprint, nrest, iter, tol) );
      PROCOUT(fc[i].real());
    }
    else {
      ComplexMatrix xm = args(6).complex_matrix_value ();
      xc = xm.fortran_vec ();
      F77_XFCN (zagmg, FORTSUB,
               (n, ac, ja, ia, fc, xc, ijob, iprint, nrest, iter, tol) );
      PROCOUT(fc[i].real());
    }
  }
}
