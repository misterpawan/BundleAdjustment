! COPYRIGHT (c) 2012-2018 Yvan Notay - ULB
!
! This file is part of AGMG software package
! Release 3.3.4-aca built on "Oct 22 2018" by Yvan Notay
!
! ALL USAGE OF AGMG IS SUBJECT TO LICENSE. PLEASE REFER TO THE FILE "LICENSE".
! IF YOU OBTAINED A SCOPY OF THIS SOFTWARE WITHOUT THIS FILE,
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
! This file provides sagmgpar and dependencies; sagmgpar is a
! parallel (MPI) implementation for real matrices in single precision
! of the method presented in [1], where the used algorithms are described
! in detail. From realease 3.x, the coarsening has been modified according
! to the results in [2,3], see the release notes in the README file
! for more details.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! SUBROUTINE sagmgpar (MAIN DRIVER): see bottom of this file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! PARAMETERS DEFINITON -  may be tuned by expert users
  MODULE sagmgpar_mem
    SAVE
! INTEGER
!  maxlev  is the maximal number of levels
!          (should be large enough - much larger than needed is armless).
!  real_len is the length of 1 REAL(kind(0.0e0)) in byte
!        (used only to display information on memory usage).
!  nsmooth  is the number of pre- and post-smoothing steps
!       (see smoothtype for details).
!  smoothtype indicates which smoother is used:
!    if smoothtype==1, the smoother is based on Gauss-Seidel sweeps;
!    if smoothtype==0, the smoother is based on SOR sweeps with automatic
!       estimation of the relaxation parameter (often reduces to Gauss-Seidel);
!    if smoothtype==-1, the smoother is based on SOR sweeps with relaxation
!       parameter specified in constant omega defined below;
!     Scheme used in all three cases: nsmooth sweeps for both pre- and
!     post-smoothing, with:
!        pre-smoothing: Forward sweep, then Backward sweep, then Forward, etc
!       post-smoothing: Backward sweep, then Forward sweep, then Backward, etc;
!    if smoothtype==2, the smoother is ILU(0), with
!      (nsmooth/2) pre- and (nsmooth/2+mod(nsmooth,1)) post-smoothing steps.
!  nstep   is the maximum number of coarsening steps;
!          nstep<0 means that coarsening is stopped only according to
!          the matrix size, see parameter maxcoarsesize.
!  nlvcyc  is the number of coarse levels (from bottom) on which V-cycle
!          formulation is enforced (Rmk: K-cycle always allowed at first
!          coarse level).
!  npass   is the maximal number of pairwise aggregation passes for each
!          coarsening step, according to the algorithms in [2,3].
!  maxcoarsesize: when the (global) size of the coarse grid matrix is less
!                 than or equal to maxcoarsesize*(N/#Proc.)^(1/3),
!                 it is factorized exactly and coarsening is stopped;
!         maxcoarsesizeslow: in case of slow coarsening,
!                 exact factorization can be performed when the size of
!                 the coarse grid matrix is less than or equal to
!                 maxcoarsesizeslow*(N/#Proc.)^(1/3).
!         (N is the global number of rows of the input matrix)
!         Remark: these parameters can be tuned in subroutine sagmgpar_init
!                 according to the number of concurrent tasks NPROC.
    INTEGER, PARAMETER :: maxlev=40, real_len=4
    INTEGER, PARAMETER :: nsmooth=1, smoothtype=0, nstep=-1, nlvcyc=0
    INTEGER, PARAMETER :: npass=3,maxcoarsesize=40,maxcoarsesizeslow=400
! REAL
!  omega is the relaxation parameter for SOR sweeps (smoothtype=-1)
!  resi is the threshold for the relative residual error in inner FCG & GCR
!       iterations, see Algorithm 3.2 in [1]
!  trspos is a threshold: if a row has a positive offdiagonal entry larger
!         than trspos times the diagonal entry, the corresponding node is
!         transferred unaggregated to the coarse grid
!  kaptg_ is the threshold used to accept or not a tentative aggregate
!         when applying the coarsening algorithms from [2,3];
!         kaptg_blocdia is used for control based on bloc diagonal smoother [2];
!         kaptg_dampJac is used for control based on Jacobi smoother [3].
!  checkdd is the threshold to keep outside aggregation nodes where
!         the matrix is strongly diagonally dominant (based on mean of row
!         and column);
!         In fact, AGMG uses the maximum of |checkdd| and of the default value
!            according to kaptg_ as indicated in [2,3]
!            (hence |checkdd| < 1 ensures that one uses this default value)
!         checkdd <0 : consider |checkdd|, but base the test on minus the
!               sum of offdiagonal elements, without taking the absolute value.
!  targetcoarsefac is the target coarsening factor (parameter tau in the main
!         coarsening algorithm in [2,3]): further pairwise aggregation passes
!         are omitted once the number of nonzero entries has been reduced by a
!         factor of at least targetcoarsefac.
!  fracnegrcsum: if, at some level, more than fracnegrcsum*nl nodes,
!         where nl is the total number of nodes at that level, have
!         negative mean row and column sum, then the aggregation algorithm
!         of [2,3] is modified, exchanging all diagonal entries for the mean
!         row and column sum (that is, the algorithm is applied to a
!         modified matrix with mean row and column sum enforced to be zero).
! K_ItGramS is the K factor in the rule to decide if a second iteration of
!         Gram Schmidt orthogonalization is needed (used within GCR)
    REAL(kind(0.0e0)), PARAMETER :: omega=0.8, resi=0.2, trspos=0.45
    REAL(kind(0.0e0)), PARAMETER :: kaptg_blocdia=4+2*npass, kaptg_dampJac=6+2*npass
    REAL(kind(0.0e0)), PARAMETER :: checkdd=0.5
    REAL(kind(0.0e0)), PARAMETER :: targetcoarsefac=2.0**npass
    REAL(kind(0.0e0)), PARAMETER :: fracnegrcsum=0.25
    REAL(kind(0.0e0)), PARAMETER :: K_ItGramS=10.0
!!!!!!!!!!!!!!!!!!!! END of PARAMETERS DEFINITION -----------------
!!!!!!!!!!!!!!!!!!! Internal variables declaration
!
    TYPE InnerData
       REAL(kind(0.0e0)), DIMENSION(:), POINTER :: a
       INTEGER, DIMENSION(:), POINTER :: ja
       INTEGER, DIMENSION(:), POINTER :: ia
       INTEGER, DIMENSION(:), POINTER :: il
       INTEGER, DIMENSION(:), POINTER :: iu
       REAL(kind(0.0e0)), DIMENSION(:), POINTER :: p
       INTEGER, DIMENSION(:), POINTER :: idiag
       INTEGER, DIMENSION(:), POINTER :: ind
       INTEGER, DIMENSION(:), POINTER :: iext
       INTEGER, DIMENSION(:), POINTER :: ilstout
       INTEGER, DIMENSION(:), POINTER :: lstout
       INTEGER, DIMENSION(:), POINTER :: ilstin
       INTEGER, DIMENSION(:), POINTER :: lstin
       INTEGER, DIMENSION(:), POINTER :: iblockl
    END TYPE InnerData
!
    TYPE(InnerData) :: dt(maxlev)
    REAL(kind(0.0e0)), ALLOCATABLE :: scald(:)
    INTEGER :: nn(maxlev),kstat(2,maxlev)=0,innermax(maxlev)
    INTEGER :: nlev,nwrkcum,iout,nrst,nbblock
    INTEGER :: maxcoarset,maxcoarseslowt,smoothtp
    REAL(kind(0.0e0)) :: wcplex(4),fracnz(maxlev),omeg=1.0e0,flop=0.0e0
    LOGICAL :: spd,wfo,wff,woo,allzero,trans,transint,zerors,gcrin
    LOGICAL :: coasing,coarsit,spdcoarse
    REAL(kind(0.0e0)), PARAMETER :: cplxmax=3.0, xsi=0.6
    REAL(kind(0.0e0)), PARAMETER :: repsmach=SQRT(EPSILON(1.0e0)),epsmach=EPSILON(1.0e0)
    REAL(kind(0.0e0)), PARAMETER :: UNP=1.0e0+10*epsmach
    INTEGER :: nlc(2),nlcp(2),nlc1(2),icum,npassr
    INTEGER :: imult
    REAL(kind(0.0e0)) :: ngl(2),nglp(2),nlctot(2),ngl1(2),ngltot(2),RELRESL1
    INTEGER, DIMENSION(:), POINTER :: iextL1,ilstoutL1,lstoutL1,ilstinL1
    INTEGER, DIMENSION(:), POINTER :: lstinL1,ipermL1
    LOGICAL :: cvtetodo,converged,timexitl,timexitg
    INTEGER(selected_int_kind(8)),ALLOCATABLE  :: ireqi(:),ireqo(:)
    INTEGER, ALLOCATABLE :: ineigh(:)
    REAL(kind(0.0e0)), ALLOCATABLE :: buffi(:),buffo(:)
    INTEGER :: IRANK,ICOMM,NPROC,nneigh,nneighi,nneigho=0,neffn
    INCLUDE 'mpif.h'
    REAL(kind(0.0e0)) :: cputm=0.0e0,eltm=0.0e0,eltmp=0.0e0
    LOGICAL :: preprocessed=.FALSE., firstcall=.TRUE., solve=.FALSE.
    REAL(kind(0.0e0)), PARAMETER ::    &
         checkddJ=MAX(ABS(checkdd),kaptg_dampJac/(kaptg_dampJac-2))
    REAL(kind(0.0e0)), PARAMETER ::    &
         checkddB=MAX(ABS(checkdd),(kaptg_blocdia+1)/(kaptg_blocdia-1))
  END MODULE sagmgpar_mem
!!!!!!!!!!!!!!!!!!! END of Internal variables declaration
!!!!!!!!!!!!!!!!!!! SOME INITIALIZATION for the parallel version
  SUBROUTINE sagmgpar_init(MPI_COMM_AGMG)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: MPI_COMM_AGMG
    INTEGER :: ierr
    ICOMM=MPI_COMM_AGMG
    CALL MPI_COMM_RANK(MPI_COMM_AGMG,IRANK,ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_AGMG,NPROC,ierr)
  ! tuning of maxcoarse & maxcoarseslow according NPROC
    maxcoarset=maxcoarsesize
    maxcoarseslowt=maxcoarsesizeslow
    RETURN
  END SUBROUTINE sagmgpar_init
!!!!!!!!!!!!!!!!!!! END of SOME INITIALIZATION for parallel version
!!!!!!!!!!!!!!!!!!! TIMING
  SUBROUTINE sagmgpar_mestime(id,cputm,eltm)
    IMPLICIT NONE
    INTEGER, SAVE :: cpt_init(10)=-1,cpt_fin,cpt_max,freq,cpt
    REAL, SAVE :: t1(10), t2
    REAL(kind(0.0e0)) :: cputm,eltm
    INTEGER :: id
    IF (id.GT.0) THEN
       !Next line may be uncommented if FORTRAN 95 function
       !CPU_TIME is implemented
       !   CALL CPU_TIME(t2)
       CALL SYSTEM_CLOCK(cpt_fin,freq,cpt_max)
       cpt = cpt_fin - cpt_init(id)
       IF (cpt_fin.LT.cpt_init(id)) cpt = cpt + cpt_max
       eltm = real(cpt) / freq
       cputm = real(t2 - t1(id))
    ELSE
       CALL SYSTEM_CLOCK(cpt_init(-id),freq,cpt_max)
       !Next line may be uncommented if FORTRAN 95 function
       !CPU_TIME is implemented
       !   CALL CPU_TIME(t1(-id))
    END IF
    RETURN
  END SUBROUTINE sagmgpar_mestime
!!!!!!!!!!!!!!!!!!! END of TIMING
  MODULE sagmgpar_PARDISSO
    TYPE PARDISO_HANDLE
       INTEGER(KIND=8) DUMMY
    END TYPE PARDISO_HANDLE
  END MODULE sagmgpar_PARDISSO
  MODULE sagmgpar_PARDISO
    USE sagmgpar_PARDISSO
    INTERFACE
       SUBROUTINE PARDISO( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA  &
            , PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
         USE sagmgpar_PARDISSO
         TYPE(PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
         INTEGER,          INTENT(IN)    :: MAXFCT
         INTEGER,          INTENT(IN)    :: MNUM
         INTEGER,          INTENT(IN)    :: MTYPE
         INTEGER,          INTENT(IN)    :: PHASE
         INTEGER,          INTENT(IN)    :: N
         INTEGER,          INTENT(IN)    :: IA(*)
         INTEGER,          INTENT(IN)    :: JA(*)
         INTEGER,          INTENT(IN)    :: PERM(*)
         INTEGER,          INTENT(IN)    :: NRHS
         INTEGER,          INTENT(INOUT) :: IPARM(*)
         INTEGER,          INTENT(IN)    :: MSGLVL
         INTEGER,          INTENT(OUT)   :: ERROR
         REAL(kind(0.0e0)), INTENT(IN)    :: A(*)
         REAL(kind(0.0e0)), INTENT(INOUT) :: B(*)
         REAL(kind(0.0e0)), INTENT(OUT)   :: X(*)
       END SUBROUTINE PARDISO
    END INTERFACE
  END MODULE sagmgpar_PARDISO
MODULE sagmgpar_ALLROUTINES
PRIVATE
PUBLIC :: sagmgpar_setupL1, sagmgpar_relmem,  sagmgpar_smoothsetup
PUBLIC :: sagmgpar_applyprec, sagmgpar_GCR, sagmgpar_FlexCG
PUBLIC :: sagmgpar_setsendrecL1
PUBLIC :: sagmgpar_partroword
PUBLIC :: sagmgpar_restorelstrank
CONTAINS
    SUBROUTINE PARDISO( PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA    &
         , PERM, NRHS, IPARM, MSGLVL, B, X, ERROR )
      USE sagmgpar_PARDISSO
      TYPE(PARDISO_HANDLE), INTENT(INOUT) :: PT(*)
      INTEGER,          INTENT(IN)    :: MAXFCT
      INTEGER,          INTENT(IN)    :: MNUM
      INTEGER,          INTENT(IN)    :: MTYPE
      INTEGER,          INTENT(IN)    :: PHASE
      INTEGER,          INTENT(IN)    :: N
      INTEGER,          INTENT(IN)    :: IA(*)
      INTEGER,          INTENT(IN)    :: JA(*)
      INTEGER,          INTENT(IN)    :: PERM(*)
      INTEGER,          INTENT(IN)    :: NRHS
      INTEGER,          INTENT(INOUT) :: IPARM(*)
      INTEGER,          INTENT(IN)    :: MSGLVL
      INTEGER,          INTENT(OUT)   :: ERROR
      REAL(kind(0.0e0)), INTENT(IN)    :: A(*)
      REAL(kind(0.0e0)), INTENT(INOUT) :: B(*)
      REAL(kind(0.0e0)), INTENT(OUT)   :: X(*)
      ERROR=1
      X(1:N)=0.0e0
      STOP
    END SUBROUTINE PARDISO
  SUBROUTINE sagmgpar_relmem(ijb)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: l, ier, ijb
    REAL(kind(0.0e0)) :: fdum(1)
    DO l=1,nlev-1
       IF(nn(l).GT.0) THEN
            IF (nn(l+1).GT.0)     &
                      DEALLOCATE(dt(l)%ind)
            DEALLOCATE(dt(l)%a,dt(l)%ja,dt(l)%il,dt(l)%iu)
            IF (smoothtp.NE.1) DEALLOCATE(dt(l)%p)
       END IF
    END DO
    IF (nn(nlev).GT.0 .AND. (nlev.GT.1 .OR. ijb.LT.0) )        &
              DEALLOCATE(dt(nlev)%a,dt(nlev)%ja,dt(nlev)%ia)
    IF (ijb.GE.0) THEN
       IF (nlev.GT.1) THEN
          DEALLOCATE(iextL1)
       ELSE
          NULLIFY(iextL1)
       END IF
    END IF
    DO l=1,nlev
       IF (nn(l).GT.0) DEALLOCATE(dt(l)%lstout,dt(l)%lstin)
       IF (nn(l).GT.0)                  &
          DEALLOCATE(dt(l)%iext,dt(l)%ilstin,dt(l)%ilstout)
    END DO
    IF (nneigho.GT.0) CALL MPI_WAITALL(nneigho,ireqo,MPI_STATUSES_IGNORE,ier)
    DEALLOCATE (ireqi,ireqo,ineigh)
    DEALLOCATE (buffi,buffo)
    nneigho=0
    IF (.NOT.allzero)                              &
            CALL sagmgpar_MUMPSpar(nn(nlev),fdum,-2)
    RETURN
  END SUBROUTINE sagmgpar_relmem
  SUBROUTINE sagmgpar_applyprec( N,f,X,ijb,a,ja,ia)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: N
    REAL(kind(0.0e0)), OPTIONAL :: f(N), X(N)
    REAL(kind(0.0e0)), OPTIONAL :: a(*)
    INTEGER, OPTIONAL :: ja(*), ia(N+1), ijb
    REAL(kind(0.0e0)), ALLOCATABLE, SAVE :: S(:)
    REAL(kind(0.0e0)) :: tau=1.35
    INTEGER :: iter, nvec, i
       iter=ijb-2
       nvec=min(2,iter)
       ALLOCATE (S(nwrkcum+nvec*N))
    IF (iter.GT.1) THEN
       CALL sagmgpar_prec_matv(1,f,X,S(N+1),S(2*N+1),.TRUE.)
       X(1:N)=tau*X(1:N)
       DO i=2,iter-1
          CALL sagmgpar_prec_matv(1,f,S,S(N+1),S(2*N+1),.TRUE.)
          X(1:N)=X(1:N)+tau*S(1:N)
          f(1:N)=f(1:N)-tau*S(N+1:2*N)
       END DO
       CALL sagmgpar_prec_matv(1,f,S,S(N+1),S(2*N+1),.FALSE.)
       X(1:N)=X(1:N)+tau*S(1:N)
    ELSE
       CALL sagmgpar_prec_matv(1,f,X,S,S(N+1),.FALSE.)
    END IF
    kstat(2,1)=kstat(2,1)+iter
    DEALLOCATE(S)
    RETURN
  END SUBROUTINE sagmgpar_applyprec
  SUBROUTINE sagmgpar_partroword(n, a, ja, ia, idiag, iext, nzext)
    IMPLICIT NONE
    INTEGER :: n, ia(n+1), idiag(n)
    INTEGER, OPTIONAL :: iext(*), nzext
    INTEGER :: ja(*)
    INTEGER, ALLOCATABLE, TARGET, DIMENSION(:)  :: iw, jcol
    REAL(kind(0.0e0)) :: a(*), p
    REAL(kind(0.0e0)), ALLOCATABLE, TARGET :: w(:)
    INTEGER :: i, j, k, ipos, nzu, nze, ila, mal, mil, ifi
    INTEGER, POINTER, DIMENSION(:) :: iwe
    REAL(kind(0.0e0)), POINTER, DIMENSION(:) :: we
    LOGICAL :: wnoloc
    INTEGER, PARAMETER :: ifim1=0
    wnoloc = .TRUE.
    ALLOCATE (w(2*n),iw(2*n),jcol(n))
    iwe => iw(n+1:2*n)
    we  =>  w(n+1:2*n)
    nzext=0
    jcol=0
    ipos=ia(1)
    DO i=1,n
       nzu=0
       nze= 0
       p=0.0e0
       ifi=ia(i)
       ia(i)=ipos
       DO k=ifi,ia(i+1)-1
          j=ja(k)
          IF (j.GT.n) THEN
             nze=nze+1
             IF (nze.GT.n .AND. wnoloc) THEN
                ALLOCATE(iwe(ia(n+1)-ia(1)),we(ia(n+1)-ia(1)))
                iwe(1:n)=iw(n+1:2*n)
                we(1:n)=w(n+1:2*n)
                wnoloc=.FALSE.
             END IF
             we(nze)=a(k)
             iwe(nze)=j
             CYCLE
          END IF
          j=j-ifim1
          IF (j.GT.i) THEN
             IF (jcol(j).EQ.0) THEN
               nzu=nzu+1
               w(nzu)=a(k)
               iw(nzu)=j
               jcol(j)=nzu
             ELSE
               w(jcol(j))=w(jcol(j))+a(k)
             END IF
          ELSE IF (j.LT.i) THEN
            IF (jcol(j).EQ.0) THEN
               a(ipos)=a(k)
               ja(ipos)=j
               jcol(j)=ipos
               ipos=ipos+1
            ELSE
               a(jcol(j))=a(jcol(j))+a(k)
            END IF
          ELSE
             p=p+a(k)
          END IF
       END DO
       DO k=ifi,ipos-1
          jcol(ja(k))=0
       END DO
       idiag(i)=ipos
       a(ipos)=p
       ja(ipos)=i
       ipos=ipos+1
       DO k=1,nzu
          a(ipos)=w(k)
          ja(ipos)=iw(k)
          jcol(iw(k))=0
          ipos=ipos+1
       END DO
       nzext=nzext+nze
       iext(i)=ipos
       DO k=1,nze
          a(ipos)=we(k)
          ja(ipos)=iwe(k)
          ipos=ipos+1
       END DO
    END DO
    ia(n+1)=ipos
    IF(.NOT.wnoloc) DEALLOCATE(iwe,we)
    DEALLOCATE (w,iw,jcol)
    RETURN
  END SUBROUTINE sagmgpar_partroword
  SUBROUTINE sagmgpar_csrdlu(n,a,ja,ia,idiag,ao,jao,iup,iext,iextop)
    IMPLICIT NONE
    INTEGER :: n, ja(*), jao(n+1:*)
    INTEGER, TARGET :: ia(n+1), idiag(n+1)
    INTEGER, OPTIONAL, TARGET :: iup(n+1), iext(n+1),iextop(n+1)
    REAL(kind(0.0e0)) :: a(*),ao(*)
    INTEGER :: i,j,ili,iui,iei,ipos,ili0,iui0,iei0,nl,nu,k,next
    INTEGER, POINTER :: il(:), iu(:), iexto(:)
    nl=0
    nu=0
    il => idiag(1:n+1)
    IF (PRESENT(iextop)) THEN
       iexto => iextop(1:n+1)
       iu => iup(1:n+1)
       DO i=1,n
          ili=idiag(i)-ia(i)
          iui=iext(i)-idiag(i)-1
          iei=ia(i+1)-iext(i)
          iextop(i)=iei
          iup(i)=iui
          idiag(i)=ili
          nl=nl+ili
          nu=nu+iui
       END DO
    ELSE
       iexto => iext(1:n+1)
       iu => ia(1:n+1)
       DO i=1,n
          ili=idiag(i)-ia(i)
          iui=iext(i)-idiag(i)-1
          iei=ia(i+1)-iext(i)
          iext(i)=iei
          ia(i)=iui
          idiag(i)=ili
          nl=nl+ili
          nu=nu+iui
       END DO
    END IF
    ipos=0
    ili=n+1
    iui=ili+nl
    iei=iui+nu
       DO i=1,n
          ili0=ili
          iui0=iui
          DO j=1,il(i)
             ipos=ipos+1
             ao(ili)=a(ipos)
             jao(ili)=ja(ipos)
             ili=ili+1
          END DO
          ipos=ipos+1
          ao(i)=a(ipos)
          DO j=1,iu(i)
             ipos=ipos+1
             ao(iui)=a(ipos)
             jao(iui)=ja(ipos)
             iui=iui+1
          END DO
          iei0=iei
          DO j=1,iexto(i)
             ipos=ipos+1
             ao(iei)=a(ipos)
             jao(iei)=ja(ipos)
             iei=iei+1
          END DO
          iexto(i)=iei0
          iu(i)=iui0
          il(i)=ili0
       END DO
       iexto(n+1)=iei
       iu(n+1)=iui
       il(n+1)=ili
    RETURN
  END SUBROUTINE sagmgpar_csrdlu
  SUBROUTINE sagmgpar_csrdluT(n,a,ja,ia,idiag,ao,jao,il,iu)
    IMPLICIT NONE
    INTEGER :: n, ja(*), ia(n+1), idiag(n), jao(n+1:*), il(n+1), iu(n+1)
    REAL(kind(0.0e0)) :: a(*),ao(*)
    INTEGER :: i,j,ili,iui,iei,ipos,ili0,iui0,iei0,nl,nu,k,next
    nl=0
    nu=0
       il(2:n+1)=0
       iu(2:n+1)=0
       DO i=1,n
          iui=idiag(i)-ia(i)
          ili=ia(i+1)-idiag(i)-1
          nl=nl+ili
          nu=nu+iui
          DO k=ia(i),idiag(i)-1
             iu(ja(k)+1)=iu(ja(k)+1)+1
          END DO
          DO k=idiag(i)+1,ia(i+1)-1
             il(ja(k)+1)=il(ja(k)+1)+1
          END DO
       END DO
    ipos=0
    ili=n+1
    iui=ili+nl
    iei=iui+nu
       il(1)=ili
       iu(1)=iui
       DO i=1,n
          il(i+1) = il(i) + il(i+1)
          iu(i+1) = iu(i) + iu(i+1)
       END DO
       DO i=1,n
          DO k=ia(i),idiag(i)-1
             j = ja(k)
             next = iu(j)
             ao(next) = a(k)
             jao(next) = i
             iu(j) = next+1
          END DO
          ao(i)=a(idiag(i))
          DO k=idiag(i)+1,ia(i+1)-1
             j = ja(k)
             next = il(j)
             ao(next) = a(k)
             jao(next) = i
             il(j) = next+1
          END DO
       END DO
       DO i=n,1,-1
          il(i+1) = il(i)
          iu(i+1) = iu(i)
       END DO
       il(1)=ili
       iu(1)=iui
    RETURN
  END SUBROUTINE sagmgpar_csrdluT
  SUBROUTINE sagmgpar_csrmv(n, a, ja, ifja, ia, x, ifx, y, iad, flop, lstin)
    IMPLICIT NONE
    INTEGER :: n, ifja, ja(ifja:*), ia(n+1), iad, ifx, i, k, j
    INTEGER, OPTIONAL :: lstin(0:*)
    REAL(kind(0.0e0)) :: a(*), x(ifx:*), y(n), t
    REAL(kind(0.0e0)) :: flop
    IF (n.LE.0) RETURN
    IF (.NOT.PRESENT(lstin)) THEN
       IF (iad.LT.-1) THEN
          DO i=1,n
             t=y(i)
             DO k=ia(i),ia(i+1)-1
                t = t - a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
       ELSE IF (iad.GT.1) THEN
          DO i=1,n
             t=y(i)
             DO k=ia(i),ia(i+1)-1
                t = t + a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
       ELSE IF (iad.EQ.-1) THEN
          DO i=1,n
             t=0.0e0
             DO k=ia(i),ia(i+1)-1
                t = t - a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
       ELSE
          DO i=1,n
             t=0.0e0
             DO k=ia(i),ia(i+1)-1
                t = t + a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
       END IF
   ELSE
       IF (iad.LT.-1) THEN
          DO j=1,lstin(0)
             i=lstin(j)
             t=y(i)
             DO k=ia(i),ia(i+1)-1
                t = t - a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
       ELSE IF (iad.GT.1) THEN
          DO j=1,lstin(0)
             i=lstin(j)
             t=y(i)
             DO k=ia(i),ia(i+1)-1
                t = t + a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
       ELSE IF (iad.EQ.-1) THEN
          DO j=1,lstin(0)
             i=lstin(j)
             t=0.0e0
             DO k=ia(i),ia(i+1)-1
                t = t - a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
       ELSE
          DO j=1,lstin(0)
             i=lstin(j)
             t=0.0e0
             DO k=ia(i),ia(i+1)-1
                t = t + a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
       END IF
    END IF
    flop=flop+real(2*(ia(n+1)-ia(1)))
    RETURN
  END SUBROUTINE sagmgpar_csrmv
  SUBROUTINE sagmgpar_mcsrmv(n, a, ja, ifja, ia, x, ifx, y, iad, flop)
    IMPLICIT NONE
    INTEGER :: n, ifja, ja(ifja:*), ia(n+1), iad, ifx, i, k, j
    REAL(kind(0.0e0)) :: a(*), x(ifx:*), y(n), t
    REAL(kind(0.0e0)) :: flop
    IF (n.LE.0) RETURN
       IF (iad.LT.-1) THEN
          DO i=1,n
             t=y(i)-a(i)*x(i)
             DO k=ia(i),ia(i+1)-1
                t = t - a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
          flop=flop+real(2*(ia(n+1)-ia(1)+n))
       ELSE IF (iad.GT.1) THEN
          DO i=1,n
             t=y(i)+a(i)*x(i)
             DO k=ia(i),ia(i+1)-1
                t = t + a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
          flop=flop+real(2*(ia(n+1)-ia(1)+n))
       ELSE IF (iad.EQ.-1) THEN
          DO i=1,n
             t=-a(i)*x(i)
             DO k=ia(i),ia(i+1)-1
                t = t - a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
          flop=flop+real(2*(ia(n+1)-ia(1))+n)
       ELSE
          DO i=1,n
             t=a(i)*x(i)
             DO k=ia(i),ia(i+1)-1
                t = t + a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
          flop=flop+real(2*(ia(n+1)-ia(1))+n)
       END IF
    RETURN
  END SUBROUTINE sagmgpar_mcsrmv
  SUBROUTINE sagmgpar_csrlsolve(n, a, ja, ifja, ia, p, x, y, iunit, flop)
    IMPLICIT NONE
    INTEGER :: n, ifja, ja(ifja:*), ia(n+1), iunit, i, k
    REAL(kind(0.0e0)) :: a(*), p(n), x(n), y(n), t
    REAL(kind(0.0e0)) :: flop
    IF (n.LE.0) RETURN
    IF (iunit.LT.0) THEN
       x(1) = y(1)
       DO i=2,n
          t = y(i)
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))*p(ja(k))
          END DO
          x(i) = t
       END DO
       flop=flop+real(n+3*(ia(n+1)-ia(1)))
    ELSE IF (iunit.GT.0) THEN
       x(1) = p(1)*y(1)
       DO i=2,n
          t = y(i)
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = p(i)*t
       END DO
       flop=flop+real(n+2*(ia(n+1)-ia(1)))
    ELSE
       x(1) = y(1)
       DO i=2,n
          t = y(i)
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = t
       END DO
       flop=flop+real(2*(ia(n+1)-ia(1)))
    END IF
    RETURN
  END SUBROUTINE sagmgpar_csrlsolve
  SUBROUTINE sagmgpar_csrusolve(n, a, ja, ifja, ia, p, x, y, iunit, flop)
    IMPLICIT NONE
    INTEGER :: n, ifja, ja(ifja:*), ia(n+1), iunit, i, k
    REAL(kind(0.0e0)) :: a(*), p(n), x(n), y(n), t
    REAL(kind(0.0e0)) :: flop
    IF (n.LE.0) RETURN
    IF (iunit.LT.0) THEN
       x(n) = y(n)
       DO i=n-1,1,-1
          t = 0.0e0
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = y(i) + p(i)*t
       END DO
       flop=flop+real(n+2*(ia(n+1)-ia(1)))
    ELSE IF (iunit.GT.0) THEN
       x(n) = p(n)*y(n)
       DO i=n-1,1,-1
          t = y(i)
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = p(i)*t
       END DO
       flop=flop+real(n+2*(ia(n+1)-ia(1)))
    ELSE
       x(n) = y(n)
       DO i=n-1,1,-1
          t = y(i)
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = t
       END DO
       flop=flop+real(2*(ia(n+1)-ia(1)))
    END IF
    RETURN
  END SUBROUTINE sagmgpar_csrusolve
  SUBROUTINE sagmgpar_sendrec(vect,ivect,ishift,l)
    USE sagmgpar_mem
    IMPLICIT NONE
    REAL(kind(0.0e0)) :: vect(*)
    INTEGER :: ivect(*),ishift,l
    INTEGER :: i, ii, k, kdeb, kfin, ksh, ier, idest
    INTEGER, DIMENSION(:), POINTER :: lstout,ilstout,ilstin
    IF (l.LE.0) THEN
       lstout => lstoutL1
       ilstout => ilstoutL1
       ilstin => ilstinL1
    ELSE
       lstout => dt(l)%lstout
       ilstout => dt(l)%ilstout
       ilstin => dt(l)%ilstin
    END IF
   IF (nneigho.GT.0) CALL MPI_WAITALL(nneigho,ireqo,MPI_STATUSES_IGNORE,ier)
    ii=0
    DO i=1,nneigh
       kdeb=ilstout(i)
       kfin=ilstout(i+1)-1
       IF (kfin.GE.kdeb) THEN
          ii=ii+1
          DO k=kdeb,kfin
             IF (ishift.GE.0) THEN
                buffo(k)=real(ivect(lstout(k))+ishift)
             ELSE
                buffo(k)=vect(lstout(k))
             END IF
          END DO
          idest=ineigh(i)
          CALL MPI_ISEND(buffo(kdeb),kfin-kdeb+1,MPI_REAL,&
               idest,985,icomm,ireqo(ii),ier)
       END IF
    END DO
    nneigho=ii
    ii=0
    DO i=1,nneigh
       kdeb=ilstin(i)
       kfin=ilstin(i+1)-1
       IF (kfin.GE.kdeb) THEN
          ii=ii+1
          idest=ineigh(i)
          CALL MPI_IRECV(buffi(kdeb),kfin-kdeb+1,MPI_REAL,&
               idest,985,ICOMM,ireqi(ii),ier)
       END IF
    END DO
    nneighi=ii
    RETURN
  END SUBROUTINE sagmgpar_sendrec
  SUBROUTINE sagmgpar_ABORT
    USE sagmgpar_mem
    INTEGER :: errc,ierr
    CALL MPI_ABORT(ICOMM,errc,ierr)
    STOP
  END SUBROUTINE sagmgpar_ABORT
  SUBROUTINE sagmgpar_FlexCG(N,f,X,ITER,RESID,a,ja,ia,ijb)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER       :: N, ITER, ijb
    REAL(kind(0.0e0)) :: f(N), X(N)
    INTEGER       :: ja(*), ia(N+1)
    REAL(kind(0.0e0)) :: a(*)
    INTEGER       :: MAXIT, ierr, kk, k, idum(1), i
    REAL(kind(0.0e0)) :: TOL, BNORM, RESID, dum0, RESID0, FREL, RHO
    REAL(kind(0.0e0)) :: ALPHA, BET0, BET1, BET, dum(4), td
    REAL(kind(0.0e0)), ALLOCATABLE :: SD(:),S(:),fsc(:)
    LOGICAL :: init,intmv,intmv2
    REAL(kind(0.0e0)), external :: SDOT
    REAL(kind(0.0e0)), external :: SNRM2
    INTEGER , parameter :: IONE=1
    intmv=.FALSE.
    IF ( (MOD(ijb,10).EQ.0 .OR. (ijb.GT.201)) .AND. nlev.GT.1) intmv=.TRUE.
    intmv2=intmv .OR. ijb.GT.201
       ALLOCATE (S(2*N+1:nwrkcum+5*N),SD(N))
    flop=0.0e0
    kstat=0
    kstat(1,MAX(2,nlev))=1
    init=.FALSE.
    IF (MOD(ijb,100).GE.10) init=.TRUE.
    IF (wfo) THEN
       WRITE(iout,940) IRANK
    END IF
    IF (wff) THEN
       IF (  trans) THEN
          WRITE(iout,941)
       END IF
       IF (nlev.GT.1) THEN
          IF (smoothtp.EQ.-1) THEN
             WRITE(iout,943) nsmooth
             WRITE(iout,944) omeg
          ELSE IF (smoothtp.EQ.1) THEN
             WRITE(iout,945) nsmooth
             WRITE(iout,947)
          ELSE
             WRITE(iout,946) nsmooth+mod(nsmooth,2),nsmooth-mod(nsmooth,2)
             WRITE(iout,947)
          END IF
       END IF
    END IF
    TOL = RESID
    MAXIT = ITER
    RESID = 1.0e0
    ITER = 0
    dum(3) = SNRM2(N, f, IONE)**2
    flop=flop+real(2*N)
    IF (init) THEN
       IF (intmv) THEN
          CALL sagmgpar_sendrec(x,idum,-1,1)
          IF (smoothtp.EQ.1) THEN
             f(1:n)=f(1:n)-x(1:n)/dt(1)%a(1:n)
             CALL sagmgpar_csrmv(n,dt(1)%a,dt(1)%ja,n+1,dt(1)%il,x,1,f,-2,flop)
          ELSE
             CALL sagmgpar_mcsrmv(n,dt(1)%a,dt(1)%ja,n+1,dt(1)%il,x,1,f,-2,flop)
             f(1:n)=f(1:n)-x(1:n)/(dt(1)%p(1:n))
          END IF
          CALL sagmgpar_csrmv(n,dt(1)%a,dt(1)%ja,n+1,dt(1)%iu,x,1,f,-2,flop)
          IF (nneighi.GT.0) CALL MPI_WAITALL(nneighi,ireqi,MPI_STATUSES_IGNORE,ierr)
          CALL  sagmgpar_csrmv(n,dt(1)%a,dt(1)%ja,n+1,dt(1)%iext,buffi,n+1,f &
               ,-2,flop,dt(1)%lstin)
       ELSE
          IF (intmv2) THEN
             CALL sagmgpar_matv(N,x,SD,dt(1)%a,dt(1)%ja,dt(1)%ia,transint &
             ,dt(1)%iext,dt(1)%lstin,1)
          ELSE
             CALL sagmgpar_matv(N,x,SD,a,ja,ia,trans &
             ,iextL1,lstinL1,0)
          END IF
          f(1:n)=f(1:n)-SD(1:n)
       END IF
       dum(2) = SNRM2(N, f, IONE)**2
       flop=flop+real(2*N)
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,dum(2:3),2,MPI_REAL,  &
            MPI_SUM,ICOMM,ierr)
       BNORM=SQRT(dum(3))
       RESID=SQRT(dum(2))
       RESID0=RESID
       IF (BNORM.EQ.0.0e0) THEN
          IF(wff) THEN
             WRITE(iout,998)
             WRITE(iout,999)
          END IF
          X(1:N)=0.0e0
          DEALLOCATE(S,SD)
          IF(ALLOCATED(fsc)) DEALLOCATE(fsc)
          RETURN
       END IF
       RESID=RESID/BNORM
       IF(wff.AND. (MAXIT.LE.0 .OR. RESID.LE.TOL)) THEN
          WRITE(iout,900) 0, resid*bnorm, resid
       END IF
    END IF
    DO WHILE ( ITER.LT.MAXIT .AND. RESID.GT.TOL )
       ITER = ITER + 1
       IF (ijb.GT.201) THEN
          CALL sagmgpar_prec_matv(1,f,S(1+3*N),S(1+4*N),S(1+5*N),intmv)
       ELSE
          CALL sagmgpar_prec_matv(1,f,S(1+3*N),S(1+4*N),S(1+5*N),intmv,a,ja,ia)
       END IF
       kk=3
       IF ( ITER.GT.1 ) THEN
          kk=2
          dum(1) = - SDOT(N, SD(1:n), IONE, S(1+3*N:4*N), IONE)
          CALL MPI_ALLREDUCE(MPI_IN_PLACE,dum(1),1,MPI_REAL,  &
               MPI_SUM,ICOMM,ierr)
          BET0=dum(1)
          BET1=BET0/RHO
          CALL SAXPY(N, BET1, S(1+2*N:3*N), IONE, S(1+3*N:4*N), IONE)
          flop=flop+real(4*N)
          IF (intmv) THEN
             CALL SAXPY(N, BET1, SD(1:n), IONE, S(1+4*N:5*N), IONE)
             flop=flop+real(2*N)
          END IF
       END IF
       CALL SCOPY(N, S(1+3*N:4*N), IONE, S(1+2*N:3*N), IONE)
       IF (intmv)  THEN
          CALL SCOPY(N, S(1+4*N:5*N), IONE, SD(1:n), IONE)
       ELSE
          IF (intmv2) THEN
             CALL sagmgpar_matv(N,S(1+2*N),SD,dt(1)%a,dt(1)%ja,dt(1)%ia,trans &
             ,dt(1)%iext,dt(1)%lstin,1)
          ELSE
             CALL sagmgpar_matv(N,S(1+2*N),SD,a,ja,ia,trans &
             ,iextL1,lstinL1,0)
          END IF
       END IF
       dum(1) =  SDOT(N, S(1+2*N:3*N), IONE, SD(1:n), IONE)
       dum(2) =  SDOT(N,S(1+2*N:3*N),IONE,f,IONE)
       flop=flop+real(4*N)
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,dum(1:kk),kk,MPI_REAL,  &
               MPI_SUM,ICOMM,ierr)
       IF (ITER.EQ.1) THEN
          IF (.NOT.init) THEN
             BNORM=SQRT(dum(3))
             RESID0=BNORM
             IF (BNORM.EQ.0.0e0) THEN
                IF(wff) THEN
                   WRITE(iout,998)
                   WRITE(iout,999)
                END IF
                X(1:N)=0.0e0
                ITER=0
                DEALLOCATE(S,SD)
                IF(ALLOCATED(fsc)) DEALLOCATE(fsc)
                RETURN
             END IF
             IF(wff) THEN
                WRITE(iout,900) 0, resid*bnorm, resid
             END IF
             IF ( RESID.LE.TOL ) THEN
                ITER=0
                EXIT
             END IF
          END IF
       END IF
       RHO=dum(1)
       ALPHA=dum(2)/RHO
       IF (RHO.EQ.0.0e0) EXIT
       IF (ITER.EQ.1 .AND. (.NOT.init)) THEN
          CALL SCOPY(N,S(1+2*N:3*N),IONE,X,IONE)
          CALL SSCAL(N,ALPHA,X,IONE)
          flop=flop+real(N)
       ELSE
          CALL SAXPY(N, ALPHA, S(1+2*N:3*N), IONE, X, IONE)
          flop=flop+real(2*N)
       END IF
       CALL SAXPY(N, -ALPHA, SD(1:n), IONE, f, IONE)
       dum0 = SNRM2(N,f,IONE)**2
       flop=flop+real(4*N)
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,dum0,1,MPI_REAL,  &
            MPI_SUM,ICOMM,ierr)
       RESID=dum0
       RESID=SQRT(RESID)
       RESID=RESID/BNORM
       IF (wff) THEN
          WRITE(iout,900) iter, resid*bnorm, resid
       END IF
    END DO
    IF( resid.GT.tol ) THEN
       IF (woo) THEN
          WRITE(iout,'()')
          WRITE(iout,950) iter
          WRITE(iout,951)
          WRITE(iout,'()')
       END IF
       iter=-iter
    ELSE
       IF (wff) THEN
          WRITE(iout,952) iter
          WRITE(iout,'()')
       END IF
    END IF
    RESID=RESID*BNORM/RESID0
    kstat(2,1)=ABS(iter)
    DEALLOCATE(S,SD)
    IF(ALLOCATED(fsc)) DEALLOCATE(fsc)
    RETURN
900 FORMAT('****  Iter=',i5,'        Resid=',e9.2,                 &
         '        Relat. res.=',e9.2)
940 FORMAT(i3,                                                     &
         '*SOLUTION: flexible conjugate gradient iterations (FCG(1))')
941 FORMAT('****     Rmk: solve system with the transpose of the input matrix')
943 FORMAT('****     AMG preconditioner with',i2,             &
         ' SOR pre- and post-smoothing sweeps')
944 FORMAT('****         at each level  (relaxation factor: omega=',f5.2,')')
945 FORMAT('****     AMG preconditioner with',i2,             &
         ' Gauss-Seidel pre- and post-smoothing sweeps')
946 FORMAT('****     AMG preconditioner with',i2,             &
         ' ILU(0) pre- plus',i2,' ILU(0) post-smoothing step(s)')
947 FORMAT('****         at each level')
950 FORMAT('**** !!!   WARNING:',I5,' ITERATIONS WERE')
951 FORMAT('**** INSUFFICIENT TO ACHIEVE THE DESIRED ACCURACY')
952 FORMAT('****  - Convergence reached in',I5,' iterations -')
998 FORMAT('**** The norm of the right hand side is zero:')
999 FORMAT('****     set x equal to the zero vector and exit')
  END SUBROUTINE sagmgpar_FlexCG
  SUBROUTINE sagmgpar_GCR(N,f,X,ITER,RESID,a,ja,ia,ijb)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER       :: N, ITER, ijb
    REAL(kind(0.0e0)) ::  f(N), X(N)
    INTEGER       :: ja(*), ia(N+1)
    REAL(kind(0.0e0)) :: a(*)
    INTEGER   :: MAXIT,i,irst,ierr,j,k,kk,idum(1),ii
    REAL(kind(0.0e0)) ::  ALPHA, BET0, tmp2
    REAL(kind(0.0e0)) ::  RESID, RESID2, BNORM2, RESID0
    REAL(kind(0.0e0)) ::  TOL, TOL2BNORM2, TOLT, dum0, t1, t2, rho, tmp1
    REAL(kind(0.0e0)) :: dum(3),xd,t
    REAL(kind(0.0e0)), ALLOCATABLE :: Sc(:),Siv(:),SiR(:)
    REAL(kind(0.0e0)), ALLOCATABLE :: Su(:),R(:),fsc(:)
    LOGICAL :: init,intmv,intmv2,scndp
    REAL(kind(0.0e0)), external :: SDOT
    REAL(kind(0.0e0)), external :: SNRM2
    INTEGER , parameter :: IONE=1
    INTEGER  :: itm1, m, info, itm
    CHARACTER*1 :: CHA
    intmv=.FALSE.
    IF ( (MOD(ijb,10).EQ.0 .OR. (ijb.GT.201)) .AND. nlev.GT.1) intmv=.TRUE.
    intmv2=intmv .OR. ijb.GT.201
    flop=0.0e0
    kstat=0
    kstat(1,MAX(2,nlev))=1
    init=.FALSE.
    IF (MOD(ijb,100).GE.10) init=.TRUE.
    IF (wfo) THEN
       WRITE(iout,938) IRANK,nrst
    END IF
    IF (wff) THEN
       IF (  trans) THEN
          WRITE(iout,941)
       END IF
       IF (nlev.GT.1) THEN
          IF (smoothtp.EQ.-1) THEN
             WRITE(iout,943) nsmooth
             WRITE(iout,944) omeg
          ELSE IF (smoothtp.EQ.1) THEN
             WRITE(iout,945) nsmooth
             WRITE(iout,947)
          ELSE
             WRITE(iout,946) nsmooth/2,(nsmooth+1)/2
             WRITE(iout,947)
          END IF
       END IF
    END IF
    TOL = RESID
    TOL2BNORM2 = TOL
    RESID2= 1.0e0
    MAXIT = ITER
    ITER = 0
    m=MIN(nrst,MAXIT)
    m=MAX(m,2)
    ALLOCATE (Su(N*m),Sc(N*m),SiR(((m+1)*m)/2),Siv(2*m+1),R(nwrkcum))
    dum(3) = SNRM2(N, f, IONE)**2
    flop=flop+real(2*N)
    IF (init) THEN
       IF (intmv) THEN
          CALL sagmgpar_sendrec(x,idum,-1,1)
          IF (smoothtp.EQ.1) THEN
             f(1:n)=f(1:n)-x(1:n)/dt(1)%a(1:n)
             CALL sagmgpar_csrmv(n,dt(1)%a,dt(1)%ja,n+1,dt(1)%il,x,1,f,-2,flop)
          ELSE
             CALL sagmgpar_mcsrmv(n,dt(1)%a,dt(1)%ja,n+1,dt(1)%il,x,1,f,-2,flop)
             f(1:n)=f(1:n)-x(1:n)/(dt(1)%p(1:n))
          END IF
          CALL sagmgpar_csrmv(n,dt(1)%a,dt(1)%ja,n+1,dt(1)%iu,x,1,f,-2,flop)
          IF (nneighi.GT.0) CALL MPI_WAITALL(nneighi,ireqi,MPI_STATUSES_IGNORE,ierr)
          CALL  sagmgpar_csrmv(n,dt(1)%a,dt(1)%ja,n+1,dt(1)%iext,buffi,n+1,f &
               ,-2,flop,dt(1)%lstin)
       ELSE
          IF (intmv2) THEN
             CALL sagmgpar_matv(N,x,Sc,dt(1)%a,dt(1)%ja,dt(1)%ia,trans &
             ,dt(1)%iext,dt(1)%lstin,1)
          ELSE
             CALL sagmgpar_matv(N,x,Sc,a,ja,ia,trans &
             ,iextL1,lstinL1,0)
          END IF
          f(1:n)=f(1:n)-Sc(1:N)
       END IF
       dum(2) = SNRM2(N, f, IONE)**2
       flop=flop+real(2*N)
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,dum(2:3),2,MPI_REAL,  &
            MPI_SUM,ICOMM,ierr)
       BNORM2=dum(3)
       RESID2=dum(2)
       IF (BNORM2.EQ.0.0e0) THEN
          IF(wff) THEN
             WRITE(iout,998)
             WRITE(iout,999)
          END IF
          X(1:N)=0.0e0
          DEALLOCATE (Su,Sc,R,Siv,SiR)
          IF(ALLOCATED(fsc)) DEALLOCATE(fsc)
          RETURN
       END IF
       TOL2BNORM2 = TOL*TOL*BNORM2
       IF (wff.AND. (MAXIT.LE.0 .OR. RESID2.LE.TOL2BNORM2)) THEN
          WRITE(iout,900) 0, 0,' ',SQRT(resid2),SQRT(resid2/bnorm2)
       END IF
       RESID0=SQRT(RESID2)
    END IF
    itm  = -1
    irst = 0
    CHA=' '
    DO WHILE ( ITER.LT.MAXIT .AND. RESID2.GT.TOL2BNORM2 )
       itm  = itm  + 1
       ITER = ITER + 1
       IF (itm.EQ.m) THEN
          CALL STPTRS('U','N','U',m,IONE,SiR,Siv,m,info)
          IF (irst.EQ.0 .AND. (.NOT.init)) THEN
             CALL SGEMV('N',N,m,1.0e0,Su,       &
                  N,Siv,IONE,0.0e0,X,IONE)
             flop=flop+real(2*m*N+m*(m+1))
          ELSE
             CALL SGEMV('N',N,m,1.0e0,Su,        &
                  N,Siv,IONE,1.0e0,X,IONE)
             flop=flop+real((2*m+1)*N+m*(m+1))
          END IF
          itm=0
          irst=irst+1
       END IF
       IF (ijb.GT.201) THEN
          CALL sagmgpar_prec_matv(1,f,Su(1+itm*N),Sc(1+itm*N),R,intmv)
       ELSE
          CALL sagmgpar_prec_matv(1,f,Su(1+itm*N),Sc(1+itm*N),R,intmv,a,ja,ia)
       END IF
       IF (.NOT.intmv) THEN
          IF (intmv2) THEN
             CALL sagmgpar_matv(N, Su(1+itm*N), Sc(1+itm*N)             &
                  , dt(1)%a, dt(1)%ja, dt(1)%ia, trans &
             ,dt(1)%iext,dt(1)%lstin,1)
          ELSE
             CALL sagmgpar_matv(N, Su(1+itm*N), Sc(1+itm*N)             &
                  , a, ja, ia, trans &
             ,iextL1,lstinL1,0)
          END IF
       END IF
       CHA=' '
       scndp=.FALSE.
       flop=flop+real(2*N)
10     CONTINUE
       IF (itm.GT.0)       &
          CALL SGEMV('C',N,itm,1.0e0,Sc,N,Sc(1+itm*N:N+itm*N), &
               IONE,0.0e0,Siv(m+2:m+itm+1),IONE)
       Siv(m+1)=SNRM2(N, Sc(1+itm*N:N+itm*N), IONE)**2
       Siv(m+itm+2)=SDOT(N, Sc(1+itm*N:(itm+1)*N), IONE, f, IONE)
       flop=flop+real((4+4*itm)*N+4*itm)
       kk=0
       IF (ITER.EQ.1) THEN
             Siv(m)=dum(3)
             kk=1
       END IF
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,Siv(m+1-kk:m+itm+2),itm+2+kk,MPI_REAL,   &
                MPI_SUM,ICOMM,ierr)
       IF (ITER.EQ.1) THEN
        IF (.NOT.init) THEN
          BNORM2=Siv(m)
          RESID2=BNORM2
          IF (BNORM2.EQ.0.0e0) THEN
             IF(wff) THEN
                WRITE(iout,998)
                WRITE(iout,999)
             END IF
             X(1:N)=0.0e0
             ITER=0
             DEALLOCATE (Su,Sc,R,Siv,SiR)
             IF(ALLOCATED(fsc)) DEALLOCATE(fsc)
             RETURN
          END IF
          TOL2BNORM2 = TOL*TOL*BNORM2
          RESID0=SQRT(RESID2)
          IF ( RESID2.LE.TOL2BNORM2 ) THEN
             ITER=0
             EXIT
          END IF
        END IF
        IF (wff) THEN
           WRITE(iout,900) 0, 0,' ',SQRT(resid2),SQRT(resid2/bnorm2)
        END IF
        TOLT = repsmach*RESID2
       END IF
       t1=Siv(m+1)
       t2=0.0e0
       DO i=0,itm-1
          t2=t2+Siv(m+2+i)*Siv(m+2+i)/SiR(1+i+(i*(i+1))/2)
          Siv(m+2+i)=-Siv(m+2+i)/SiR(1+i+(i*(i+1))/2)
       END DO
       rho=t1-t2
       alpha=Siv(m+itm+2)
       IF (itm.GT.0) THEN
          CALL SGEMV('N',N,itm,1.0e0,Sc,N,Siv(m+2:m+itm+1),IONE,    &
                 1.0e0,Sc(1+itm*N:N+itm*N),IONE)
          IF (scndp) THEN
             SiR(1+(itm*(itm+1))/2:itm+(itm*(itm+1))/2)=                   &
               SiR(1+(itm*(itm+1))/2:itm+(itm*(itm+1))/2)-Siv(m+2:m+itm+1)
          ELSE
             SiR(1+(itm*(itm+1))/2:itm+(itm*(itm+1))/2)=-Siv(m+2:m+itm+1)
             IF ( rho*(K_ItGramS**2).LE.t1 ) THEN
                scndp=.TRUE.
                GOTO 10
             END IF
          END IF
       END IF
       IF (rho.LE.0.0e0) THEN
          itm=itm-1
          EXIT
       END IF
       SiR(1+itm+(itm*(itm+1))/2)=rho
       bet0=alpha/rho
       Siv(1+itm)=bet0
       CALL SAXPY(N, -bet0, Sc(1+itm*N:(itm+1)*N), IONE, f, IONE)
       RESID2 = RESID2 - (ABS(alpha)**2)/rho
       IF (RESID2.LE.TOLT) THEN
          RESID2 = SNRM2(N,f,IONE)**2
          flop=flop+real(2*N)
          dum0 = RESID2
          CALL MPI_ALLREDUCE(MPI_IN_PLACE,dum0,1,MPI_REAL,  &
               MPI_SUM,ICOMM,ierr)
          RESID2=dum0
          TOLT = repsmach*RESID2
       END IF
       IF (wff)THEN
          WRITE(iout,900) iter,irst,CHA,SQRT(resid2),SQRT(resid2/bnorm2)
       END IF
    END DO
    IF (itm.GE.0) THEN
       itm1=itm+1
       CALL STPTRS('U','N','U',itm1, IONE,SiR,Siv,m,info)
       IF (irst.EQ.0 .AND. (.NOT.init)) THEN
          CALL SGEMV('N',N,itm1,1.0e0,Su,        &
               N,Siv,IONE,0.0e0,X,IONE)
          flop=flop+real(2*(itm+1)*N+(itm+1)*(itm+2))
       ELSE
          CALL SGEMV('N',N,itm1,1.0e0,Su,        &
               N,Siv,IONE,1.0e0,X,IONE)
          flop=flop+real((2*(itm+1)+1)*N+(itm+1)*(itm+2))
       END IF
    END IF
    RESID=SQRT(RESID2/BNORM2)
    IF( resid.GT.tol ) THEN
       IF (woo) THEN
          WRITE(iout,'()')
          WRITE(iout,950) iter
          WRITE(iout,951)
          WRITE(iout,'()')
       END IF
       iter=-iter
    ELSE
       IF (wff) THEN
          WRITE(iout,952) iter
          WRITE(iout,'()')
       END IF
    END IF
    RESID=RESID*SQRT(BNORM2)/RESID0
    DEALLOCATE (Su,Sc,R,Siv,SiR)
    IF(ALLOCATED(fsc)) DEALLOCATE(fsc)
    kstat(2,1)=ABS(iter)
    RETURN
900 FORMAT('****  Iter=',i5,' (',i2,' rest.)',A1,'       Resid=',e9.2,    &
         '        Relat. res.=',e9.2)
938 FORMAT(i3,'*SOLUTION: GCR iterations (GCR(',i2,'))')
941 FORMAT('****     Rmk: solve system with the transpose of the input matrix')
943 FORMAT('****     AMG preconditioner with',i2,             &
         ' SOR pre- and post-smoothing sweeps')
944 FORMAT('****         at each level  (relaxation factor: omega=',f5.2,')')
945 FORMAT('****     AMG preconditioner with',i2,             &
         ' Gauss-Seidel pre- and post-smoothing sweeps')
946 FORMAT('****     AMG preconditioner with',i2,             &
         ' ILU(0) pre- plus',i2,' ILU(0) post-smoothing step(s)')
947 FORMAT(  '****         at each level')
950 FORMAT('**** !!!   WARNING:',I5,' ITERATIONS WERE')
951 FORMAT('**** INSUFFICIENT TO ACHIEVE THE DESIRED ACCURACY')
952 FORMAT('****  - Convergence reached in',I5,' iterations -')
998 FORMAT('**** The norm of the right hand side is zero:')
999 FORMAT('****     set x equal to the zero vector and exit')
  END SUBROUTINE sagmgpar_GCR
  SUBROUTINE sagmgpar_matv(n, x, y, a, ja, ia, transpose,         &
       iext, lstin, l)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: n, ja(*), ia(n+1), i, kk, k1, k2, ier, idum(1)
    INTEGER, OPTIONAL :: iext(n), lstin(0:*), l
    REAL(kind(0.0e0)) :: y(n), a(*), t, xx
    REAL(kind(0.0e0)) :: x(n)
    LOGICAL :: transpose
    CALL sagmgpar_sendrec(x,idum,-1,l)
    IF (.NOT.transpose) THEN
       DO i = 1, n
          k1 = ia(i)
          xx = x(ja(k1))
          t = a(k1) * xx
          k2 = iext(i)-1
          DO kk = k1+1, k2
             xx = x(ja(kk))
             t = t + a(kk)*xx
          END DO
          y(i) = t
       END DO
    ELSE
       y(1:n)=0.0e0
       DO i = 1, n
          xx = x(i)
          DO kk = ia(i), ia(i+1)-1
             y(ja(kk)) = y(ja(kk)) + a(kk)*xx
          END DO
       END DO
    END IF
    IF (nneighi.GT.0) CALL MPI_WAITALL(nneighi,ireqi,MPI_STATUSES_IGNORE,ier)
    DO k1=1,lstin(0)
       i=lstin(k1)
       t = y(i)
       DO kk = iext(i), ia(i+1) -1
          xx = buffi(ja(kk)-n)
          t = t + a(kk)*xx
       END DO
       y(i) = t
    END DO
    flop = flop + real(2 *(ia(n+1)-ia(1))-n)
    RETURN
  END SUBROUTINE sagmgpar_matv
  RECURSIVE SUBROUTINE sagmgpar_prec_matv(l,B,X,AX,R,matv,a,ja,ia)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: l
    REAL(kind(0.0e0)), OPTIONAL ::  B(*), X(*), AX(*), R(*), A(*)
    LOGICAL, OPTIONAL :: matv
    INTEGER, OPTIONAL :: ja(*), ia(*)
    REAL(kind(0.0e0)) ::  dum(1)
    LOGICAL :: update
    INTEGER :: is,n,nnext,XN,BN,AXN,RN,idum(1),i,j,XNe,BNe,nnexte,nnext1
    INTEGER, PARAMETER :: smoothilu=max(smoothtype,1)-1
    INTEGER, PARAMETER :: npostsmooth=nsmooth+smoothilu*mod(nsmooth,2)
    INTEGER, PARAMETER ::  npresmooth=nsmooth-smoothilu*mod(nsmooth,2)
    n=nn(l)
    nnext=nn(l+1)
    IF (n.EQ.0) THEN
       IF (l+1.EQ.nlev) THEN
          IF (.NOT.allzero)                                      &
               CALL sagmgpar_MUMPSpar(nnext,dum,2)
       ELSE IF ( innermax(l+1).LE.1 ) THEN
          CALL sagmgpar_prec_matv(l+1)
          kstat(1,l+1)=MAX(kstat(1,l+1),1)
          kstat(2,l+1)=kstat(2,l+1)+1
       ELSE
          CALL sagmgpar_inner_fool(l+1)
       END IF
       RETURN
    END IF
    IF (l.EQ.nlev) THEN
       X(1:N)=B(1:N)
       IF (PRESENT(a)) THEN
          CALL sagmgpar_MUMPSpar(n,X,2,a,ja,ia)
       ELSE
          CALL sagmgpar_MUMPSpar(n,X,2)
       END IF
       RETURN
    END IF
    IF (npresmooth.EQ.0) THEN
       X(1:n)=0.0e0
    ELSE IF (npresmooth.EQ.1) THEN
       CALL sagmgpar_fwGS(l,B,X,AX,.FALSE.,R)
    ELSE
       update=.FALSE.
       DO is=2,npresmooth,2
          CALL sagmgpar_fwbwsolve1(l,B,X,AX,update,R)
          update=.TRUE.
       END DO
       IF (MOD(npresmooth,2).EQ.1) THEN
          CALL sagmgpar_fwGS(l,B,X,AX,.TRUE.,R)
       END IF
    END IF
    IF (nnext.GT.0) THEN
       XN=1
       BN=XN+nnext
       IF (l+1.EQ.nlev) BN=XN
       IF (npresmooth.EQ.0) THEN
         CALL sagmgpar_restag(N,nnext,B,R(BN),dt(l)%ind,flop)
       ELSE
         CALL sagmgpar_restag(N,nnext,AX,R(BN),dt(l)%ind,flop)
       END IF
       IF (l+1.EQ.nlev) THEN
          CALL sagmgpar_MUMPSpar(nnext,R(BN),2)
       ELSE IF ( innermax(l+1).LE.1 ) THEN
          AXN=BN+nnext
          RN=AXN+nnext
          CALL sagmgpar_prec_matv(l+1,R(BN),R(XN),R(AXN),R(RN),.FALSE.)
          kstat(1,l+1)=MAX(kstat(1,l+1),1)
          kstat(2,l+1)=kstat(2,l+1)+1
       ELSE
          CALL sagmgpar_inner_iter(nnext,R(XN),R(BN),l+1)
       END IF
       CALL sagmgpar_prolag(N,nnext,X,R(XN),dt(l)%ind,flop)
    ELSE
       IF (l+1.EQ.nlev) THEN
          IF (.NOT.allzero)                                      &
               CALL sagmgpar_MUMPSpar(nnext,dum,2)
       ELSE IF ( innermax(l+1).LE.1 ) THEN
          CALL sagmgpar_prec_matv(l+1,B,X,AX,R,.FALSE.)
          kstat(1,l+1)=MAX(kstat(1,l+1),1)
          kstat(2,l+1)=kstat(2,l+1)+1
       ELSE
          CALL sagmgpar_inner_fool(l+1)
       END IF
    END IF
    IF (npostsmooth.EQ.0) THEN
    ELSE IF (npostsmooth.EQ.1) THEN
       CALL sagmgpar_bwGS(l,B,X,AX,matv)
    ELSE
       IF (MOD(npostsmooth,2).EQ.1) THEN
          CALL sagmgpar_bwGS(l,B,X,AX,.FALSE.)
       END IF
       update=.FALSE.
       DO is=2,npostsmooth,2
          IF (is.GE.npostsmooth-1) update=matv
          CALL sagmgpar_fwbwsolve2(l,B,X,AX,update,R)
       END DO
    END IF
    RETURN
  END SUBROUTINE sagmgpar_prec_matv
  SUBROUTINE sagmgpar_fwGS(l,B,X,AX,update,R)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: l, n, idum(1), ier
    REAL(kind(0.0e0)) ::  B(*), X(*), AX(*), R(*)
    LOGICAL :: update
    n=nn(l)
    IF (update) THEN
       IF (smoothtp.EQ.1) THEN
         CALL sagmgpar_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%a   &
                              ,R,AX,1,flop)
       ELSE
         CALL sagmgpar_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%p   &
                              ,R,AX,1,flop)
       END IF
       CALL sagmgpar_sendrec(R,idum, -1,l)
       X(1:n)=X(1:n)+R(1:n)
       flop=flop+real(n)
       IF (smoothtp.EQ.1) THEN
         CALL  sagmgpar_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,R,1,AX,-1,flop)
       ELSE
         CALL  sagmgpar_mcsrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,R,1,AX,-1,flop)
       END IF
    ELSE
       IF (smoothtp.EQ.1) THEN
         CALL sagmgpar_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%a,  &
                               X,B,1,flop)
       ELSE
         CALL sagmgpar_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%p,  &
                               X,B,1,flop)
       END IF
       CALL sagmgpar_sendrec(X,idum,-1,l)
       IF (smoothtp.EQ.1) THEN
         CALL  sagmgpar_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,X,1,AX,-1,flop)
       ELSE
         CALL  sagmgpar_mcsrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,X,1,AX,-1,flop)
       END IF
    END IF
    IF (nneighi.GT.0) CALL MPI_WAITALL(nneighi,ireqi,MPI_STATUSES_IGNORE,ier)
    CALL  sagmgpar_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iext,buffi,n+1,AX,-2,flop &
                      ,dt(l)%lstin)
    RETURN
  END SUBROUTINE sagmgpar_fwGS
  SUBROUTINE sagmgpar_bwGS(l,B,X,AX,matv)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: l, n, idum(1), ier
    REAL(kind(0.0e0)) ::  B(*), X(*), AX(*)
    LOGICAL :: matv
    n=nn(l)
    CALL sagmgpar_sendrec(X,idum,-1,l)
    AX(1:n)=B(1:n)
    IF (smoothtp.EQ.1) THEN
      CALL  sagmgpar_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,-2,flop)
    ELSE
      CALL  sagmgpar_mcsrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,-2,flop)
    END IF
    IF (nneighi.GT.0) CALL MPI_WAITALL(nneighi,ireqi,MPI_STATUSES_IGNORE,ier)
    CALL  sagmgpar_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iext,buffi,n+1,AX,-2,flop &
         ,dt(l)%lstin)
    IF (smoothtp.EQ.1) THEN
      CALL sagmgpar_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,dt(l)%a,X,AX,1,flop)
    ELSE
      CALL sagmgpar_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,dt(l)%p,X,AX,1,flop)
    END IF
    IF (.NOT.matv) RETURN
    CALL sagmgpar_sendrec(X,idum,-1,l)
    IF (smoothtp.EQ.1) THEN
      CALL  sagmgpar_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,2,flop)
    ELSE
      CALL  sagmgpar_mcsrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,2,flop)
    END IF
    IF (nneighi.GT.0) CALL MPI_WAITALL(nneighi,ireqi,MPI_STATUSES_IGNORE,ier)
    CALL  sagmgpar_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iext,buffi,n+1,AX,2,flop &
         ,dt(l)%lstin)
    RETURN
  END SUBROUTINE sagmgpar_bwGS
  SUBROUTINE sagmgpar_fwbwsolve1(l,B,X,AX,update,R)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: l, n, idum(1), ier
    REAL(kind(0.0e0)) ::  B(*), X(*), AX(*), R(*)
    LOGICAL :: update
    n=nn(l)
    IF (.NOT.update) THEN
       IF (smoothtp.EQ.1) THEN
         CALL sagmgpar_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%a    &
                              ,R,B,-1,flop)
       ELSE
         CALL sagmgpar_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%p    &
                              ,R,B,-1,flop)
       END IF
       IF (smoothtp.EQ.1) THEN
         CALL sagmgpar_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,dt(l)%a    &
                              ,X,R,1,flop)
       ELSE
         CALL sagmgpar_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,dt(l)%p    &
                              ,X,R,1,flop)
       END IF
       CALL sagmgpar_sendrec(X,idum,-1,l)
       AX(1:n)=B(1:n)-R(1:n)
       flop=flop+real(n)
       IF (smoothtp.EQ.1) THEN
         CALL  sagmgpar_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,-2,flop)
       ELSE
         CALL  sagmgpar_mcsrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,-2,flop)
       END IF
       IF (nneighi.GT.0) CALL MPI_WAITALL(nneighi,ireqi,MPI_STATUSES_IGNORE,ier)
       CALL  sagmgpar_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iext,buffi,n+1,AX   &
                         ,-2,flop,dt(l)%lstin)
    ELSE
       IF (smoothtp.EQ.1) THEN
         CALL sagmgpar_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%a    &
                              ,R,AX,-1,flop)
       ELSE
         CALL sagmgpar_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%p    &
                              ,R,AX,-1,flop)
       END IF
       IF (smoothtp.EQ.1) THEN
         CALL sagmgpar_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,dt(l)%a    &
                              ,R(n+1),R,1,flop)
       ELSE
         CALL sagmgpar_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,dt(l)%p    &
                              ,R(n+1),R,1,flop)
       END IF
       X(1:n)=X(1:n)+R(n+1:2*n)
       CALL sagmgpar_sendrec(R(n+1),idum,-1,l)
       AX(1:n)=AX(1:n)-R(1:n)
       flop=flop+real(2*n)
       IF (smoothtp.EQ.1) THEN
         CALL  sagmgpar_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,R(n+1),1,AX  &
                           ,-2,flop)
       ELSE
         CALL  sagmgpar_mcsrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,R(n+1),1,AX  &
                           ,-2,flop)
       END IF
       IF (nneighi.GT.0) CALL MPI_WAITALL(nneighi,ireqi,MPI_STATUSES_IGNORE,ier)
       CALL  sagmgpar_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iext,buffi,n+1,AX   &
                         ,-2,flop,dt(l)%lstin)
    END IF
    RETURN
  END SUBROUTINE sagmgpar_fwbwsolve1
  SUBROUTINE sagmgpar_fwbwsolve2(l,B,X,AX,matv,R)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: l, n, idum(1), ier
    REAL(kind(0.0e0)) ::  B(*), X(*), AX(*), R(*)
    LOGICAL :: matv
    n=nn(l)
    CALL sagmgpar_sendrec(X,idum,-1,l)
    IF (smoothtp.EQ.1) THEN
      CALL  sagmgpar_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,X,1,AX,0,flop)
    ELSE
      CALL  sagmgpar_mcsrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,X,1,AX,0,flop)
    END IF
    R(1:n)=B(1:n)-AX(1:n)
    IF (nneighi.GT.0) CALL MPI_WAITALL(nneighi,ireqi,MPI_STATUSES_IGNORE,ier)
    CALL  sagmgpar_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iext,buffi,n+1,R,-2,flop &
                      ,dt(l)%lstin)
    IF (smoothtp.EQ.1) THEN
      CALL sagmgpar_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%a    &
                           ,R(n+1),R,1,flop)
    ELSE
      CALL sagmgpar_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%p    &
                           ,R(n+1),R,1,flop)
    END IF
    IF (matv) THEN
      IF (smoothtp.EQ.1) THEN
         AX(1:n)=AX(1:n)+R(n+1:2*n)/dt(l)%a(1:n)
      ELSE
         AX(1:n)=AX(1:n)+R(n+1:2*n)/dt(l)%p(1:n)
      END IF
      flop=flop+real(2*n)
    END IF
    R(n+1:2*n) = R(n+1:2*n) - X(1:n)
    IF (smoothtp.EQ.1) THEN
      CALL sagmgpar_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,dt(l)%a    &
                           ,R,R(n+1),-1,flop)
    ELSE
      CALL sagmgpar_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,dt(l)%p    &
                           ,R,R(n+1),-1,flop)
    END IF
    X(1:n)=X(1:n) + R(1:n)
    flop=flop+real(3*n)
    IF (.NOT.matv) RETURN
    CALL sagmgpar_sendrec(X,idum,-1,l)
    CALL  sagmgpar_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,2,flop)
    IF (smoothtp.NE.1) AX(1:n)=AX(1:n)+dt(l)%a(1:n)*R(1:n)
    IF (nneighi.GT.0) CALL MPI_WAITALL(nneighi,ireqi,MPI_STATUSES_IGNORE,ier)
    CALL  sagmgpar_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iext,buffi,n+1,AX,2,flop &
                      ,dt(l)%lstin)
    RETURN
  END SUBROUTINE sagmgpar_fwbwsolve2
  RECURSIVE SUBROUTINE sagmgpar_inner_iter(n,X,R,l)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER   :: N, ITER, l, ierr, k, i
    REAL(kind(0.0e0)) ::  RESID, BNORM, det
    REAL(kind(0.0e0)) :: X(N), R(MAX(N,1),*)
    REAL(kind(0.0e0)) :: alpha1,alpha2,bet0,bet1,bet2,rho1,rho2,gamm0,gamm1,zeta
    REAL(kind(0.0e0)) :: dum(5),y1,y2
    REAL(kind(0.0e0)), external :: SDOT
    REAL(kind(0.0e0)), external :: SNRM2
    INTEGER , parameter :: IONE=1
    ITER = 1
    CALL sagmgpar_prec_matv(l,R,X,R(1,2),R(1,3),.TRUE.)
    dum(1) = SDOT(N,X,IONE,R(1:N,2),IONE)
    dum(2) = SDOT(N,X,IONE,R(1:N,1),IONE)
    dum(3) = SNRM2(N, R(1:N,1), IONE)**2
    dum(4) = SNRM2(N, R(1:N,2), IONE)**2
    dum(5) = SDOT(N, R(1:N,1), IONE, R(1:N,2), IONE )
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,dum(1:5),5,MPI_REAL,  &
            MPI_SUM,ICOMM,ierr)
    IF (dum(3).EQ.0.0e0) THEN
       X(1:N)=0.0e0
       flop=flop+real(10*N)
       RETURN
    END IF
    BNORM = dum(3)
    rho1=dum(1)
    alpha1=dum(2)
    bet0=alpha1/rho1
    IF (ABS(dum(1)).LE.repsmach*ABS(dum(2))) THEN
       dum(1)=dum(4)
       dum(2)=dum(5)
       flop=flop+real(4*N)
       GOTO 110
    END IF
    RESID = BNORM + (ABS(bet0)**2)*dum(4)-2*bet0*dum(5)
    IF (RESID.LE.resi*resi*BNORM) THEN
       CALL SSCAL( N, bet0, X, IONE )
       flop=flop+real(11*N)
       kstat(1,l)=MAX(kstat(1,l),iter)
       kstat(2,l)=kstat(2,l)+iter
       RETURN
    END IF
    CALL SAXPY(N, -bet0, R(1:N,2), IONE, R(1:N,1), IONE)
    ITER = 2
    CALL sagmgpar_prec_matv(l,R,R(1,3),R(1,4),R(1,5),.TRUE.)
    dum(1) = SDOT(N,R(1:N,3),IONE,R(1:N,2),IONE)
    dum(2) = SDOT(N,R(1:N,3),IONE,R(1:N,1),IONE)
    dum(3) = SDOT(N,R(1:N,3),IONE,R(1:N,4),IONE)
    IF (.NOT.spd) dum(4) = SDOT(N,X,IONE,R(1:N,4),IONE)
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,dum(1:4),4,MPI_REAL,  &
            MPI_SUM,ICOMM,ierr)
    gamm0 =dum(1)
    alpha2=dum(2)
    rho2  =dum(3)
    IF (spd) THEN
       gamm1 = gamm0
    ELSE
       gamm1 = dum(4)
    END IF
    bet1=rho2-gamm0*gamm1/rho1
    bet2=(alpha1-alpha2*gamm1/bet1)/rho1
    zeta=alpha2/bet1
    IF (ABS(bet1).LE.repsmach*ABS(alpha2) .OR.   &
        ABS(bet2)*repsmach.GE.1.0e0)            THEN
       flop=flop+real(6*N)
       IF (.NOT.spd) flop=flop+real(2*N)
      GOTO 200
    END IF
    CALL SSCAL(N, bet2, X, IONE)
    CALL SAXPY(N, zeta, R(1:N,3), IONE, X, IONE)
    kstat(1,l)=MAX(kstat(1,l),iter)
    kstat(2,l)=kstat(2,l)+iter
    flop=flop+real(21*N)
    IF (.NOT.spd) flop=flop+real(2*N)
    RETURN
100 CONTINUE
    dum(1)=SNRM2(N, R(1:N,2), IONE)**2
    dum(2)=SDOT(N, R(1:N,2), IONE, R(1:N,1), IONE )
    dum(3) = SNRM2(N, R(1:N,1), IONE)**2
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,dum(1:3),3,MPI_REAL,  &
            MPI_SUM,ICOMM,ierr)
    IF (dum(3).EQ.0.0e0) THEN
       X(1:N)=0.0e0
       flop=flop+real(10*N)
       RETURN
    END IF
    BNORM = dum(3)
110 CONTINUE
    rho1=dum(1)
    alpha1=dum(2)
    bet0=alpha1/rho1
    RESID=BNORM-alpha1*bet0
    IF (RESID.LE.resi*resi*BNORM) THEN
       CALL SSCAL( N, bet0, X, IONE )
       flop=flop+real(7*N)
       kstat(1,l)=MAX(kstat(1,l),iter)
       kstat(2,l)=kstat(2,l)+iter
       RETURN
    END IF
    CALL SAXPY(N, -bet0, R(1:N,2), IONE, R(1:N,1), IONE)
    ITER = 2
    CALL sagmgpar_prec_matv(l,R,R(1,3),R(1,4),R(1,5),.TRUE.)
    dum(1) = SDOT(N,R(1:N,4),IONE,R(1:N,2),IONE)
    dum(2) = SDOT(N,R(1:N,4),IONE,R(1:N,1),IONE)
    dum(3) = SNRM2(N,R(1:N,4),IONE)**2
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,dum(1:3),3,MPI_REAL,  &
            MPI_SUM,ICOMM,ierr)
    gamm0 =dum(1)
    alpha2=dum(2)
    rho2  =dum(3)
    gamm1 = gamm0
    bet1=rho2-gamm0*gamm1/rho1
    bet2=(alpha1-alpha2*gamm1/bet1)/rho1
    zeta=alpha2/bet1
    CALL SSCAL(N, bet2, X, IONE)
    CALL SAXPY(N, zeta, R(1:N,3), IONE, X, IONE)
    kstat(1,l)=MAX(kstat(1,l),iter)
    kstat(2,l)=kstat(2,l)+iter
    flop=flop+real(19*N)
    RETURN
200 CONTINUE
    dum(1) = SNRM2(N,R(1:N,2),IONE)**2
    dum(2) = SDOT(N,R(1:N,4),IONE,R(1:N,2),IONE)
    dum(3) = SNRM2(N,R(1:N,4),IONE)**2
    dum(4) = SDOT(N,R(1:N,2),IONE,R(1:N,1),IONE)
    dum(5) = SDOT(N,R(1:N,4),IONE,R(1:N,1),IONE)
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,dum(1:5),5,MPI_REAL,  &
            MPI_SUM,ICOMM,ierr)
     dum(4) = dum(4)+bet0*dum(1)
     dum(5) = dum(5)+bet0*dum(2)
     det = dum(1)*dum(3)-dum(2)*dum(2)
     y1  = (dum(3)*dum(4)-dum(2)*dum(5))/det
     y2  = (-dum(2)*dum(4)+dum(1)*dum(5))/det
     CALL SSCAL( N, y1, X, IONE )
     CALL SAXPY( N, y2, R(1:N,3), IONE, X, IONE )
    kstat(1,l)=MAX(kstat(1,l),iter)
    kstat(2,l)=kstat(2,l)+iter
    flop=flop+real(25*N)
    RETURN
  END SUBROUTINE sagmgpar_inner_iter
  SUBROUTINE sagmgpar_prolag(n, nc, V, B, ind, flop)
    IMPLICIT NONE
    INTEGER :: n, nc, ind(n), k, i
    REAL(kind(0.0e0)) :: V(n), B(nc)
    REAL(kind(0.0e0)) :: flop
    IF (nc.LE.0) RETURN
    DO i = 1, n
       k = ind(i)
       IF (k.GT.0) THEN
          V(i) = V(i)+B(k)
       END IF
    END DO
    flop = flop + real(n)
    RETURN
  END SUBROUTINE sagmgpar_prolag
  SUBROUTINE sagmgpar_restag(n, nc, V, B, ind, flop)
    IMPLICIT NONE
    INTEGER :: n, nc, ind(n), k, i
    REAL(kind(0.0e0)) :: V(n), B(nc)
    REAL(kind(0.0e0)) :: flop
    IF (nc.LE.0) RETURN
    B(1:nc)=0.0e0
    DO i = 1, n
       k = ind(i)
       IF (k.GT.0) B(k) = B(k) + V(i)
    END DO
    flop = flop + real(n)
    RETURN
  END SUBROUTINE sagmgpar_restag
  RECURSIVE SUBROUTINE sagmgpar_inner_fool(l)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: ITER, l, ierr
    REAL(kind(0.0e0)) :: dum(5)
    REAL(kind(0.0e0)) :: RESID, BNORM, det
    REAL(kind(0.0e0)) :: alpha1,alpha2,bet0,bet1,bet2,rho1,rho2,gamm0,gamm1,zeta
    ITER = 1
    CALL sagmgpar_prec_matv(l)
    dum(1:5) = 0.0e0
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,dum(1:5),5,MPI_REAL,  &
            MPI_SUM,ICOMM,ierr)
    IF (dum(3).EQ.0.0e0) RETURN
    BNORM = dum(3)
    rho1=dum(1)
    alpha1=dum(2)
    bet0=alpha1/rho1
    IF (ABS(dum(1)).LE.repsmach*ABS(dum(2))) THEN
       dum(1)=dum(4)
       dum(2)=dum(5)
       GOTO 110
    END IF
    RESID = BNORM + bet0*bet0*dum(4)-2*bet0*dum(5)
    IF (RESID.LE.resi*resi*BNORM) THEN
       kstat(1,l)=MAX(kstat(1,l),iter)
       kstat(2,l)=kstat(2,l)+iter
       RETURN
    END IF
    ITER = 2
    CALL sagmgpar_prec_matv(l)
    dum(1:4) =  0.0e0
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,dum(1:4),4,MPI_REAL,  &
            MPI_SUM,ICOMM,ierr)
    gamm0 =dum(1)
    alpha2=dum(2)
    rho2  =dum(3)
    IF (spd) THEN
       gamm1 = gamm0
    ELSE
       gamm1 = dum(4)
    END IF
    bet1=rho2-gamm0*gamm1/rho1
    bet2=(alpha1-alpha2*gamm1/bet1)/rho1
    zeta=alpha2/bet1
    IF (ABS(bet1).LE.repsmach*ABS(alpha2) .OR.   &
        ABS(bet2)*repsmach.GE.1.0e0)            THEN
      GOTO 200
    END IF
    kstat(1,l)=MAX(kstat(1,l),iter)
    kstat(2,l)=kstat(2,l)+iter
    RETURN
100 CONTINUE
    dum(1:3)=0.0e0
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,dum(1:3),3,MPI_REAL,  &
            MPI_SUM,ICOMM,ierr)
    IF (dum(3).EQ.0.0e0)  RETURN
    BNORM = dum(3)
110 CONTINUE
    rho1=dum(1)
    alpha1=dum(2)
    bet0=alpha1/rho1
    RESID=BNORM-alpha1*bet0
    IF (RESID.LE.resi*resi*BNORM) THEN
       kstat(1,l)=MAX(kstat(1,l),iter)
       kstat(2,l)=kstat(2,l)+iter
       RETURN
    END IF
    ITER = 2
    CALL sagmgpar_prec_matv(l)
    dum(1:3)=0.0e0
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,dum(1:3),3,MPI_REAL,  &
            MPI_SUM,ICOMM,ierr)
    kstat(1,l)=MAX(kstat(1,l),iter)
    kstat(2,l)=kstat(2,l)+iter
    RETURN
200 CONTINUE
    dum(1:5)=0.0e0
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,dum(1:5),5,MPI_REAL,  &
            MPI_SUM,ICOMM,ierr)
    kstat(1,l)=MAX(kstat(1,l),iter)
    kstat(2,l)=kstat(2,l)+iter
    RETURN
  END SUBROUTINE sagmgpar_inner_fool
  SUBROUTINE sagmgpar_setupL1(n,a,ja,ia)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: n
    INTEGER :: ja(*),ia(n+1)
    REAL(kind(0.0e0)) :: a(*)
    INTEGER :: nc,ierr,i,j,k,nz,nza,icmm
    LOGICAL :: slcoarse
    INTEGER, POINTER, DIMENSION(:) :: jap
    REAL(kind(0.0e0)), POINTER, DIMENSION(:) :: ap
    REAL(kind(0.0e0)) :: eta,dum(3)
    CHARACTER(len=13) :: prtint
    REAL (kind(0.0e0)) :: fff(1)
    nn(1)=n
    nlc(1)=n
    nlc(2)=ia(n+1)-ia(1)
    nza=nlc(2)
    ngl=nlc
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,ngl,2,MPI_REAL,  &
            MPI_SUM,ICOMM,ierr)
    maxcoarset=FLOOR(maxcoarset*((ngl(1)/NPROC)**(1.0e0/3)))
    maxcoarseslowt=FLOOR(maxcoarseslowt*((ngl(1)/NPROC)**(1.0e0/3)))
    IF ( 0.EQ.nstep  .OR. 1.EQ.maxlev .OR. ngl(1).LE.maxcoarset ) nlev=1
       coasing=.FALSE.
       wcplex=1.0e0
       wcplex(4)=0.0e0
       nlctot=nlc
       ngltot=ngl
       ngl1=ngl
       nlc1=nlc
       fracnz(1)=1.0e0
       allzero=.FALSE.
       nglp=0.0e0
       IF (wfo) THEN
          IF(wff) THEN
             IF (ngl(1).GT.9.9e10) THEN
                WRITE(prtint(1:12),'(1pe12.5)') ngl(1)
             ELSE
                WRITE(prtint,'(f13.0)') ngl(1)
             END IF
             WRITE(iout,'()')
             WRITE(iout,905) prtint(1:12)
          END IF
          IF (n.GT.9.9e10) THEN
             WRITE(prtint(1:12),'(1pe12.5)') real(n)
          ELSE
             WRITE(prtint(1:12),'(i12)') n
          END IF
          WRITE(iout,906) IRANK,prtint(1:12)
          IF(wff) THEN
             IF (ngl(2).GT.9.9e10) THEN
                WRITE(prtint(1:12),'(1pe12.5)') ngl(2)
             ELSE
                WRITE(prtint,'(f13.0)') ngl(2)
             END IF
             WRITE(iout,907)  prtint(1:12),ngl(2)/ngl(1)
          END IF
          IF (nlc(2).GT.9.9e10) THEN
             WRITE(prtint(1:12),'(1pe12.5)') real(nlc(2))
          ELSE
             WRITE(prtint(1:12),'(i12)') nlc(2)
          END IF
          WRITE(iout,908) IRANK,prtint(1:12),real(nlc(2))/real(n)
          WRITE(iout,'()')
       END IF
    nlcp=nlc
    nglp=ngl
    IF (1.NE.nlev) THEN
       CALL sagmgpar_aggregation(1,n,a,ja,ia,nc)
       CALL sagmgpar_setup(2,nc)
    ELSE
       DEALLOCATE(dt(1)%idiag)
       CALL sagmgpar_MUMPSpar(n,fff,1,a,ja,ia)
       nwrkcum=1
       IF (ngl(1).GT.maxcoarset) wcplex(4)=log10(real(ngl(1))/maxcoarset)
       RETURN
    END IF
       nz=ia(n+1)-ia(1)
       IF (transint) THEN
          ALLOCATE(ap(nz),jap(nz-n),dt(1)%il(n+1),dt(1)%iu(n+1))
          CALL sagmgpar_csrdluT(n,a,ja,ia,dt(1)%idiag              &
               ,ap,jap,dt(1)%il,dt(1)%iu)
          DEALLOCATE(dt(1)%idiag)
       ELSE
          ALLOCATE(ap(nz),jap(nz-n),dt(1)%iu(n+1))
          dt(1)%il => dt(1)%idiag
          NULLIFY(dt(1)%iext)
          ALLOCATE(dt(1)%iext(n+1))
          CALL sagmgpar_csrdlu(n,a,ja,ia,dt(1)%idiag              &
               ,ap,jap,dt(1)%iu &
               ,iextL1,dt(1)%iext)
          NULLIFY(dt(1)%idiag)
       END IF
       dt(1)%a  => ap
       dt(1)%ja => jap
       NULLIFY(ap,jap)
       innermax(nlev)=0
       innermax(1)=1
       eta=xsi/((1-xsi)*(cplxmax-1))
       icmm=1
       DO i=2,nlev-1
          innermax(i)=min(2,floor(xsi**(i-1)/(eta*fracnz(i)*icmm)))
          IF (nlev-i.LE.nlvcyc .AND. i.GT.2) innermax(i)=1
          icmm=icmm*innermax(i)
          wcplex(2)=wcplex(2)+icmm*fracnz(i)
          wcplex(3)=wcplex(3)+(2**(i-1))*fracnz(i)
       END DO
       wcplex(2)=wcplex(2)+(2**(nlev-1))*fracnz(nlev)
       wcplex(3)=wcplex(3)+(2**(nlev-1))*fracnz(nlev)
       IF (nsmooth.GT.1 .OR. smoothtype.GT.1)  THEN
          nwrkcum=2*MAX(nn(nlev-1),1)
       ELSE
          nwrkcum=MAX(nn(nlev),1)
       END IF
       DO i=nlev-2,1,-1
          nwrkcum=nwrkcum+3*MAX(nn(i+1),1)
          IF (innermax(i+1).GT.1) nwrkcum=nwrkcum+2*MAX(nn(i+1),1)
          IF (nsmooth.GT.1 .OR. smoothtype.GT.1)  nwrkcum=max(2*nn(i),nwrkcum)
       END DO
             IF (wfo) THEN
                WRITE(iout,'()')
             END IF
             IF (wff) THEN
                WRITE(iout,950)  ngltot(1)/ngl1(1)
                WRITE(iout,951)  ngltot(2)/ngl1(2)
             END IF
             IF (wfo) THEN
                WRITE(iout,952)  IRANK,nlctot(1)/real(n)
                WRITE(iout,953)  IRANK,nlctot(2)/real(nza)
             END IF
       IF (wff) THEN
          WRITE(iout,956) wcplex(3)
          WRITE(iout,957) wcplex(2)
          WRITE(iout,'()')
       END IF
    RETURN
905       FORMAT('****',' Global number of unknowns:', A12)
906       FORMAT(i3,'*','      Number of local rows:', A12)
907       FORMAT('****',' Global number of nonzeros:',A12,                 &
               ' (per row:',f7.2,')')
908       FORMAT(i3,'*','    Nonzeros in local rows:',A12,                 &
               ' (per row:',f7.2,')')
950 FORMAT('****','           Global grid complexity:',f9.2)
951 FORMAT('****','       Global Operator complexity:',f9.2)
952 FORMAT(i3,'*','            Local grid complexity:',f9.2)
953 FORMAT(i3,'*','        Local Operator complexity:',f9.2)
956 FORMAT('****','  Theoretical Weighted complexity:',f9.2, &
                ' (K-cycle at each level)' )
957 FORMAT('****','    Effective Weighted complexity:',f9.2, &
                ' (V-cycle enforced where needed)' )
958 FORMAT('****','              Weighted complexity:',f9.2)
  END SUBROUTINE sagmgpar_setupL1
  RECURSIVE SUBROUTINE sagmgpar_setup(l,n)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: l,n
    INTEGER :: nc,ierr,i,j,k,nz
    LOGICAL :: slcoarse
    INTEGER, POINTER, DIMENSION(:) :: jap
    REAL(kind(0.0e0)), POINTER, DIMENSION(:) :: ap
    LOGICAL, SAVE :: slowcoarse
    REAL(kind(0.0e0)) :: fw,eta,dum(3),aa
    CHARACTER(len=13) :: prtint
    REAL (kind(0.0e0)) :: fff(1)
    nn(l)=n
    nlc(1)=n
    nlc(2)=0
    IF (n.GT.0) THEN
          nlc(2)=dt(l)%ia(n+1)-dt(l)%ia(1)
    END IF
    ngl=nlc
    IF (l.EQ.2) slowcoarse=.FALSE.
       CALL MPI_ALLREDUCE(MPI_IN_PLACE,ngl,2,MPI_REAL,  &
            MPI_SUM,ICOMM,ierr)
    slcoarse = 2*nglp(1).LT.3*ngl(1) .AND. 2*nglp(2).LT.3*ngl(2)
    IF( l.EQ.nstep+1  .OR. l.EQ.maxlev                        &
         .OR. ( ngl(1).LE.maxcoarset)                         &
         .OR. ( nglp(1).LT.2*ngl(1) .AND. nglp(2).LT.2*ngl(2)   &
         .AND. ngl(1).LE.maxcoarseslowt )  &
         .OR. ( slowcoarse .AND. slcoarse )                   &
         .OR. nglp(1).EQ.ngl(1) )                       THEN
       nlev=l
    END IF
    IF (l.EQ.2) wcplex(1)=ngl1(2)/ngl(2)
    slowcoarse=slcoarse
    fracnz(l)=ngl(2)/ngl1(2)
    nlctot=nlctot+nlc
    ngltot=ngltot+ngl
    IF (wfo) THEN
       WRITE(iout,'()')
       WRITE(iout,914) IRANK,l
       IF(wff) THEN
          IF (ngl(1).GT.0.5e0) THEN
             IF (ngl(1).GT.9.9e10) THEN
                WRITE(prtint(1:12),'(1pe12.5)') ngl(1)
             ELSE
                WRITE(prtint,'(f13.0)') ngl(1)
             END IF
             WRITE(iout,930) prtint(1:12),nglp(1)/ngl(1)
          ELSE
             WRITE(iout,934) 0
          END IF
       END IF
       IF (n.GT.0) THEN
          IF (n.GT.9.9e10) THEN
             WRITE(prtint(1:12),'(1pe12.5)') real(n)
          ELSE
             WRITE(prtint(1:12),'(i12)') n
          END IF
          WRITE(iout,931) IRANK,prtint,real(nlcp(1))/real(n)
       ELSE
          WRITE(iout,935) IRANK,n
       END IF
       IF(wff .AND. ngl(1).GT.0.5e0) THEN
          IF (ngl(2).GT.9.9e10) THEN
             WRITE(prtint(1:12),'(1pe12.5)') ngl(2)
          ELSE
             WRITE(prtint,'(f13.0)') ngl(2)
          END IF
          WRITE(iout,932)  prtint(1:12),ngl(2)/ngl(1),nglp(2)/ngl(2)
       END IF
       IF (n.GT.0) THEN
          IF (nlc(2).GT.9.9e10) THEN
             WRITE(prtint(1:12),'(1pe12.5)') real(nlc(2))
          ELSE
             WRITE(prtint(1:12),'(i12)') nlc(2)
          END IF
          WRITE(iout,933) IRANK,prtint(1:12),real(nlc(2))/real(n),  &
               real(nlcp(2))/real(nlc(2))
       END IF
    END IF
    nlcp=nlc
    nglp=ngl
    IF (n.EQ.0) THEN
       nc=0
       allzero=.TRUE.
       IF (ngl(1).GT.0.5e0) allzero=.FALSE.
       IF (nlev /= l) THEN
          CALL sagmgpar_setup(l+1,nc)
       ELSE IF (.NOT.allzero) THEN
          CALL sagmgpar_MUMPSpar(n,fff,1)
       END IF
       RETURN
    END IF
    IF (l /= nlev) THEN
       CALL sagmgpar_aggregation(l,n,dt(l)%a,dt(l)%ja,dt(l)%ia,nc)
       CALL sagmgpar_setup(l+1,nc)
    ELSE
       IF (ngl(1).GT.maxcoarset) wcplex(4)=log10(real(ngl(1))/maxcoarset)
          DEALLOCATE(dt(l)%idiag)
          CALL sagmgpar_MUMPSpar(n,fff,1,dt(l)%a,dt(l)%ja,dt(l)%ia)
          IF (coasing .AND. wfo) THEN
             WRITE(iout,912) IRANK
          END IF
          RETURN
    END IF
    IF (n.EQ.0) RETURN
    nz=dt(l)%ia(n+1)-dt(l)%ia(1)
    IF (transint) THEN
       ALLOCATE(ap(nz),jap(nz-n),dt(l)%il(n+1),dt(l)%iu(n+1))
       CALL sagmgpar_csrdluT(n,dt(l)%a,dt(l)%ja,dt(l)%ia,dt(l)%idiag &
            ,ap,jap,dt(l)%il,dt(l)%iu)
       DEALLOCATE(dt(l)%a,dt(l)%ja)
       DEALLOCATE(dt(l)%idiag,dt(l)%ia)
    ELSE
       ALLOCATE(ap(nz),jap(nz-n))
       dt(l)%iu => dt(l)%ia
       dt(l)%il => dt(l)%idiag
       CALL sagmgpar_csrdlu(n,dt(l)%a,dt(l)%ja,dt(l)%ia,dt(l)%idiag &
            ,ap,jap &
       ,dt(l)%iu,dt(l)%iext)
       DEALLOCATE(dt(l)%a,dt(l)%ja)
       NULLIFY(dt(l)%idiag,dt(l)%ia)
    END IF
    dt(l)%a  => ap
    dt(l)%ja => jap
    NULLIFY(ap,jap)
    RETURN
912 FORMAT(i3,'*','        Warning: coarsest grid matrix treated as singular')
914 FORMAT(i3,'*','                      Level:',I12)
930 FORMAT('****',' Global number of variables:',A12,              &
         '          (reduction ratio:',f5.2,')')
931 FORMAT(i3,'*','       Number of local rows:',A12,              &
         '          (reduction ratio:',f5.2,')')
932 FORMAT('****','  Global number of nonzeros:',A12,              &
         ' (per row:',f4.1,  &
         '; red. ratio:',f5.2,')')
933 FORMAT(i3,'*','     Nonzeros in local rows:',A12,              &
         ' (per row:',f4.1,  &
         '; red. ratio:',f5.2,')')
934 FORMAT('****',' Global number of variables:',i12)
935 FORMAT(i3,'*','       Number of local rows:',i12)
  END SUBROUTINE sagmgpar_setup
  SUBROUTINE sagmgpar_smoothsetup
    USE sagmgpar_mem
    IMPLICIT NONE
    REAL(kind(0.0e0)) :: unmominv
    INTEGER :: l,i,llev
    smoothtp=smoothtype
    IF (smoothtype.EQ.0) THEN
       IF (spd) THEN
          smoothtp=1
       ELSE IF (omeg.LE.0.97) THEN
          smoothtp=-1
       ELSE
          smoothtp=1
       END IF
    END IF
    llev=nlev-1
    IF (smoothtp.EQ.-1) THEN
       IF (smoothtype.EQ.-1) omeg=omega
       unmominv=1.0e0-1.0e0/omeg
       DO l=1,llev
          IF (nn(l).GT.0) THEN
            ALLOCATE(dt(l)%p(nn(l)))
            DO i=1,nn(l)
               dt(l)%p(i)=omeg/dt(l)%a(i)
               dt(l)%a(i)=unmominv*dt(l)%a(i)
            END DO
          END IF
       END DO
    ELSE IF (smoothtp.EQ.1) THEN
       DO l=1,llev
          IF (nn(l).GT.0) THEN
            DO i=1,nn(l)
               dt(l)%a(i)=1.0e0/dt(l)%a(i)
            END DO
          END IF
       END DO
    ELSE IF (smoothtp.EQ.2) THEN
       DO l=1,llev
          IF (nn(l).GT.0) THEN
            ALLOCATE(dt(l)%p(nn(l)))
            CALL sagmgpar_ilu0(nn(l),dt(l)%a,dt(l)%ja,dt(l)%iu,dt(l)%il,dt(l)%p)
          END IF
       END DO
    END IF
    RETURN
  END SUBROUTINE sagmgpar_smoothsetup
  SUBROUTINE sagmgpar_ilu0(n,a,ja,iu,il,p)
    IMPLICIT NONE
    INTEGER :: n, ja(n+1:*), il(n+1), iu(n+1)
    REAL(kind(0.0e0)) :: a(*), p(n), t
    INTEGER :: i,j,k,l
    DO i=1,n
      t=a(i)
      DO j=il(i),il(i+1)-1
         k=ja(j)
         DO l=iu(k),iu(k+1)-1
            IF (ja(l).EQ.i) THEN
               t=t-a(j)*a(l)*p(k)
               EXIT
            END IF
         END DO
      END DO
      a(i)=a(i)-t
      p(i)=1.0e0/t
    END DO
    RETURN
  END SUBROUTINE sagmgpar_ilu0
  SUBROUTINE sagmgpar_omegaest(n,a,ja,ia,idiag,si,w,iext)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n)
    INTEGER, OPTIONAL :: iext(*)
    REAL(kind(0.0e0)) :: a(*)
    REAL(kind(0.0e0)) :: si(n),w(2*n),t,tm
    INTEGER :: i,j,j1,j2,jd,ll,l2,kl,kk
    INTEGER :: ifirst,ilast,kb
    w(1:2*n)=0.0e0
    ifirst=1
    ilast=n
     DO i=ilast,ifirst,-1
       j1 =ia(i)
       jd=idiag(i)
       j2=iext(i)-1
       DO kk=j1,jd-1
        j=ja(kk)
          l2=iext(j)-1
          kl=0
          DO ll=idiag(j)+1,l2
            IF (ja(ll).EQ.i) THEN
              kl=ll
              EXIT
            END IF
          END DO
          IF (kl.EQ.0) THEN
             w(j)=w(j)+ABS(a(kk))
          ELSE
             w(j)=w(j)+ABS(a(kk)-a(kl))-ABS(a(kl))
          END IF
       END DO
       t=w(i)
       DO kk=jd+1,j2
          t=t+ABS(a(kk))
       END DO
       w(i)=t/(max(a(idiag(i)),si(i)))
     END DO
     tm=0.0e0
    ifirst=1
    ilast=n
     DO i=ifirst,ilast
       j1 =ia(i)
       jd=idiag(i)
       j2=iext(i)-1
       DO kk=jd+1,j2
        j=ja(kk)
          kl=0
          DO ll=ia(j),idiag(j)-1
            IF (ja(ll).EQ.i) THEN
               kl=ll
               EXIT
            END IF
          END DO
          IF (kl.EQ.0) THEN
             w(n+j)=w(n+j)+ABS(a(kk))*w(i)
          ELSE
             w(n+j)=w(n+j)+(ABS(a(kk)-a(kl))-ABS(a(kl)))*w(i)
          END IF
       END DO
       t=w(n+i)
       DO kk=j1,jd-1
          t=t+ABS(a(kk))*w(ja(kk))
       END DO
       t=t/(max(a(idiag(i)),si(i)))
       tm=max(t,tm)
     END DO
    IF (tm.LE.1.0e0) THEN
      omeg=1.0e0
    ELSE
      omeg=2/(1.0e0+SQRT(tm))
    END IF
    RETURN
  END SUBROUTINE sagmgpar_omegaest
  SUBROUTINE sagmgpar_aggregation(l,n,a,ja,ia,nc)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: l,n,nc
    INTEGER :: ja(*),ia(n+1)
    REAL (kind(0.0e0)) :: a(*), dum
    INTEGER :: ier,i,j,k,maxdg,np,kpass,nzc,m1,ndd,nzp,isize,nddp,npass1,nz,i0
    LOGICAL :: skipass
    INTEGER, POINTER, DIMENSION(:) :: jan,ian,idiagn,iextn,ind2,lcg,lcgn
    REAL(kind(0.0e0)), POINTER, DIMENSION(:) :: an
    REAL(kind(0.0e0)), POINTER, DIMENSION(:) :: sinn
    INTEGER, POINTER, DIMENSION(:) :: jap,iap,idiagp,iextp,lcg1
    REAL(kind(0.0e0)), POINTER, DIMENSION(:) :: ap
    REAL(kind(0.0e0)), POINTER, DIMENSION(:) :: sip
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ldd,iw,iperm,riperm
    REAL(kind(0.0e0)), ALLOCATABLE, DIMENSION(:) :: si1,w
    REAL(kind(0.0e0)), ALLOCATABLE, DIMENSION(:) :: wc
    IF (l.EQ.1) THEN
       IF (wfo) THEN
          WRITE (iout,901) IRANK
       END IF
       IF (wff) THEN
          IF (spd) THEN
             WRITE (iout,906)
          ELSE IF (  transint) THEN
             WRITE (iout,908)
          END IF
          IF (.not.spd) then
             WRITE (iout,902) 'Jacobi',kaptg_dampJac,checkddJ
          ELSE
             WRITE (iout,902) 'BlockD',kaptg_blocdia,checkddB
          END IF
          IF (checkdd.LT.0) THEN
             WRITE (iout,904)
          END IF
          WRITE(iout,903) npass,targetcoarsefac
          WRITE (iout,905) trspos
       END IF
    END IF
    IF (l.EQ.1) THEN
       ALLOCATE(si1(n),ind2(n),iperm(n),riperm(n))
       CALL sagmgpar_setCMK(n,ja,ia,dt(l)%idiag,riperm,iperm,dt(l)%iext)
    ELSE
       ALLOCATE(si1(n),ind2(n),iperm(n))
       iperm(1:n)=1
    END IF
    CALL sagmgpar_prepareagg(n,a,ja,ia,dt(l)%idiag,ind2,iperm,si1,ndd,l &
            ,dt(l)%iext)
    IF (ndd.EQ.n) THEN
       nc=0
       nzc=0
       DEALLOCATE(si1,iperm,ind2)
       IF (l.EQ.1) THEN
          DEALLOCATE(riperm)
          IF (smoothtype.EQ.0) omeg=1.0e0
       END IF
       GOTO 999
    END IF
    IF (smoothtype.EQ.0 .AND. l.EQ.1 .AND. (.NOT.spd)) THEN
        ALLOCATE(w(2*n))
        CALL sagmgpar_omegaest(n,a,ja,ia,dt(l)%idiag,si1,w &
            ,dt(l)%iext)
        DEALLOCATE(w)
    END IF
    ALLOCATE(ldd(ndd),lcg(2*(n-ndd)))
    IF (real(n).GT.targetcoarsefac*(n-ndd)) THEN
       skipass=.TRUE.
       npass1=npass+1
    ELSE
       skipass=.FALSE.
       npass1=npass
    END IF
    IF (l.GT.1) THEN
       IF (spd) THEN
          CALL sagmgpar_findpairs_SI(l,n,a,ja,ia,dt(l)%idiag,si1  &
               ,ind2,lcg,nc,ndd,ldd,skipass,iperm &
               ,dt(l)%iext)
       ELSE
          CALL sagmgpar_findpairs_GI(l,n,a,ja,ia,dt(l)%idiag,si1  &
               ,ind2,lcg,nc,ndd,ldd,skipass,iperm &
               ,dt(l)%iext)
       END IF
       DEALLOCATE(iperm)
    ELSE
       IF (spd) THEN
          CALL sagmgpar_findpairs_SI1(l,n,a,ja,ia,dt(l)%idiag,si1  &
               ,ind2,lcg,nc,ndd,ldd,skipass,riperm,iperm &
               ,dt(l)%iext)
       ELSE
          CALL sagmgpar_findpairs_GI1(l,n,a,ja,ia,dt(l)%idiag,si1  &
               ,ind2,lcg,nc,ndd,ldd,skipass,riperm,iperm &
               ,dt(l)%iext)
       END IF
       DEALLOCATE(iperm,riperm)
    END IF
10  CONTINUE
    nz=ia(n+1)-ia(1)
    IF (npass1.GT.1) THEN
       isize=nc
    ELSE
       isize=1
    END IF
    ALLOCATE(an(nz-2*(n-nc)+ndd),jan(nz-2*(n-nc)+ndd)             &
         ,ian(nc+1),idiagn(nc+1),iextn(nc+1),sinn(isize),wc(nc),iw(2*nc))
    CALL sagmgpar_setcg(n,a,ja,ia,dt(l)%idiag,si1,ind2,lcg,nc,an,jan,ian      &
         ,idiagn,sinn,npass1.GT.1,maxdg,iw,wc &
               ,dt(l)%iext,iextn)
    DEALLOCATE(wc,iw)
    nzc=ian(nc+1)-1
    IF (real(nz).GT.targetcoarsefac*nzc .OR. npass1.LE.1) THEN
       DEALLOCATE(si1,sinn,lcg)
       IF(ALLOCATED(ldd)) DEALLOCATE(ldd)
       dt(l)%ind  => ind2
       dt(l+1)%a    => an
       dt(l+1)%ja   => jan
       dt(l+1)%ia   => ian
       dt(l+1)%idiag=> idiagn
       NULLIFY(ind2,an,jan,ian,idiagn)
       dt(l+1)%iext=> iextn
       NULLIFY(iextn)
       GOTO 999
    END IF
    DEALLOCATE(ind2)
    lcg1 => lcg
    NULLIFY(lcg)
    m1=1
    DO kpass=2,npass1
       m1=2*m1
       np  = nc
       nzp = nzc
       ap     => an
       jap    => jan
       iap    => ian
       idiagp => idiagn
       sip    => sinn
       NULLIFY(an,jan,ian,idiagn,sinn)
       iextp  => iextn
       NULLIFY(iextn)
       ALLOCATE(lcg(2*np),ind2(np),w(maxdg),iw(maxdg))
       ind2(1:np)=-1
       IF (spd) THEN
          CALL sagmgpar_findpairs_SF(l,np,ap,jap,iap,idiagp,sip  &
               ,ind2,lcg,nc                                  &
               ,m1,lcg1,a,ja,ia,dt(l)%idiag,si1,w,iw &
               ,iextp,dt(l)%iext)
       ELSE
          CALL sagmgpar_findpairs_GF(l,np,ap,jap,iap,idiagp,sip     &
               ,ind2,lcg,nc                                     &
               ,m1,lcg1,a,ja,ia,dt(l)%idiag,si1,w,iw &
               ,iextp,dt(l)%iext)
       END IF
       DEALLOCATE(w,iw)
       IF (kpass.NE.npass1) THEN
          isize=nc
       ELSE
          isize=1
          DEALLOCATE(si1)
       END IF
       ALLOCATE(an(nzp-2*(np-nc)),jan(nzp-2*(np-nc))                 &
            ,ian(nc+1),idiagn(nc+1),iextn(nc+1),sinn(isize),wc(nc),iw(2*nc))
       CALL sagmgpar_setcg(np,ap,jap,iap,idiagp,sip,ind2,lcg,nc,an     &
            ,jan,ian,idiagn,sinn,kpass.NE.npass1,maxdg,iw,wc &
               ,iextp,iextn)
       DEALLOCATE(ap,jap,iap,idiagp,iextp,sip,ind2,wc,iw)
       ALLOCATE(lcgn(2*m1*nc))
       CALL sagmgpar_lcgmix(nc,m1,lcg1,lcg,lcgn)
       DEALLOCATE(lcg,lcg1)
       lcg1 => lcgn
       NULLIFY(lcgn)
       nzc=ian(nc+1)-1
       IF ( kpass.NE.npass1 .AND. real(nz).GT.targetcoarsefac*nzc ) THEN
          DEALLOCATE(si1)
          EXIT
       END IF
    END DO
    DEALLOCATE(sinn)
    ALLOCATE(dt(l)%ind(n))
    CALL sagmgpar_setind(nc,ndd,ldd,lcg1,2*m1,dt(l)%ind)
    DEALLOCATE(lcg1,ldd)
       dt(l+1)%a    => an
       dt(l+1)%ja   => jan
       dt(l+1)%ia   => ian
       dt(l+1)%idiag=> idiagn
       NULLIFY(an,jan,ian,idiagn)
       dt(l+1)%iext => iextn
       NULLIFY(iextn)
999 CONTINUE
    isize=MAX(nc,nzc)
    isize=MAX(isize,1)
    ALLOCATE (jan(isize))
    CALL sagmgpar_setextnum(nc,l,jan,nzp)
    IF (nc.EQ.0) THEN
        DEALLOCATE(jan)
        RETURN
    ENDIF
    nzp=dt(l+1)%ilstin(nneigh+1)-1
    ALLOCATE (ian(nzp),an(nzc))
    CALL sagmgpar_setcgext(l,nc,dt(l+1)%a,dt(l+1)%ja,dt(l+1)%ia   &
            ,dt(l+1)%idiag,dt(l+1)%iext,n,an,jan,ian,nzp)
    DEALLOCATE(dt(l+1)%a,dt(l+1)%ja,ian)
    dt(l+1)%a    => an
    dt(l+1)%ja   => jan
    RETURN
901 FORMAT(i3,'*SETUP: Coarsening by multiple pairwise aggregation')
902 FORMAT('****       Quality threshold (',A6,'):',f6.2, &
         ' ;  Strong diag. dom. trs:',f5.2)
903 FORMAT('****         Maximal number of passes:',i3,     &
         '  ; Target coarsening factor:',f5.2)
904 FORMAT('****           Diag. dom. checked w.r.t. sum of offdiag', &
          ' (no absolute vaues)')
905 FORMAT('****',22x,'Threshold for rows with large pos. offdiag.:',f5.2)
906 FORMAT('****  Rmk: Setup performed assuming the matrix symmetric')
908 FORMAT('****  Rmk: Setup performed for the transpose of the input matrix')
  END SUBROUTINE sagmgpar_aggregation
  SUBROUTINE sagmgpar_setCMK(n,ja,ia,idiag,riperm,iperm,iext)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),riperm(*),iperm(n)
    INTEGER, OPTIONAL :: iext(*)
    LOGICAL :: exc
    INTEGER :: i,j,jj,jk,jd,k,kk,j1,j2,i1,i2,ijs,ijs1,ijs2,dg,mindg,kdim
    INTEGER :: ifirst,ilast,kb,i0,i2rec,i2rec0
    ifirst=1
    ilast=n
    i2=ifirst
    mindg=n+1
    DO i = ifirst,ilast
       dg=ia(i+1)-ia(i)
       IF (dg.GT.1) THEN
          iperm(i)=-dg
          IF (dg.LT.mindg) THEN
             mindg=dg
             jj=i
          END IF
       ELSE
          riperm(i2)=i
          iperm(i)=i2
          i2=i2+1
       END IF
    ENDDO
    ijs=ifirst-1
    i1=i2
15  CONTINUE
    IF (i2.LE.ilast) THEN
      riperm(i2)=jj
      iperm(jj)=i2
    END IF
    DO WHILE (i1.LE.i2 .AND. i2.LT.ilast)
       i=riperm(i1)
       ijs1=i2+1
       j1 =ia(i)
       jd=idiag(i)
       j2 = iext(i)-1
       DO kk = j1,jd-1
          j=ja(kk)
          IF (iperm(j).LT.0) THEN
             i2=i2+1
             riperm(i2)=j
          END IF
       ENDDO
       DO kk = jd+1,j2
          j=ja(kk)
          IF (iperm(j).LT.0) THEN
             i2=i2+1
             riperm(i2)=j
          END IF
       ENDDO
       ijs2=i2
       exc=.TRUE..AND. ijs2.GT.ijs1
       DO WHILE(exc)
          exc=.FALSE.
          DO kk=ijs1+1,ijs2
             IF( iperm(riperm(kk)).GT.iperm(riperm(kk-1)) )THEN
                j=riperm(kk)
                riperm(kk)=riperm(kk-1)
                riperm(kk-1)=j
                exc=.TRUE.
             END IF
          END DO
       END DO
       DO kk=ijs1,ijs2
          iperm(riperm(kk))=kk
       END DO
       i1=i1+1
    END DO
    IF (i2.LT.ilast) THEN
       jj=0
       DO WHILE (jj.EQ.0)
          ijs=ijs+1
          IF (ijs.GT.ilast) THEN
             mindg=mindg+1
             ijs=ifirst
          END IF
          ijs1=ijs
          IF (iperm(ijs1).LT.0 .AND. ia(ijs1+1)-ia(ijs1).EQ.mindg) &
               jj=ijs1
       END DO
       i2=i2+1
       GOTO 15
    END IF
    RETURN
  END SUBROUTINE sagmgpar_setCMK
  SUBROUTINE sagmgpar_prepareagg(n,a,ja,ia,idiag,ind2,iperm,si,ndd,l,iext)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),ind2(n),iperm(n)
    INTEGER, OPTIONAL :: iext(*)
    INTEGER :: ndd, l
    REAL(kind(0.0e0)) :: a(*)
    REAL(kind(0.0e0)) , TARGET :: si(n)
    REAL(kind(0.0e0)) :: checkddl,oda,odm,ods,vald
    INTEGER :: i,j,jj,jk,jd,k,kk,j1,j2,i1,i2,ijs,ijs1,ijs2,dg,kdim,nnegrcs
    INTEGER :: ifirst,ilast,i0,kb
    REAL(kind(0.0e0)), POINTER, DIMENSION(:) :: odmax,odabs,osi
    REAL(kind(0.0e0)), ALLOCATABLE, DIMENSION(:,:), TARGET :: odd, buff2
    INTEGER :: ntot,n2,ii,kdeb,kfin,idest,ier
    IF (.NOT.spd) THEN
       checkddl=checkddJ
    ELSE
       checkddl=checkddB
    END IF
    IF (spd .AND. n.EQ.0) THEN
         ndd=0
         RETURN
    END IF
    IF (.NOT.spd) THEN
       ntot=n+dt(l)%ilstin(nneigh+1)-1
       n2=dt(l)%ilstout(nneigh+1)-1
       kdim=2
       IF (checkdd.GT.0) kdim=3
       ALLOCATE(odd(kdim,ntot),buff2(kdim,n2))
       osi   => odd(1,1:ntot)
       odmax => odd(2,1:ntot)
       IF (checkdd.GT.0) odabs => odd(3,1:ntot)
       odd=0.0e0
       ifirst=1
       ilast=n
        DO i=ilast,ifirst,-1
          j =ia(i)
          jd=idiag(i)
          jj=ia(i+1)-1
          DO k = j,jd-1
            jk=ja(k)
             osi(jk)=osi(jk)+a(k)
             odmax(jk)=max(odmax(jk),a(k))
             IF (checkdd.GT.0) odabs(jk)=odabs(jk)+abs(a(k))
          ENDDO
          DO k = jd+1,jj
            jk=ja(k)
             osi(jk)=osi(jk)+a(k)
             odmax(jk)=max(odmax(jk),a(k))
             IF (checkdd.GT.0) odabs(jk)=odabs(jk)+abs(a(k))
          ENDDO
        ENDDO
       IF (nneigho.GT.0) CALL MPI_WAITALL(nneigho,ireqo,MPI_STATUSES_IGNORE,ier)
       ii=0
       DO i=1,nneigh
          kdeb=dt(l)%ilstin(i)
          kfin=dt(l)%ilstin(i+1)-1
          IF (kfin.GE.kdeb) THEN
             ii=ii+1
             idest=ineigh(i)
             CALL MPI_ISEND(odd(1,n+kdeb),kdim*(kfin-kdeb+1)     &
                  ,MPI_REAL,idest,986,ICOMM,ireqi(ii),ier)
          END IF
       END DO
       nneigho=ii
       ii=0
       DO i=1,nneigh
          kdeb=dt(l)%ilstout(i)
          kfin=dt(l)%ilstout(i+1)-1
          IF (kfin.GE.kdeb) THEN
             ii=ii+1
             idest=ineigh(i)
             CALL MPI_IRECV(buff2(1,kdeb),kdim*(kfin-kdeb+1)    &
                  ,MPI_REAL,idest,986,icomm,ireqo(ii),ier)
          END IF
       END DO
       nneighi=ii
       IF (nneighi.GT.0) CALL MPI_WAITALL(nneighi,ireqo,MPI_STATUSES_IGNORE,ier)
      DO k=1,dt(l)%ilstout(nneigh+1)-1
         jk=dt(l)%lstout(k)
         osi(jk)=osi(jk)+buff2(1,k)
         odmax(jk)=max(odmax(jk),buff2(2,k))
         IF (checkdd.GT.0) odabs(jk)=odabs(jk)+buff2(3,k)
      END DO
    END IF
    ndd=0
    nnegrcs=0
    ifirst=1
    ilast=n
     DO i=ifirst,ilast
       j1 =ia(i)
       jd=idiag(i)
       j2 = ia (i+1)-1
       vald = a(jd)
       odm=0.0e0
       oda=0.0e0
       ods=0.0e0
       DO kk = j1,jd-1
          ods=ods+a(kk)
          odm=max(odm,a(kk))
          IF (checkdd.GT.0) oda=oda+abs(a(kk))
       ENDDO
       DO kk = jd+1,j2
          ods=ods+a(kk)
          odm=max(odm,a(kk))
          IF (checkdd.GT.0) oda=oda+abs(a(kk))
       ENDDO
       IF (.NOT.spd) THEN
          ods=(osi(i)+ods)/2
          odm=max(odm,odmax(i))
          IF (checkdd.GT.0) oda=(oda+odabs(i))/2
       END IF
       IF ((vald+ods).LT.-repsmach*ABS(vald)) nnegrcs=nnegrcs+1
       si(i)=-ods
       IF ( (checkdd.GT.0 .AND. vald.GT.checkddl*oda)     &
            .OR. (checkdd.LT.0 .AND. vald.GT.checkddl*abs(ods)) ) THEN
          ind2(i)=0
          ndd=ndd+1
       ELSE
          ind2(i)=-1
          IF (odm.GT.trspos*vald) iperm(i)=0
       ENDIF
     END DO
    zerors=.FALSE.
    IF (nnegrcs.GT.fracnegrcsum*n) THEN
       zerors=.TRUE.
       ndd=0
       ind2(1:n)=-1
    END IF
100 CONTINUE
    IF (.NOT.spd) THEN
       IF (nneigho.GT.0) CALL MPI_WAITALL(nneigho,ireqi,MPI_STATUSES_IGNORE,ier)
       nneigho=0
       DEALLOCATE(odd,buff2)
    END IF
    RETURN
  END SUBROUTINE sagmgpar_prepareagg
  SUBROUTINE sagmgpar_findpairs_GF(l,n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,m1,lcg1,a1,ja1,ia1,idiag1,si1,rtent,jtent &
    ,iext,iext1)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: l,n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    INTEGER :: iext(*)
    REAL(kind(0.0e0)) :: a(*)
    REAL(kind(0.0e0)) :: si(n)
    INTEGER :: m1,ja1(*),ia1(*),idiag1(*),jtent(*),lcg1(m1,n)
    INTEGER :: iext1(*)
    REAL(kind(0.0e0)) :: a1(*)
    REAL(kind(0.0e0)) :: si1(*),rtent(*)
    REAL(kind(0.0e0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0e0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2,nt,kb,i0
    LOGICAL :: acc,ifirst
    REAL(kind(0.0e0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
    ifirst=.TRUE.
  1 CONTINUE
    kaptg=imult*kaptg_dampJac
    bndmum1m1=1.0e0/(kaptg-1.0e0)
    dbndmum1=2*1.0e0/kaptg
    umdbndmum1=1.0e0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
      nt=n
      DO WHILE (nmark.LT.nt)
       isel=ijs
       ijs=ijs+1
       IF (ind(isel).GE.0) CYCLE
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       IF (lcg1(2,isel).EQ.0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       ntentleft=0
       i2=iext(isel)-1
       DO i = ia(isel),i2
          IF (i.EQ.idiag(isel)) CYCLE
          j = ja (i)
          IF(ind(j).GE.0) CYCLE
          IF(lcg1(2,j).EQ.0) CYCLE
          kk=0
          IF (i.LT.idiag(isel)) THEN
             j2=iext(j)-1
             DO jk=idiag(j)+1,j2
                IF (ja(jk).EQ.isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ELSE
             DO jk=ia(j),idiag(j)-1
                IF (ja(jk).EQ.isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ENDIF
          vals=-a(i)/2
          IF(kk.NE.0) vals=vals-a(kk)/2
          IF (zerors) THEN
             rsi=0.0e0
             rsj=0.0e0
             eta1=2*si(isel)
             eta2=2*si(j)
          ELSE
             rsi=-si(isel)+a(idiag(isel))
             rsj=-si(j)+a(idiag(j))
             eta1=2*a(idiag(isel))
             eta2=2*a(idiag(j))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          IF (sig1.GT.0.0e0) THEN
             del1=rsi
          ELSE
             del1=rsi+2*sig1
          END IF
          IF (sig2.GT.0.0e0) THEN
             del2=rsj
          ELSE
             del2=rsj+2*sig2
          END IF
          IF (vals.GT.0.0e0) THEN
             epsr=repsmach*vals
             IF (ABS(del1).LT.epsr .AND. ABS(del2).LT.epsr) THEN
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1).LT.epsr) THEN
                IF (del2.LT.-epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2).LT.epsr) THEN
                IF (del1.LT.-epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12.LT.-epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp.LT.0.0e0) CYCLE
                valp=((eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1.LE.0.0e0 .OR. del2.LE.0.0e0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp.LT.0.0e0) CYCLE
             vals=(eta1*eta2)/(eta1+eta2)
             valp=vals/valp
          END IF
          IF (valp.GT.kaptg) CYCLE
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          ntentleft=ntentleft+1
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. j.LT.ipair)  GOTO 9
          rtent(ntentleft)=tent
          jtent(ntentleft)=j
          CYCLE
9         CONTINUE
          rtent(ntentleft)=val
          jtent(ntentleft)=ipair
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       IF (ipair.EQ.0) GOTO 25
20     CONTINUE
       CALL sagmgpar_checktentagg_GF
       IF (.NOT.acc) THEN
          ipair = 0
          IF (ntentleft.GT.0) THEN
             i=1
             j=1
             DO WHILE (i.LE.ntentleft)
                IF (jtent(j).GT.0) THEN
                   tent=rtent(j)
                   IF (ipair.EQ.0) GOTO 22
                   IF (16*(tent-val).LT.-1) GOTO 22
                   IF (16*(tent-val).LT.1 .AND. j.LT.ipair) GOTO 22
                   GOTO 23
22                 CONTINUE
                   val=tent
                   ipair=jtent(j)
                   ijtent=j
23                 CONTINUE
                   i=i+1
                END IF
                j=j+1
             END DO
             ntentleft=ntentleft-1
             jtent(ijtent)=0
             GOTO 20
          END IF
       END IF
25     CONTINUE
       IF (ipair.EQ.0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
      ENDDO
    RETURN
  CONTAINS
  SUBROUTINE sagmgpar_checktentagg_GF
    INTEGER , PARAMETER :: mm=max(2**(npass+1),8)
    REAL(kind(0.0e0)) :: W(mm,mm), sig(mm), AGe(mm), v(mm)
    REAL(kind(0.0e0)) :: alpha, alp, tmp, beta, beta0, f1, f2
    INTEGER :: j,jj,k,l,m, setdim, l2, k2
    INTEGER  :: info, setdim1
    INTEGER :: set(mm), l1, wdthT
    REAL(kind(0.0e0)) :: T
    LOGICAL :: exc
    IF (m1.eq.2) THEN
       IF (lcg1(2,isel).LT.0) THEN
          IF (lcg1(2,ipair).LT.0) THEN
             set(1)=lcg1(1,isel)
             set(2)=lcg1(1,ipair)
             setdim=2
          ELSE
             set(1)=lcg1(1,isel)
             set(2)=lcg1(1,ipair)
             set(3)=lcg1(2,ipair)
             setdim=3
          END IF
          l1=1
       ELSE
          IF (lcg1(2,ipair).LT.0) THEN
             set(1)=lcg1(1,isel)
             set(2)=lcg1(2,isel)
             set(3)=lcg1(1,ipair)
             setdim=3
          ELSE
             set(1)=lcg1(1,isel)
             set(2)=lcg1(2,isel)
             set(3)=lcg1(1,ipair)
             set(4)=lcg1(2,ipair)
             setdim=4
          END IF
          l1=2
       END IF
    ELSE
       l1=m1
       IF (lcg1(m1,isel).LT.0) l1=-lcg1(m1,isel)
       set(1:l1)=lcg1(1:l1,isel)
       l2=m1
       IF (lcg1(m1,ipair).LT.0) l2=-lcg1(m1,ipair)
       set(l1+1:l1+l2)=lcg1(1:l2,ipair)
       setdim=l1+l2
    END IF
    exc=.TRUE.
    DO WHILE(exc)
       exc=.FALSE.
       DO l=2,SetDim
          IF( set(l).LT.set(l-1) )THEN
             jj=set(l)
             set(l)=set(l-1)
             set(l-1)=jj
             exc=.TRUE.
          END IF
       END DO
    END DO
    DO j=1,SetDim
       jj=Set(j)
       sig(j)=si1(jj)
       IF (zerors) THEN
          W(j,j)=sig(j)
          AGe(j)=0.0e0
       ELSE
          W(j,j)=a1(idiag1(jj))
          AGe(j)=W(j,j)-sig(j)
       END IF
       l2=j+1
       DO l=l2,SetDim
          W(j,l)=0.0e0
          W(l,j)=0.0e0
       END DO
       k2=iext1(jj)-1
       DO k=idiag1(jj)+1,k2
          DO l=l2,SetDim
             m=Set(l)
             IF(ja1(k).EQ.m)THEN
                alpha=a1(k)/2
                W(j,l)=alpha
                W(l,j)=alpha
                EXIT
             END IF
          END DO
       END DO
       DO k=ia1(jj),idiag1(jj)-1
          DO l=1,j-1
             m=Set(l)
             IF(ja1(k).EQ.m)THEN
                alpha=a1(k)/2
                W(j,l)=W(j,l)+alpha
                W(l,j)=W(j,l)
                EXIT
             END IF
          END DO
       END DO
    END DO
    DO j=1,SetDim
       DO k=1,SetDim
          IF (j.ne.k) THEN
             sig(j)=sig(j)+W(j,k)
          END IF
       ENDDO
       IF (sig(j).LT.0.0e0)  AGe(j)=AGe(j)+2*sig(j)
       v(j)=W(j,j)
       W(j,j)=umdbndmum1*W(j,j)-abs(sig(j))
       IF (j.EQ.1) THEN
          beta0=v(j)
          alp=abs(AGe(j))
       ELSE
          beta0=beta0+v(j)
          alp=max(alp,abs(AGe(j)))
       END IF
    END DO
    beta=dbndmum1/beta0
    DO j=1,SetDim
       DO k=1,SetDim
          W(j,k)=W(j,k)+beta*v(j)*v(k)
       END DO
    END DO
    IF (alp.LT.repsmach*beta0) THEN
       SetDim1=SetDim-1
    ELSE
       SetDim1=SetDim
    END IF
    acc=.FALSE.
    SELECT CASE (SetDim1)
    CASE (1)
       GOTO 11
    CASE (2)
       GOTO 12
    CASE (3)
       GOTO 13
    CASE (4)
       GOTO 14
    CASE (5)
       GOTO 15
    CASE (6)
       GOTO 16
    CASE (7)
       GOTO 17
    CASE (8)
       GOTO 18
    CASE DEFAULT
       CALL SPOTRF('U',SetDim1,W,mm,info)
       IF (info.NE.0) RETURN
       GOTO 10
    END SELECT
18  CONTINUE
    IF (W(8,8).LE.0.0e0) RETURN
    W(7,7) = W(7,7) - (W(7,8)/W(8,8)) * W(7,8)
    T = W(6,8)/W(8,8)
    W(6,7) = W(6,7) - T * W(7,8)
    W(6,6) = W(6,6) - T * W(6,8)
    T = W(5,8)/W(8,8)
    W(5,7) = W(5,7) - T * W(7,8)
    W(5,6) = W(5,6) - T * W(6,8)
    W(5,5) = W(5,5) - T * W(5,8)
    T = W(4,8)/W(8,8)
    W(4,7) = W(4,7) - T * W(7,8)
    W(4,6) = W(4,6) - T * W(6,8)
    W(4,5) = W(4,5) - T * W(5,8)
    W(4,4) = W(4,4) - T * W(4,8)
    T = W(3,8)/W(8,8)
    W(3,7) = W(3,7) - T * W(7,8)
    W(3,6) = W(3,6) - T * W(6,8)
    W(3,5) = W(3,5) - T * W(5,8)
    W(3,4) = W(3,4) - T * W(4,8)
    W(3,3) = W(3,3) - T * W(3,8)
    T = W(2,8)/W(8,8)
    W(2,7) = W(2,7) - T * W(7,8)
    W(2,6) = W(2,6) - T * W(6,8)
    W(2,5) = W(2,5) - T * W(5,8)
    W(2,4) = W(2,4) - T * W(4,8)
    W(2,3) = W(2,3) - T * W(3,8)
    W(2,2) = W(2,2) - T * W(2,8)
    T = W(1,8)/W(8,8)
    W(1,7) = W(1,7) - T * W(7,8)
    W(1,6) = W(1,6) - T * W(6,8)
    W(1,5) = W(1,5) - T * W(5,8)
    W(1,4) = W(1,4) - T * W(4,8)
    W(1,3) = W(1,3) - T * W(3,8)
    W(1,2) = W(1,2) - T * W(2,8)
    W(1,1) = W(1,1) - T * W(1,8)
17  CONTINUE
    IF (W(7,7).LE.0.0e0) RETURN
    W(6,6) = W(6,6) - (W(6,7)/W(7,7)) * W(6,7)
    T = W(5,7)/W(7,7)
    W(5,6) = W(5,6) - T * W(6,7)
    W(5,5) = W(5,5) - T * W(5,7)
    T = W(4,7)/W(7,7)
    W(4,6) = W(4,6) - T * W(6,7)
    W(4,5) = W(4,5) - T * W(5,7)
    W(4,4) = W(4,4) - T * W(4,7)
    T = W(3,7)/W(7,7)
    W(3,6) = W(3,6) - T * W(6,7)
    W(3,5) = W(3,5) - T * W(5,7)
    W(3,4) = W(3,4) - T * W(4,7)
    W(3,3) = W(3,3) - T * W(3,7)
    T = W(2,7)/W(7,7)
    W(2,6) = W(2,6) - T * W(6,7)
    W(2,5) = W(2,5) - T * W(5,7)
    W(2,4) = W(2,4) - T * W(4,7)
    W(2,3) = W(2,3) - T * W(3,7)
    W(2,2) = W(2,2) - T * W(2,7)
    T = W(1,7)/W(7,7)
    W(1,6) = W(1,6) - T * W(6,7)
    W(1,5) = W(1,5) - T * W(5,7)
    W(1,4) = W(1,4) - T * W(4,7)
    W(1,3) = W(1,3) - T * W(3,7)
    W(1,2) = W(1,2) - T * W(2,7)
    W(1,1) = W(1,1) - T * W(1,7)
16  CONTINUE
    IF (W(6,6).LE.0.0e0) RETURN
    W(5,5) = W(5,5) - (W(5,6)/W(6,6)) * W(5,6)
    T = W(4,6)/W(6,6)
    W(4,5) = W(4,5) - T * W(5,6)
    W(4,4) = W(4,4) - T * W(4,6)
    T = W(3,6)/W(6,6)
    W(3,5) = W(3,5) - T * W(5,6)
    W(3,4) = W(3,4) - T * W(4,6)
    W(3,3) = W(3,3) - T * W(3,6)
    T = W(2,6)/W(6,6)
    W(2,5) = W(2,5) - T * W(5,6)
    W(2,4) = W(2,4) - T * W(4,6)
    W(2,3) = W(2,3) - T * W(3,6)
    W(2,2) = W(2,2) - T * W(2,6)
    T = W(1,6)/W(6,6)
    W(1,5) = W(1,5) - T * W(5,6)
    W(1,4) = W(1,4) - T * W(4,6)
    W(1,3) = W(1,3) - T * W(3,6)
    W(1,2) = W(1,2) - T * W(2,6)
    W(1,1) = W(1,1) - T * W(1,6)
15  CONTINUE
    IF (W(5,5).LE.0.0e0) RETURN
    W(4,4) = W(4,4) - (W(4,5)/W(5,5)) * W(4,5)
    T = W(3,5)/W(5,5)
    W(3,4) = W(3,4) - T * W(4,5)
    W(3,3) = W(3,3) - T * W(3,5)
    T = W(2,5)/W(5,5)
    W(2,4) = W(2,4) - T * W(4,5)
    W(2,3) = W(2,3) - T * W(3,5)
    W(2,2) = W(2,2) - T * W(2,5)
    T = W(1,5)/W(5,5)
    W(1,4) = W(1,4) - T * W(4,5)
    W(1,3) = W(1,3) - T * W(3,5)
    W(1,2) = W(1,2) - T * W(2,5)
    W(1,1) = W(1,1) - T * W(1,5)
14  CONTINUE
    IF (W(4,4).LE.0.0e0) RETURN
    W(3,3) = W(3,3) - (W(3,4)/W(4,4)) * W(3,4)
    T = W(2,4)/W(4,4)
    W(2,3) = W(2,3) - T * W(3,4)
    W(2,2) = W(2,2) - T * W(2,4)
    T = W(1,4)/W(4,4)
    W(1,3) = W(1,3) - T * W(3,4)
    W(1,2) = W(1,2) - T * W(2,4)
    W(1,1) = W(1,1) - T * W(1,4)
13  CONTINUE
    IF (W(3,3).LE.0.0e0) RETURN
    W(2,2) = W(2,2) - (W(2,3)/W(3,3)) * W(2,3)
    T = W(1,3)/W(3,3)
    W(1,2) = W(1,2) - T * W(2,3)
    W(1,1) = W(1,1) - T * W(1,3)
12  CONTINUE
    IF (W(2,2).LE.0.0e0) RETURN
    W(1,1) = W(1,1) - (W(1,2)/W(2,2)) * W(1,2)
11  CONTINUE
    IF (W(1,1).LE.0.0e0) RETURN
10  CONTINUE
    acc=.TRUE.
    RETURN
  END SUBROUTINE sagmgpar_checktentagg_GF
  END SUBROUTINE sagmgpar_findpairs_GF
  SUBROUTINE sagmgpar_findpairs_SF(l,n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,m1,lcg1,a1,ja1,ia1,idiag1,si1,rtent,jtent &
    ,iext,iext1)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: l,n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    INTEGER :: iext(*)
    REAL(kind(0.0e0)) :: a(*)
    REAL(kind(0.0e0)) :: si(n)
    INTEGER :: m1,ja1(*),ia1(*),idiag1(*),jtent(*),lcg1(m1,n)
    INTEGER :: iext1(*)
    REAL(kind(0.0e0)) :: a1(*)
    REAL(kind(0.0e0)) :: si1(*),rtent(*)
    REAL(kind(0.0e0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0e0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2,nt,kb,i0
    LOGICAL :: acc,ifirst
    REAL(kind(0.0e0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
    ifirst=.TRUE.
  1 CONTINUE
    kaptg=imult*kaptg_blocdia
    bndmum1m1=1.0e0/(kaptg-1.0e0)
    dbndmum1=2*1.0e0/kaptg
    umdbndmum1=1.0e0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
      nt=n
      DO WHILE (nmark.LT.nt)
       isel=ijs
       ijs=ijs+1
       IF (ind(isel).GE.0) CYCLE
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       IF (lcg1(2,isel).EQ.0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       ntentleft=0
       i2=iext(isel)-1
       DO i = ia(isel),i2
          IF (i.EQ.idiag(isel)) CYCLE
          j = ja (i)
          IF(ind(j).GE.0) CYCLE
          IF(lcg1(2,j).EQ.0) CYCLE
          vals=-a(i)
          IF (zerors) THEN
             rsi=0.0e0
             rsj=0.0e0
          ELSE
             rsi=-si(isel)+a(idiag(isel))
             rsj=-si(j)+a(idiag(j))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          epsr=repsmach*ABS(vals)
          IF (sig1.GT.0.0e0) THEN
             del1=rsi
             eta1=rsi+2*sig1
          ELSE
             del1=rsi+2*sig1
             eta1=rsi
          END IF
          IF (eta1.LT.-epsr) CYCLE
          IF (sig2.GT.0.0e0) THEN
             del2=rsj
             eta2=rsj+2*sig2
          ELSE
             del2=rsj+2*sig2
             eta2=rsj
          END IF
          IF (eta2.LT.-epsr) CYCLE
          IF (vals.GT.0.0e0) THEN
             IF (ABS(del1).LT.epsr .AND. ABS(del2).LT.epsr) THEN
               valp=1.0e0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1).LT.epsr) THEN
                IF (del2.LT.-epsr) CYCLE
                valp=1.0e0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2).LT.epsr) THEN
                IF (del1.LT.-epsr) CYCLE
                valp=1.0e0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12.LT.-epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp.LT.0.0e0) CYCLE
                valp=(vals+(eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1.LE.0.0e0 .OR. del2.LE.0.0e0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp.LT.0.0e0) CYCLE
             vals=vals+(eta1*eta2)/(eta1+eta2)
             IF (vals.LT.0.0e0) CYCLE
             valp=vals/valp
          END IF
          IF (valp.GT.kaptg) CYCLE
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          ntentleft=ntentleft+1
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. j.LT.ipair)  GOTO 9
          rtent(ntentleft)=tent
          jtent(ntentleft)=j
          CYCLE
9         CONTINUE
          rtent(ntentleft)=val
          jtent(ntentleft)=ipair
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       IF (ipair.EQ.0) GOTO 25
20     CONTINUE
       CALL sagmgpar_checktentagg_SF
       IF (.NOT.acc) THEN
          ipair = 0
          IF (ntentleft.GT.0) THEN
             i=1
             j=1
             DO WHILE (i.LE.ntentleft)
                IF (jtent(j).GT.0) THEN
                   tent=rtent(j)
                   IF (ipair.EQ.0) GOTO 22
                   IF (16*(tent-val).LT.-1) GOTO 22
                   IF (16*(tent-val).LT.1 .AND. j.LT.ipair) GOTO 22
                   GOTO 23
22                 CONTINUE
                   val=tent
                   ipair=jtent(j)
                   ijtent=j
23                 CONTINUE
                   i=i+1
                END IF
                j=j+1
             END DO
             ntentleft=ntentleft-1
             jtent(ijtent)=0
             GOTO 20
          END IF
       END IF
25     CONTINUE
       IF (ipair.EQ.0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
      ENDDO
    RETURN
  CONTAINS
  SUBROUTINE sagmgpar_checktentagg_SF
    INTEGER , PARAMETER :: mm=max(2**(npass+1),8)
    REAL(kind(0.0e0)) :: W(mm,mm), sig(mm), AGe(mm), v(mm)
    REAL(kind(0.0e0)) :: alpha, alp, tmp, beta, beta0, f1, f2
    INTEGER :: j,jj,k,l,m, setdim, l2, k2
    INTEGER  :: info, setdim1
    INTEGER :: set(mm), l1, wdthT
    REAL(kind(0.0e0)) :: T
    LOGICAL :: exc
    IF (m1.eq.2) THEN
       IF (lcg1(2,isel).LT.0) THEN
          IF (lcg1(2,ipair).LT.0) THEN
             set(1)=lcg1(1,isel)
             set(2)=lcg1(1,ipair)
             setdim=2
          ELSE
             set(1)=lcg1(1,isel)
             set(2)=lcg1(1,ipair)
             set(3)=lcg1(2,ipair)
             setdim=3
          END IF
          l1=1
       ELSE
          IF (lcg1(2,ipair).LT.0) THEN
             set(1)=lcg1(1,isel)
             set(2)=lcg1(2,isel)
             set(3)=lcg1(1,ipair)
             setdim=3
          ELSE
             set(1)=lcg1(1,isel)
             set(2)=lcg1(2,isel)
             set(3)=lcg1(1,ipair)
             set(4)=lcg1(2,ipair)
             setdim=4
          END IF
          l1=2
       END IF
    ELSE
       l1=m1
       IF (lcg1(m1,isel).LT.0) l1=-lcg1(m1,isel)
       set(1:l1)=lcg1(1:l1,isel)
       l2=m1
       IF (lcg1(m1,ipair).LT.0) l2=-lcg1(m1,ipair)
       set(l1+1:l1+l2)=lcg1(1:l2,ipair)
       setdim=l1+l2
    END IF
    exc=.TRUE.
    DO WHILE(exc)
       exc=.FALSE.
       DO l=2,SetDim
          IF( set(l).LT.set(l-1) )THEN
             jj=set(l)
             set(l)=set(l-1)
             set(l-1)=jj
             exc=.TRUE.
          END IF
       END DO
    END DO
    DO j=1,SetDim
       jj=Set(j)
       sig(j)=si1(jj)
       IF (zerors) THEN
          W(j,j)=sig(j)
          AGe(j)=0.0e0
       ELSE
          W(j,j)=a1(idiag1(jj))
          AGe(j)=W(j,j)-sig(j)
       END IF
       l2=j+1
       DO l=l2,SetDim
          W(j,l)=0.0e0
          W(l,j)=0.0e0
       END DO
       k2=iext1(jj)-1
       DO k=idiag1(jj)+1,k2
          DO l=l2,SetDim
             m=Set(l)
             IF(ja1(k).EQ.m)THEN
                alpha=a1(k)
                W(j,l)=alpha
                W(l,j)=alpha
                EXIT
             END IF
          END DO
       END DO
    END DO
    DO j=1,SetDim
       DO k=1,SetDim
          IF (j.ne.k) THEN
             sig(j)=sig(j)+W(j,k)
          END IF
       ENDDO
       IF (sig(j).LT.0.0e0)  AGe(j)=AGe(j)+2*sig(j)
       W(j,j)=W(j,j)-abs(sig(j))
       tmp=2*abs(sig(j))
       W(j,j)=W(j,j)-bndmum1m1*tmp
       v(j)=tmp+AGe(j)
       IF (j.EQ.1) THEN
          beta0=v(j)
          alp=abs(AGe(j))
       ELSE
          beta0=beta0+v(j)
          alp=max(alp,abs(AGe(j)))
       END IF
    END DO
    beta=bndmum1m1/beta0
    DO j=1,SetDim
       DO k=1,SetDim
          W(j,k)=W(j,k)+beta*v(j)*v(k)
       END DO
    END DO
    IF (alp.LT.repsmach*beta0) THEN
       SetDim1=SetDim-1
    ELSE
       SetDim1=SetDim
    END IF
    acc=.FALSE.
    SELECT CASE (SetDim1)
    CASE (1)
       GOTO 11
    CASE (2)
       GOTO 12
    CASE (3)
       GOTO 13
    CASE (4)
       GOTO 14
    CASE (5)
       GOTO 15
    CASE (6)
       GOTO 16
    CASE (7)
       GOTO 17
    CASE (8)
       GOTO 18
    CASE DEFAULT
       CALL SPOTRF('U',SetDim1,W,mm,info)
       IF (info.NE.0) RETURN
       GOTO 10
    END SELECT
18  CONTINUE
    IF (W(8,8).LE.0.0e0) RETURN
    W(7,7) = W(7,7) - (W(7,8)/W(8,8)) * W(7,8)
    T = W(6,8)/W(8,8)
    W(6,7) = W(6,7) - T * W(7,8)
    W(6,6) = W(6,6) - T * W(6,8)
    T = W(5,8)/W(8,8)
    W(5,7) = W(5,7) - T * W(7,8)
    W(5,6) = W(5,6) - T * W(6,8)
    W(5,5) = W(5,5) - T * W(5,8)
    T = W(4,8)/W(8,8)
    W(4,7) = W(4,7) - T * W(7,8)
    W(4,6) = W(4,6) - T * W(6,8)
    W(4,5) = W(4,5) - T * W(5,8)
    W(4,4) = W(4,4) - T * W(4,8)
    T = W(3,8)/W(8,8)
    W(3,7) = W(3,7) - T * W(7,8)
    W(3,6) = W(3,6) - T * W(6,8)
    W(3,5) = W(3,5) - T * W(5,8)
    W(3,4) = W(3,4) - T * W(4,8)
    W(3,3) = W(3,3) - T * W(3,8)
    T = W(2,8)/W(8,8)
    W(2,7) = W(2,7) - T * W(7,8)
    W(2,6) = W(2,6) - T * W(6,8)
    W(2,5) = W(2,5) - T * W(5,8)
    W(2,4) = W(2,4) - T * W(4,8)
    W(2,3) = W(2,3) - T * W(3,8)
    W(2,2) = W(2,2) - T * W(2,8)
    T = W(1,8)/W(8,8)
    W(1,7) = W(1,7) - T * W(7,8)
    W(1,6) = W(1,6) - T * W(6,8)
    W(1,5) = W(1,5) - T * W(5,8)
    W(1,4) = W(1,4) - T * W(4,8)
    W(1,3) = W(1,3) - T * W(3,8)
    W(1,2) = W(1,2) - T * W(2,8)
    W(1,1) = W(1,1) - T * W(1,8)
17  CONTINUE
    IF (W(7,7).LE.0.0e0) RETURN
    W(6,6) = W(6,6) - (W(6,7)/W(7,7)) * W(6,7)
    T = W(5,7)/W(7,7)
    W(5,6) = W(5,6) - T * W(6,7)
    W(5,5) = W(5,5) - T * W(5,7)
    T = W(4,7)/W(7,7)
    W(4,6) = W(4,6) - T * W(6,7)
    W(4,5) = W(4,5) - T * W(5,7)
    W(4,4) = W(4,4) - T * W(4,7)
    T = W(3,7)/W(7,7)
    W(3,6) = W(3,6) - T * W(6,7)
    W(3,5) = W(3,5) - T * W(5,7)
    W(3,4) = W(3,4) - T * W(4,7)
    W(3,3) = W(3,3) - T * W(3,7)
    T = W(2,7)/W(7,7)
    W(2,6) = W(2,6) - T * W(6,7)
    W(2,5) = W(2,5) - T * W(5,7)
    W(2,4) = W(2,4) - T * W(4,7)
    W(2,3) = W(2,3) - T * W(3,7)
    W(2,2) = W(2,2) - T * W(2,7)
    T = W(1,7)/W(7,7)
    W(1,6) = W(1,6) - T * W(6,7)
    W(1,5) = W(1,5) - T * W(5,7)
    W(1,4) = W(1,4) - T * W(4,7)
    W(1,3) = W(1,3) - T * W(3,7)
    W(1,2) = W(1,2) - T * W(2,7)
    W(1,1) = W(1,1) - T * W(1,7)
16  CONTINUE
    IF (W(6,6).LE.0.0e0) RETURN
    W(5,5) = W(5,5) - (W(5,6)/W(6,6)) * W(5,6)
    T = W(4,6)/W(6,6)
    W(4,5) = W(4,5) - T * W(5,6)
    W(4,4) = W(4,4) - T * W(4,6)
    T = W(3,6)/W(6,6)
    W(3,5) = W(3,5) - T * W(5,6)
    W(3,4) = W(3,4) - T * W(4,6)
    W(3,3) = W(3,3) - T * W(3,6)
    T = W(2,6)/W(6,6)
    W(2,5) = W(2,5) - T * W(5,6)
    W(2,4) = W(2,4) - T * W(4,6)
    W(2,3) = W(2,3) - T * W(3,6)
    W(2,2) = W(2,2) - T * W(2,6)
    T = W(1,6)/W(6,6)
    W(1,5) = W(1,5) - T * W(5,6)
    W(1,4) = W(1,4) - T * W(4,6)
    W(1,3) = W(1,3) - T * W(3,6)
    W(1,2) = W(1,2) - T * W(2,6)
    W(1,1) = W(1,1) - T * W(1,6)
15  CONTINUE
    IF (W(5,5).LE.0.0e0) RETURN
    W(4,4) = W(4,4) - (W(4,5)/W(5,5)) * W(4,5)
    T = W(3,5)/W(5,5)
    W(3,4) = W(3,4) - T * W(4,5)
    W(3,3) = W(3,3) - T * W(3,5)
    T = W(2,5)/W(5,5)
    W(2,4) = W(2,4) - T * W(4,5)
    W(2,3) = W(2,3) - T * W(3,5)
    W(2,2) = W(2,2) - T * W(2,5)
    T = W(1,5)/W(5,5)
    W(1,4) = W(1,4) - T * W(4,5)
    W(1,3) = W(1,3) - T * W(3,5)
    W(1,2) = W(1,2) - T * W(2,5)
    W(1,1) = W(1,1) - T * W(1,5)
14  CONTINUE
    IF (W(4,4).LE.0.0e0) RETURN
    W(3,3) = W(3,3) - (W(3,4)/W(4,4)) * W(3,4)
    T = W(2,4)/W(4,4)
    W(2,3) = W(2,3) - T * W(3,4)
    W(2,2) = W(2,2) - T * W(2,4)
    T = W(1,4)/W(4,4)
    W(1,3) = W(1,3) - T * W(3,4)
    W(1,2) = W(1,2) - T * W(2,4)
    W(1,1) = W(1,1) - T * W(1,4)
13  CONTINUE
    IF (W(3,3).LE.0.0e0) RETURN
    W(2,2) = W(2,2) - (W(2,3)/W(3,3)) * W(2,3)
    T = W(1,3)/W(3,3)
    W(1,2) = W(1,2) - T * W(2,3)
    W(1,1) = W(1,1) - T * W(1,3)
12  CONTINUE
    IF (W(2,2).LE.0.0e0) RETURN
    W(1,1) = W(1,1) - (W(1,2)/W(2,2)) * W(1,2)
11  CONTINUE
    IF (W(1,1).LE.0.0e0) RETURN
10  CONTINUE
    acc=.TRUE.
    RETURN
  END SUBROUTINE sagmgpar_checktentagg_SF
  END SUBROUTINE sagmgpar_findpairs_SF
  SUBROUTINE sagmgpar_findpairs_GI(l,n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,nddl,ldd,skipass                                      &
    ,ipc &
    ,iext)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: l,n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    INTEGER :: iext(*)
    REAL(kind(0.0e0)) :: a(*)
    REAL(kind(0.0e0)) :: si(n)
    INTEGER :: ldd(nddl),nddl
    LOGICAL :: skipass
    INTEGER :: ipc(n)
    REAL(kind(0.0e0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0e0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2,nt,kb,i0
    LOGICAL :: acc,ifirst
    REAL(kind(0.0e0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
    ifirst=.TRUE.
  1 CONTINUE
    kaptg=imult*kaptg_dampJac
    bndmum1m1=1.0e0/(kaptg-1.0e0)
    dbndmum1=2*1.0e0/kaptg
    umdbndmum1=1.0e0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
      nt=n
      DO WHILE (nmark.LT.nt)
       isel=ijs
       ijs=ijs+1
       IF (ind(isel).EQ.0) THEN
          idd=idd+1
          nmark=nmark+1
          ldd(idd)=isel
          CYCLE
       END IF
       IF (ind(isel).GE.0) CYCLE
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       IF (ipc(isel).EQ.0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       IF (skipass) THEN
          lcg(2,nc) = -1
          CYCLE
       END IF
       i2=iext(isel)-1
       DO i = ia(isel),i2
          IF (i.EQ.idiag(isel)) CYCLE
          j = ja (i)
          IF(ind(j).GE.0) CYCLE
          IF(ipc(j).EQ.0) CYCLE
          kk=0
          IF (i.LT.idiag(isel)) THEN
             j2=iext(j)-1
             DO jk=idiag(j)+1,j2
                IF (ja(jk).EQ.isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ELSE
             DO jk=ia(j),idiag(j)-1
                IF (ja(jk).EQ.isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ENDIF
          vals=-a(i)/2
          IF(kk.NE.0) vals=vals-a(kk)/2
          IF (zerors) THEN
             rsi=0.0e0
             rsj=0.0e0
             eta1=2*si(isel)
             eta2=2*si(j)
          ELSE
             rsi=-si(isel)+a(idiag(isel))
             rsj=-si(j)+a(idiag(j))
             eta1=2*a(idiag(isel))
             eta2=2*a(idiag(j))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          IF (sig1.GT.0.0e0) THEN
             del1=rsi
          ELSE
             del1=rsi+2*sig1
          END IF
          IF (sig2.GT.0.0e0) THEN
             del2=rsj
          ELSE
             del2=rsj+2*sig2
          END IF
          IF (vals.GT.0.0e0) THEN
             epsr=repsmach*vals
             IF (ABS(del1).LT.epsr .AND. ABS(del2).LT.epsr) THEN
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1).LT.epsr) THEN
                IF (del2.LT.-epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2).LT.epsr) THEN
                IF (del1.LT.-epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12.LT.-epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp.LT.0.0e0) CYCLE
                valp=((eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1.LE.0.0e0 .OR. del2.LE.0.0e0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp.LT.0.0e0) CYCLE
             vals=(eta1*eta2)/(eta1+eta2)
             valp=vals/valp
          END IF
          IF (valp.GT.kaptg) CYCLE
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. j.LT.ipair)  GOTO 9
          CYCLE
9         CONTINUE
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       IF (ipair.EQ.0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
      ENDDO
    IF (ifirst .AND. 3*nc.GT.2*n .AND. ngl(1).GT.maxcoarseslowt) THEN
       imult=2*imult
       ifirst=.FALSE.
       IF (wfo) THEN
          WRITE (iout,901) IRANK,imult*kaptg_dampJac
       END IF
901    FORMAT(i3, &
       '*   Coarsening too slow: quality threshold increased to',es9.2)
       ind(1:n)=-1
       DO i=1,nddl
          ind(ldd(i))=0
       END DO
       GOTO 1
    ENDIF
    RETURN
  END SUBROUTINE sagmgpar_findpairs_GI
  SUBROUTINE sagmgpar_findpairs_SI(l,n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,nddl,ldd,skipass                                      &
    ,ipc &
    ,iext)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: l,n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    INTEGER :: iext(*)
    REAL(kind(0.0e0)) :: a(*)
    REAL(kind(0.0e0)) :: si(n)
    INTEGER :: ldd(nddl),nddl
    LOGICAL :: skipass
    INTEGER :: ipc(n)
    REAL(kind(0.0e0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0e0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2,nt,kb,i0
    LOGICAL :: acc,ifirst
    REAL(kind(0.0e0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
    ifirst=.TRUE.
  1 CONTINUE
    kaptg=imult*kaptg_blocdia
    bndmum1m1=1.0e0/(kaptg-1.0e0)
    dbndmum1=2*1.0e0/kaptg
    umdbndmum1=1.0e0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
      nt=n
      DO WHILE (nmark.LT.nt)
       isel=ijs
       ijs=ijs+1
       IF (ind(isel).EQ.0) THEN
          idd=idd+1
          nmark=nmark+1
          ldd(idd)=isel
          CYCLE
       END IF
       IF (ind(isel).GE.0) CYCLE
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       IF (ipc(isel).EQ.0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       IF (skipass) THEN
          lcg(2,nc) = -1
          CYCLE
       END IF
       i2=iext(isel)-1
       DO i = ia(isel),i2
          IF (i.EQ.idiag(isel)) CYCLE
          j = ja (i)
          IF(ind(j).GE.0) CYCLE
          IF(ipc(j).EQ.0) CYCLE
          vals=-a(i)
          IF (zerors) THEN
             rsi=0.0e0
             rsj=0.0e0
          ELSE
             rsi=-si(isel)+a(idiag(isel))
             rsj=-si(j)+a(idiag(j))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          epsr=repsmach*ABS(vals)
          IF (sig1.GT.0.0e0) THEN
             del1=rsi
             eta1=rsi+2*sig1
          ELSE
             del1=rsi+2*sig1
             eta1=rsi
          END IF
          IF (eta1.LT.-epsr) CYCLE
          IF (sig2.GT.0.0e0) THEN
             del2=rsj
             eta2=rsj+2*sig2
          ELSE
             del2=rsj+2*sig2
             eta2=rsj
          END IF
          IF (eta2.LT.-epsr) CYCLE
          IF (vals.GT.0.0e0) THEN
             IF (ABS(del1).LT.epsr .AND. ABS(del2).LT.epsr) THEN
               valp=1.0e0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1).LT.epsr) THEN
                IF (del2.LT.-epsr) CYCLE
                valp=1.0e0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2).LT.epsr) THEN
                IF (del1.LT.-epsr) CYCLE
                valp=1.0e0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12.LT.-epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp.LT.0.0e0) CYCLE
                valp=(vals+(eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1.LE.0.0e0 .OR. del2.LE.0.0e0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp.LT.0.0e0) CYCLE
             vals=vals+(eta1*eta2)/(eta1+eta2)
             IF (vals.LT.0.0e0) CYCLE
             valp=vals/valp
          END IF
          IF (valp.GT.kaptg) CYCLE
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. j.LT.ipair)  GOTO 9
          CYCLE
9         CONTINUE
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       IF (ipair.EQ.0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
      ENDDO
    IF (ifirst .AND. 3*nc.GT.2*n .AND. ngl(1).GT.maxcoarseslowt) THEN
       imult=2*imult
       ifirst=.FALSE.
       IF (wfo) THEN
          WRITE (iout,901) IRANK,imult*kaptg_blocdia
       END IF
901    FORMAT(i3, &
       '*   Coarsening too slow: quality threshold increased to',es9.2)
       ind(1:n)=-1
       DO i=1,nddl
          ind(ldd(i))=0
       END DO
       GOTO 1
    ENDIF
    RETURN
  END SUBROUTINE sagmgpar_findpairs_SI
  SUBROUTINE sagmgpar_findpairs_GI1(l,n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,nddl,ldd,skipass                                      &
    ,riperm,iperm &
    ,iext)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: l,n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    INTEGER :: iext(*)
    REAL(kind(0.0e0)) :: a(*)
    REAL(kind(0.0e0)) :: si(n)
    INTEGER :: ldd(nddl),nddl
    LOGICAL :: skipass
    INTEGER :: iperm(n),riperm(n)
    REAL(kind(0.0e0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0e0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2,nt,kb,i0
    LOGICAL :: acc,ifirst
    REAL(kind(0.0e0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
    ifirst=.TRUE.
  1 CONTINUE
    kaptg=imult*kaptg_dampJac
    bndmum1m1=1.0e0/(kaptg-1.0e0)
    dbndmum1=2*1.0e0/kaptg
    umdbndmum1=1.0e0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
      nt=n
      DO WHILE (nmark.LT.nt)
       isel=ijs
       isel=riperm(ijs)
       ijs=ijs+1
       IF (ind(isel).EQ.0) THEN
          idd=idd+1
          nmark=nmark+1
          ldd(idd)=isel
          CYCLE
       END IF
       IF (ind(isel).GE.0) CYCLE
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       IF (iperm(isel).EQ.0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       IF (skipass) THEN
          lcg(2,nc) = -1
          CYCLE
       END IF
       i2=iext(isel)-1
       DO i = ia(isel),i2
          IF (i.EQ.idiag(isel)) CYCLE
          j = ja (i)
          IF(ind(j).GE.0) CYCLE
          IF(iperm(j).EQ.0) CYCLE
          kk=0
          IF (i.LT.idiag(isel)) THEN
             j2=iext(j)-1
             DO jk=idiag(j)+1,j2
                IF (ja(jk).EQ.isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ELSE
             DO jk=ia(j),idiag(j)-1
                IF (ja(jk).EQ.isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ENDIF
          vals=-a(i)/2
          IF(kk.NE.0) vals=vals-a(kk)/2
          IF (zerors) THEN
             rsi=0.0e0
             rsj=0.0e0
             eta1=2*si(isel)
             eta2=2*si(j)
          ELSE
             rsi=-si(isel)+a(idiag(isel))
             rsj=-si(j)+a(idiag(j))
             eta1=2*a(idiag(isel))
             eta2=2*a(idiag(j))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          IF (sig1.GT.0.0e0) THEN
             del1=rsi
          ELSE
             del1=rsi+2*sig1
          END IF
          IF (sig2.GT.0.0e0) THEN
             del2=rsj
          ELSE
             del2=rsj+2*sig2
          END IF
          IF (vals.GT.0.0e0) THEN
             epsr=repsmach*vals
             IF (ABS(del1).LT.epsr .AND. ABS(del2).LT.epsr) THEN
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1).LT.epsr) THEN
                IF (del2.LT.-epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2).LT.epsr) THEN
                IF (del1.LT.-epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12.LT.-epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp.LT.0.0e0) CYCLE
                valp=((eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1.LE.0.0e0 .OR. del2.LE.0.0e0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp.LT.0.0e0) CYCLE
             vals=(eta1*eta2)/(eta1+eta2)
             valp=vals/valp
          END IF
          IF (valp.GT.kaptg) CYCLE
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. iperm(j).LT.iperm(ipair))  GOTO 9
          CYCLE
9         CONTINUE
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       IF (ipair.EQ.0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
      ENDDO
    IF (ifirst .AND. 3*nc.GT.2*n .AND. ngl(1).GT.maxcoarseslowt) THEN
       imult=2*imult
       ifirst=.FALSE.
       IF (wfo) THEN
          WRITE (iout,901) IRANK,imult*kaptg_dampJac
       END IF
901    FORMAT(i3, &
       '*   Coarsening too slow: quality threshold increased to',es9.2)
       ind(1:n)=-1
       DO i=1,nddl
          ind(ldd(i))=0
       END DO
       GOTO 1
    ENDIF
    RETURN
  END SUBROUTINE sagmgpar_findpairs_GI1
  SUBROUTINE sagmgpar_findpairs_SI1(l,n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,nddl,ldd,skipass                                      &
    ,riperm,iperm &
    ,iext)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: l,n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    INTEGER :: iext(*)
    REAL(kind(0.0e0)) :: a(*)
    REAL(kind(0.0e0)) :: si(n)
    INTEGER :: ldd(nddl),nddl
    LOGICAL :: skipass
    INTEGER :: iperm(n),riperm(n)
    REAL(kind(0.0e0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0e0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2,nt,kb,i0
    LOGICAL :: acc,ifirst
    REAL(kind(0.0e0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
    ifirst=.TRUE.
  1 CONTINUE
    kaptg=imult*kaptg_blocdia
    bndmum1m1=1.0e0/(kaptg-1.0e0)
    dbndmum1=2*1.0e0/kaptg
    umdbndmum1=1.0e0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
      nt=n
      DO WHILE (nmark.LT.nt)
       isel=ijs
       isel=riperm(ijs)
       ijs=ijs+1
       IF (ind(isel).EQ.0) THEN
          idd=idd+1
          nmark=nmark+1
          ldd(idd)=isel
          CYCLE
       END IF
       IF (ind(isel).GE.0) CYCLE
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       IF (iperm(isel).EQ.0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       IF (skipass) THEN
          lcg(2,nc) = -1
          CYCLE
       END IF
       i2=iext(isel)-1
       DO i = ia(isel),i2
          IF (i.EQ.idiag(isel)) CYCLE
          j = ja (i)
          IF(ind(j).GE.0) CYCLE
          IF(iperm(j).EQ.0) CYCLE
          vals=-a(i)
          IF (zerors) THEN
             rsi=0.0e0
             rsj=0.0e0
          ELSE
             rsi=-si(isel)+a(idiag(isel))
             rsj=-si(j)+a(idiag(j))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          epsr=repsmach*ABS(vals)
          IF (sig1.GT.0.0e0) THEN
             del1=rsi
             eta1=rsi+2*sig1
          ELSE
             del1=rsi+2*sig1
             eta1=rsi
          END IF
          IF (eta1.LT.-epsr) CYCLE
          IF (sig2.GT.0.0e0) THEN
             del2=rsj
             eta2=rsj+2*sig2
          ELSE
             del2=rsj+2*sig2
             eta2=rsj
          END IF
          IF (eta2.LT.-epsr) CYCLE
          IF (vals.GT.0.0e0) THEN
             IF (ABS(del1).LT.epsr .AND. ABS(del2).LT.epsr) THEN
               valp=1.0e0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1).LT.epsr) THEN
                IF (del2.LT.-epsr) CYCLE
                valp=1.0e0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2).LT.epsr) THEN
                IF (del1.LT.-epsr) CYCLE
                valp=1.0e0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12.LT.-epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp.LT.0.0e0) CYCLE
                valp=(vals+(eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1.LE.0.0e0 .OR. del2.LE.0.0e0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp.LT.0.0e0) CYCLE
             vals=vals+(eta1*eta2)/(eta1+eta2)
             IF (vals.LT.0.0e0) CYCLE
             valp=vals/valp
          END IF
          IF (valp.GT.kaptg) CYCLE
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. iperm(j).LT.iperm(ipair))  GOTO 9
          CYCLE
9         CONTINUE
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       IF (ipair.EQ.0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
      ENDDO
    IF (ifirst .AND. 3*nc.GT.2*n .AND. ngl(1).GT.maxcoarseslowt) THEN
       imult=2*imult
       ifirst=.FALSE.
       IF (wfo) THEN
          WRITE (iout,901) IRANK,imult*kaptg_blocdia
       END IF
901    FORMAT(i3, &
       '*   Coarsening too slow: quality threshold increased to',es9.2)
       ind(1:n)=-1
       DO i=1,nddl
          ind(ldd(i))=0
       END DO
       GOTO 1
    ENDIF
    RETURN
  END SUBROUTINE sagmgpar_findpairs_SI1
  SUBROUTINE sagmgpar_lcgmix(nc,m,lcg1,lcg,lcgn)
    INTEGER :: nc,m,lcg1(m,*),lcg(2,*),lcgn(2*m,*),i,l,l1,l2
    IF (m.eq.2) THEN
       DO i=1,nc
          IF(lcg(2,i).EQ.0) THEN
             lcgn(1,i)=lcg1(1,lcg(1,i))
             lcgn(2,i)=0
             lcgn(4,i)=-1
          ELSE IF(lcg(2,i).LT.0) THEN
             IF (lcg1(2,lcg(1,i)).LT.0) THEN
                lcgn(1,i)=lcg1(1,lcg(1,i))
                lcgn(2,i)=-1
                lcgn(4,i)=-1
             ELSE
                lcgn(1,i)=lcg1(1,lcg(1,i))
                lcgn(2,i)=lcg1(2,lcg(1,i))
                lcgn(4,i)=-2
             END IF
          ELSE
             IF (lcg1(2,lcg(1,i)).LT.0) THEN
                IF (lcg1(2,lcg(2,i)).LT.0) THEN
                   lcgn(1,i)=lcg1(1,lcg(1,i))
                   lcgn(2,i)=lcg1(1,lcg(2,i))
                   lcgn(4,i)=-2
                ELSE
                   lcgn(1,i)=lcg1(1,lcg(1,i))
                   lcgn(2,i)=lcg1(1,lcg(2,i))
                   lcgn(3,i)=lcg1(2,lcg(2,i))
                   lcgn(4,i)=-3
                END IF
             ELSE
                IF (lcg1(2,lcg(2,i)).LT.0) THEN
                   lcgn(1,i)=lcg1(1,lcg(1,i))
                   lcgn(2,i)=lcg1(2,lcg(1,i))
                   lcgn(3,i)=lcg1(1,lcg(2,i))
                   lcgn(4,i)=-3
                ELSE
                   lcgn(1,i)=lcg1(1,lcg(1,i))
                   lcgn(2,i)=lcg1(2,lcg(1,i))
                   lcgn(3,i)=lcg1(1,lcg(2,i))
                   lcgn(4,i)=lcg1(2,lcg(2,i))
                END IF
             END IF
          END IF
       END DO
    ELSE
       DO i=1,nc
          IF(lcg(2,i).EQ.0) THEN
             lcgn(1,i)=lcg1(1,lcg(1,i))
             lcgn(2,i)=0
             lcgn(2*m,i)=-1
          ELSE
             lcgn(2,i)=-1
             l1=m
             IF (lcg1(m,lcg(1,i)).LT.0) l1=-lcg1(m,lcg(1,i))
             lcgn(1:l1,i)=lcg1(1:l1,lcg(1,i))
             IF(lcg(2,i).LT.0) THEN
                l=l1
             ELSE
                l2=m
                IF (lcg1(m,lcg(2,i)).LT.0) l2=-lcg1(m,lcg(2,i))
                lcgn(l1+1:l1+l2,i)=lcg1(1:l2,lcg(2,i))
                l=l1+l2
             END IF
             IF(l.LT.2*m) lcgn(2*m,i)=-l
          END IF
       END DO
    END IF
    RETURN
  END SUBROUTINE sagmgpar_lcgmix
  SUBROUTINE sagmgpar_setind(nc,ndd,ldd,lcg,m,ind)
    INTEGER :: nc,m,lcg(m,*),nll,ldd(ndd),ind(*),i,k,l
    DO i=1,ndd
       ind(ldd(i))=0
    END DO
    DO i=1,nc
       l=m
       IF (lcg(m,i).LT.0) l=-lcg(m,i)
       DO k=1,l
          ind(lcg(k,i))=i
       END DO
    END DO
    RETURN
  END SUBROUTINE sagmgpar_setind
  SUBROUTINE sagmgpar_setcg(n,a,ja,ia,idiag,si,ind,lcg           &
       ,nc,a2,ja2,ia2,idiag2,si2,ysi,maxdg,iw,w,iext,iext2)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),ind(n),nc,lcg(2,nc)
    INTEGER :: ja2(*),ia2(nc+1),idiag2(nc),maxdg
    INTEGER, TARGET :: iw(2*nc)
    INTEGER, OPTIONAL :: iext(*),iext2(*)
    REAL(kind(0.0e0)) :: a(*),a2(*),w(nc),vald
    REAL(kind(0.0e0)) :: si(n),si2(*),sii
    LOGICAL :: ysi
    INTEGER :: nz,nzu,i,jj,jc,jcol,ki,kb,kf,jpos
    INTEGER, POINTER, DIMENSION(:) :: iw1, iw2
    iw1 => iw(1:nc)
    iw2 => iw(nc+1:2*nc)
    nz = 0
    iw1(1:nc)=0
    maxdg=0
    ia2(1)=1
    DO i = 1,nc
       sii=0.0e0
       vald=0.0e0
       nzu=0
       DO ki= 1,2
          jj = lcg(ki,i)
          IF (ki.EQ.1 .OR. jj.GT.0) THEN
             IF (ysi) sii=sii+si(jj)
             kf = iext(jj)-1
             DO kb = ia(jj),kf
                jc = ja(kb)
                jcol = ind(jc)
                IF (jcol.GT.0) THEN
                   IF (jcol.LT.i) THEN
                      jpos = iw1(jcol)
                      IF (jpos.EQ.0) THEN
                         nz = nz+1
                         ja2(nz) = jcol
                         iw1(jcol) = nz
                         a2(nz) = a(kb)
                      ELSE
                         a2(jpos) = a2(jpos) + a(kb)
                      ENDIF
                   ELSE IF (jcol.GT.i) THEN
                      jpos = iw1(jcol)
                      IF (jpos.EQ.0) THEN
                         nzu = nzu+1
                         iw2(nzu) = jcol
                         iw1(jcol) = nzu
                         w(nzu) = a(kb)
                      ELSE
                         w(jpos) = w(jpos) + a(kb)
                      ENDIF
                   ELSE
                      vald=vald+a(kb)
                      IF (ysi .AND. jc.NE.jj) sii=sii+a(kb)
                   ENDIF
                ENDIF
             ENDDO
          ENDIF
       ENDDO
       nz=nz+1
       a2(nz)=vald
       idiag2(i)=nz
       ja2(nz)=i
       a2(nz+1:nz+nzu)=w(1:nzu)
       ja2(nz+1:nz+nzu)=iw2(1:nzu)
       nz=nz+nzu
       maxdg=max(maxdg,nz-ia2(i))
       DO kb = ia2(i), nz
          iw1(ja2(kb))=0
       ENDDO
       IF (ysi) si2(i)=sii
       iext2(i)=nz+1
       DO ki= 1,2
          jj = lcg(ki,i)
          IF (ki.EQ.1 .OR. jj.GT.0) THEN
             DO kb = iext(jj),ia(jj+1)-1
                nz=nz+1
                a2(nz)=a(kb)
                ja2(nz)=ja(kb)
             ENDDO
          ENDIF
       ENDDO
       ia2(i+1)=nz+1
    ENDDO
    RETURN
  END SUBROUTINE sagmgpar_setcg
  SUBROUTINE sagmgpar_setsendrecL1(n,ja,ia,iext,nzext,news,ifl,listrank)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: n, ja(*), ia(n+1), iext(n),nzext,ifl
    INTEGER, OPTIONAL :: listrank(ifl:*)
    LOGICAL :: news
    INTEGER, ALLOCATABLE :: lstot(:),lo(:),ifo(:)
    INTEGER (SELECTED_INT_KIND(2)), ALLOCATABLE :: lproc(:)
    INTEGER :: i,j,ksh,jj,jp,k,kk,ir,ip,ipn,io,ii,ier
    INTEGER :: jmax,jmin,nrowin,llsto,nz,ifim1,maxneigh
    LOGICAL :: firstip,extrank
    maxneigh=MIN(HUGE(lproc),NPROC)
    IF (ALLOCATED(ineigh)) THEN
       ALLOCATE(ilstinL1(maxneigh+1),ilstoutL1(maxneigh+1)    &
               ,ifo(maxneigh),lstot(nzext),lo(nzext),lproc(0:NPROC-1))
    ELSE
       ALLOCATE(ineigh(maxneigh),ilstinL1(maxneigh+1),ilstoutL1(maxneigh+1)    &
               ,ifo(maxneigh),lstot(nzext),lo(nzext),lproc(0:NPROC-1))
    END IF
    IF (news) nneigh=0
    jmax=0
    jmin=HUGE(jmin)
    io=0
    ii=0
    nrowin=0
    lproc=0
    DO i=1,n
       IF (iext(i).LT.ia(i+1)) nrowin=nrowin+1
       DO kk=iext(i),ia(i+1)-1
          j=ja(kk)
          jmin=MIN(j,jmin)
          jmax=MAX(j,jmax)
          ir=listrank(j)
          firstip=.FALSE.
          IF (ir.GT.NPROC) THEN
             ip=ir-NPROC
          ELSE
             ip=lproc(ir)
             IF (ip.EQ.0) THEN
                IF (news) THEN
                   nneigh=nneigh+1
                   IF (nneigh.GT.maxneigh) THEN
                      WRITE (6,1020) IRANK,maxneigh
                      CALL sagmgpar_ABORT
                   END IF
                   ip=nneigh
                   ineigh(ip)=ir
                ELSE
                   DO k=1,nneigh
                      IF (ineigh(k).EQ.ir) EXIT
                   END DO
                   IF (k.GT.nneigh) THEN
                      WRITE (6,1010) IRANK
                      CALL sagmgpar_ABORT
                   END IF
                   ip=k
                END IF
                lproc(ir)=ip
                ilstinL1(ip+1)=1
                firstip=.TRUE.
             ELSE
                ilstinL1(ip+1)=ilstinL1(ip+1)+1
             END IF
             listrank(j)=ip+NPROC
          END IF
          IF (firstip) THEN
             io=io+1
             lstot(io)=i
             lo(io)=0
             ifo(ip)=io
             ilstoutL1(ip)=io
          ELSE IF (lstot(ilstoutL1(ip)) /= i) THEN
             io=io+1
             lstot(io)=i
             lo(io)=0
             lo(ilstoutL1(ip))=io
             ilstoutL1(ip)=io
          END IF
       END DO
    END DO
    llsto=io
    ALLOCATE(lstoutL1(llsto))
    io=1
    ilstoutL1(1)=1
    ilstinL1(1)=1
    DO i=1,nneigh
      ip=ifo(i)
      DO WHILE (ip /= 0)
          lstoutL1(io)=lstot(ip)
          ip=lo(ip)
          io=io+1
       END DO
       ilstoutL1(i+1)=io
       ilstinL1(i+1)=ilstinL1(i)+ilstinL1(i+1)
       ifo(i)=ilstinL1(i)
    END DO
    DO j=jmin,jmax
       ip=listrank(j)-NPROC
       IF (ip.GT.0 .AND. ip.LE.nneigh) THEN
          listrank(j)=ifo(ip)+n
          ifo(ip)=ifo(ip)+1
       END IF
    END DO
    DEALLOCATE(ifo,lo,lstot,lproc)
    IF (ALLOCATED(ireqi)) THEN
       ALLOCATE(lstinL1(0:nrowin))
       IF (nneigh.GT.UBOUND(ireqi,1)) THEN
          DEALLOCATE(ireqi,ireqo)
          ALLOCATE(ireqi(nneigh),ireqo(nneigh))
       END IF
       IF (ilstinL1(nneigh+1)-1.GT.UBOUND(buffi,1)) THEN
          DEALLOCATE(buffi)
          ALLOCATE(buffi(ilstinL1(nneigh+1)-1))
       END IF
       IF (ilstoutL1(nneigh+1)-1.GT.UBOUND(buffo,1)) THEN
          DEALLOCATE(buffo)
          ALLOCATE(buffo(ilstoutL1(nneigh+1)-1))
       END IF
    ELSE
       ALLOCATE(lstinL1(0:nrowin),                          &
         ireqi(nneigh),ireqo(nneigh),                       &
         buffi(ilstinL1(nneigh+1)-1),buffo(ilstoutL1(nneigh+1)-1))
    END IF
    lstinL1(0)=nrowin
    nrowin=0
    DO i=1,n
       IF (iext(i).LT.ia(i+1)) THEN
          nrowin=nrowin+1
          lstinL1(nrowin)=i
       END IF
       DO kk=iext(i),ia(i+1)-1
          ja(kk)=listrank(ja(kk))
       END DO
    END DO
    RETURN
1001 FORMAT(i3,'*',                                                   &
         ' FATAL ERROR; likely: listrank not properly defined')
1010 FORMAT(' The new input matrix has on rank',i9, 'has external'      &
         ,/,' connection(s) with processor(s) (MPI rank(s)) that were'  &
         ,/,' not connected to it while the preconditioner was built.'  &
         ,/,' THIS IS NOT ALLOWED: ABORTING !!!'                        &
         ,/,' Please consider rebuilding of the preconditioner in such cases')
1020 FORMAT(' The matrix rows on rank',i9, 'have external connections'  &
         ,/,' with more than',i5,' different processors (MPI ranks).'   &
         ,/,' This is above standard expectation and'                   &
         ,/,' NOT ALLOWED with the present version of AGMG'            &
         ,/,' ... ABORTING !!!'                                         &
         ,/,' Please contact support@agmg.eu to solve this issue')
  END SUBROUTINE sagmgpar_setsendrecL1
  SUBROUTINE sagmgpar_setextnum(nc,l,iw,imax)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: nc,l,iw(nc),imax
    INTEGER :: i,ii,j,jj,jp,ic,ip,ict,ier,idest,kdeb,kfin,i1,k,ksh
    LOGICAL :: mv
    IF (nneigho.GT.0) CALL MPI_WAITALL(nneigho,ireqo,MPI_STATUSES_IGNORE,ier)
    IF (nc.EQ.0) THEN
       ii=0
       DO i=1,nneigh
          kdeb=dt(l)%ilstout(i)
          kfin=dt(l)%ilstout(i+1)-1
          IF (kfin.GE.kdeb) THEN
             buffo(kdeb:kfin)=-1.0e0
             ii=ii+1
             idest=ineigh(i)
             CALL MPI_ISEND(buffo(kdeb),kfin-kdeb+1,MPI_REAL,&
                  idest,989,icomm,ireqo(ii),ier)
          END IF
       END DO
       nneigho=ii
    ELSE
       ALLOCATE(dt(l+1)%ilstout(nneigh+1),dt(l+1)%ilstin(nneigh+1)   &
            ,dt(l+1)%lstout(dt(l)%ilstout(nneigh+1)-1) )
       iw(1:nc)=0
       dt(l+1)%ilstout(1)=1
       ict=0
       ii=0
       DO i=1,nneigh
          ic=0
          kdeb=dt(l)%ilstout(i)
          kfin=dt(l)%ilstout(i+1)-1
          IF (kfin.GE.kdeb) THEN
             ii=ii+1
             mv=.TRUE.
             DO j=kdeb,kfin
                jj=dt(l)%ind(dt(l)%lstout(j))
                IF (jj.GT.0) THEN
                   IF (iw(jj).EQ.0) THEN
                      ic=ic+1
                      buffo(j)=real(ic)
                      iw(jj)=ic
                      dt(l+1)%lstout(ic+ict)=jj
                   ELSE
                      buffo(j)=real(iw(jj))
                   END IF
                   mv=.FALSE.
                ELSE
                   buffo(j)=0.0e0
                END IF
             END DO
             IF (mv) buffo(kfin)=-1.0e0
             idest=ineigh(i)
             CALL MPI_ISEND(buffo(kdeb),kfin-kdeb+1,MPI_REAL,&
                  idest,989,icomm,ireqo(ii),ier)
             DO j=ict+1,ict+ic
                iw(dt(l+1)%lstout(j))=0
             END DO
          END IF
          ict=ic+ict
          dt(l+1)%ilstout(i+1)=ict+1
       END DO
       nneigho=ii
    END IF
    ii=0
    DO i=1,nneigh
       kdeb=dt(l)%ilstin(i)
       kfin=dt(l)%ilstin(i+1)-1
       IF (kfin.GE.kdeb) THEN
          ii=ii+1
          idest=ineigh(i)
          CALL MPI_IRECV(buffi(kdeb),kfin-kdeb+1,MPI_REAL,&
               idest,989,ICOMM,ireqi(ii),ier)
       END IF
    END DO
    nneighi=ii
    IF (nneighi.GT.0) CALL MPI_WAITALL(nneighi,ireqi,MPI_STATUSES_IGNORE,ier)
    IF (nc.EQ.0) RETURN
    IF (nc.EQ.0) THEN
       IF (dt(l)%ilstout(1).NE.0) DEALLOCATE(dt(l+1)%lstout)
       dt(l+1)%ilstout=0
       dt(l+1)%ilstin=0
    ELSE
     ict=0
     mv=.FALSE.
     dt(l+1)%ilstin(1)=1
     ip=dt(l+1)%ilstout(1)
     DO i=1,nneigh
       ic=0
       kdeb=dt(l)%ilstin(i)
       kfin=dt(l)%ilstin(i+1)-1
       ii=0
       IF (kfin.GE.kdeb .AND. dt(l+1)%ilstout(i+1).GT.ip ) &
          THEN
          DO j=kdeb,kfin
             jj=NINT(buffi(j))
             IF (jj.GT.0) THEN
                ic=max(ic,jj)
                buffi(j)=ict+jj+nc
             END IF
          END DO
          IF (jj.LT.0) THEN
             mv=.TRUE.
          ELSE IF (mv) THEN
             ii=dt(l+1)%ilstout(i+1)-ip
          END IF
       END IF
       ict=ict+ic
       dt(l+1)%ilstin(i+1)=ict+1
!!#ifdef
       IF (mv) THEN
          i1=dt(l+1)%ilstout(i)
          dt(l+1)%lstout(i1:i1+ii-1)=dt(l+1)%lstout(ip:ip+ii-1)
          ip=dt(l+1)%ilstout(i+1)
          dt(l+1)%ilstout(i+1)=dt(l+1)%ilstout(i)+ii
       ELSE
          ip=dt(l+1)%ilstout(i+1)
       END IF
     END DO
    END IF
    RETURN
  END SUBROUTINE sagmgpar_setextnum
  SUBROUTINE sagmgpar_restorelstrank(ilstin,listrank,ifl)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: listrank(*),ilstin(nneigh+1),ifl,i,j,ir
    DO i=1,nneigh
       ir=ineigh(i)
       DO j=ilstin(i),ilstin(i+1)-1
          listrank(j)=ir
       END DO
    END DO
    RETURN
  END SUBROUTINE sagmgpar_restorelstrank
  SUBROUTINE sagmgpar_setcgext(l,n,a,ja,ia,idiag,iext,nfine             &
                            ,ao,jao,iw,nout)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER :: l,n, ja(*), ia(n+1), idiag(n), iext(n), jao(*)
    INTEGER :: nfine, nout, iw(n+1:n+nout)
    REAL(kind(0.0e0)) :: a(*), ao(*)
    INTEGER :: i, k, jj, jc, kb, jcol, ipos, iposi, jpos, nrowin
    iw(n+1:n+nout)=0
    ipos=1
    nrowin=0
    DO i = 1, n
       iposi=ipos
       idiag(i)=idiag(i)-(ia(i)-ipos)
       DO k=ia(i),iext(i)-1
          ao(ipos)=a(k)
          jao(ipos)=ja(k)
          ipos=ipos+1
       END DO
       ia(i)=iposi
       iposi=ipos
       DO k=iext(i),ia(i+1)-1
          jcol=NINT(buffi(ja(k)-nfine))
          IF (jcol.GT.0) THEN
             jpos = iw(jcol)
             IF (jpos.EQ.0) THEN
                jao(ipos) = jcol
                iw(jcol) = ipos
                ao(ipos) = a(k)
                ipos=ipos+1
             ELSE
                ao(jpos) = ao(jpos) + a(k)
             END IF
          END IF
       END DO
       DO k=iposi,ipos-1
          iw(jao(k)) = 0
       END DO
       iext(i)=iposi
       IF (ipos.GT.iposi) nrowin=nrowin+1
    END DO
    ia(n+1)=ipos
    ALLOCATE(dt(l+1)%lstin(0:nrowin))
    dt(l+1)%lstin(0)=nrowin
    nrowin=0
    DO i=1,n
       IF (iext(i).LT.ia(i+1)) THEN
          nrowin=nrowin+1
          dt(l+1)%lstin(nrowin)=i
       END IF
    END DO
    RETURN
  END SUBROUTINE sagmgpar_setcgext
  SUBROUTINE sagmgpar_MUMPSpar(n,f,ijob,a,ja,ia)
    USE sagmgpar_mem
    IMPLICIT NONE
    INCLUDE 'smumps_struc.h'
    INTEGER :: n,ijob
    REAL(kind(0.0e0)), TARGET :: f(n)
    REAL(kind(0.0e0)), OPTIONAL, TARGET :: a(*)
    INTEGER, OPTIONAL :: ia(n+1)
    INTEGER, OPTIONAL, TARGET :: ja(*)
    INTEGER(selected_int_kind(8)), SAVE :: nglobl, nnn
    REAL(kind(0.0e0)), SAVE :: iflop
    INTEGER(selected_int_kind(8)), ALLOCATABLE, SAVE :: nloc(:), idispl(:)
    TYPE(SMUMPS_STRUC), SAVE :: mumps_par
    INTEGER :: ierr, i, j, k
    nnn=n
    IF (ijob.EQ.-2) THEN
       mumps_par%JOB = -2
       CALL SMUMPS(mumps_par)
       DEALLOCATE(mumps_par%RHS,nloc,idispl,mumps_par%IRN_loc,mumps_par%JCN_loc)
       IF (nlev.EQ.1) DEALLOCATE(mumps_par%A_loc)
    ELSE IF (ijob.EQ.1) THEN
       mumps_par%COMM = ICOMM
       mumps_par%JOB = -1
       mumps_par%SYM = 0
       mumps_par%PAR = 1
       CALL SMUMPS(mumps_par)
       mumps_par%ICNTL(2)=-1
       mumps_par%ICNTL(3)=-1
       mumps_par%ICNTL(4)=0
       mumps_par%ICNTL(14)=80
       mumps_par%ICNTL(18)=3
       mumps_par%ICNTL(24)=1
       IF (n.GT.0) THEN
          mumps_par%NZ_loc=ia(n+1)-ia(1)
          ALLOCATE(nloc(NPROC),mumps_par%IRN_loc(mumps_par%NZ_loc),   &
                mumps_par%JCN_loc(mumps_par%NZ_loc))
          mumps_par%JCN_loc = ja(1:mumps_par%NZ_loc)
          CALL sagmgpar_globnum(nnn,ja,ia,dt(nlev)%iext   &
            ,mumps_par%IRN_loc,mumps_par%JCN_loc        &
            ,nlev,dt(nlev)%lstin,nloc)
          IF (nlev.GT.1) THEN
             mumps_par%A_loc => dt(nlev)%a(1:mumps_par%NZ_loc)
          ELSE
             ALLOCATE ( mumps_par%A_loc(mumps_par%NZ_loc) )
             mumps_par%A_loc = a(1:mumps_par%NZ_loc)
          END IF
       ELSE
          mumps_par%NZ_loc=0
          ALLOCATE(mumps_par%IRN_loc(0),mumps_par%JCN_loc(0)        &
                  ,mumps_par%A_loc(0),nloc(NPROC))
          CALL MPI_ALLGATHER(n,1,MPI_INTEGER,nloc,1,MPI_INTEGER,ICOMM,ierr)
         nneighi=0
       END IF
       IF (mumps_par%MYID.EQ.0) THEN
          ALLOCATE(idispl(NPROC))
          IF(IRANK /= 0) THEN
             PRINT *, 'MYID in MUMPS /= IRANK, ABORTING'
             CALL sagmgpar_ABORT
          END IF
          idispl(1)=0
          nglobl=nnn
          DO i=2,NPROC
             idispl(i)=nglobl
             nglobl=nglobl+nloc(i)
          END DO
          mumps_par%N=nglobl
          ALLOCATE( mumps_par%RHS(mumps_par%N) )
       ELSE
          ALLOCATE(mumps_par%RHS(1),idispl(1))
       END IF
       mumps_par%JOB = 4
       CALL SMUMPS(mumps_par)
       flop=mumps_par%RINFO(3)
       i=mumps_par%INFO(27)
       IF (i.GT.0) THEN
          iflop=4*real(i)-n
       ELSE
          iflop=-4.0d6*real(i)-n
       END IF
       IF (nneigho.GT.0) CALL MPI_WAITALL(nneigho,ireqo,MPI_STATUSES_IGNORE,ierr)
    ELSE IF (ijob.EQ.2) THEN
       mumps_par%JOB = 3
       CALL MPI_GATHERV(f,n,MPI_REAL,mumps_par%RHS      &
            ,nloc,idispl,MPI_REAL,0,ICOMM,ierr)
       CALL SMUMPS(mumps_par)
       CALL MPI_SCATTERV(mumps_par%RHS,nloc,idispl     &
            ,MPI_REAL,f,n,MPI_REAL,0,ICOMM,ierr)
       flop=flop+iflop
    END IF
    RETURN
  END SUBROUTINE sagmgpar_MUMPSpar
  SUBROUTINE sagmgpar_globnum(n,ja,ia,iext,irg,jcg,l,lstin,nloc)
    USE sagmgpar_mem
    IMPLICIT NONE
    INTEGER(selected_int_kind(8)) :: n,nloc(*)
    INTEGER :: ia(n+1),ja(*),iext(n),irg(*),jcg(*)
    INTEGER :: l,lstin(0:*)
    REAL(kind(0.0e0)) :: dum(1)
    INTEGER :: i,ii,iii,k,ier,ishift,jj
    CALL MPI_ALLGATHER(n,1,MPI_INTEGER,nloc,1,MPI_INTEGER,ICOMM,ier)
    IF (n.LE.0) THEN
       nneighi=0
       RETURN
    END IF
    ishift=0
    DO i=0,IRANK-1
       ishift=ishift+nloc(i+1)
    END DO
    DO i=1,n
      irg(i)=i
    END DO
    CALL sagmgpar_sendrec(dum,irg,ishift,l)
    DO i=1,n
       DO k=ia(i),iext(i)-1
          irg(k)=i+ishift
          jcg(k)=ja(k)+ishift
       END DO
    END DO
    IF (nneighi.GT.0) CALL MPI_WAITALL(nneighi,ireqi,MPI_STATUSES_IGNORE,ier)
    DO iii=1,lstin(0)
       i=lstin(iii)
       DO k=iext(i),ia(i+1)-1
          irg(k)=i+ishift
          jcg(k)=NINT(buffi(ja(k)-n))
       END DO
    END DO
    RETURN
  END SUBROUTINE sagmgpar_globnum
END MODULE sagmgpar_ALLROUTINES
!!!!!!!!!!!!!!!!!!! MAIN DRIVER
  SUBROUTINE sagmgpar( n,a,ja,ia,f,x,ijob,iprint,nrest,iter,tol &
                 , MPI_COMM_AGMG ,listrank,ifirstlistrank )
    USE sagmgpar_mem
    USE sagmgpar_ALLROUTINES
    IMPLICIT NONE
    INTEGER    :: n,ia(n+1),ijob,iprint,nrest,iter
    INTEGER  :: ja(*)
    REAL (kind(0.0e0)) :: a(*),f(n),x(n)
    REAL (kind(0.0e0)) :: tol
    INTEGER    :: MPI_COMM_AGMG ,listrank(*),ifirstlistrank
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Arguments
!  =========
!
!  N       (input) INTEGER.
!          The number of rows in the local part of the matrix.
!
!  A       (input/output) REAL (kind(0.0e0)). Numerical values of the matrix
!                                                (local rows).
!  IA      (input/output) INTEGER. Pointers for every local row.
!  JA      (input/output) INTEGER. Column indices (w.r.t. local numbering).
!
!          Significant only if IJOB==0,1,2,10,12
!
!          Detailed description of the matrix format
!
!              On input, IA(I), I=1,...,N, refers to the physical start
!              of row I. That is, the entries of row I are located
!              in A(K), where K=IA(I),...,IA(I+1)-1. JA(K) carries the
!              associated column indices. IA(N+1) must be defined in such
!              a way that the above rule also works for I=N (that is,
!              the last valid entry in arrays A,JA should correspond to
!              index K=IA(N+1)-1).
!
!          AGMG ASSUMES THAT ALL DIAGONAL ENTRIES ARE POSITIVE
!
!              That is, AGMG assumes that, for K=IA(I),...,IA(I+1)-1,
!              some of JA(K) is equal to I with corresponding A(K)>0
!              (value of the diagonal element, which must be positive).
!
!              A,IA,JA are "output" parameters because on exit the
!              entries of each row may occur in a different order
!              (The matrix is mathematically the same, but stored in
!              different way).
!
!
!  F       (input/output) REAL (kind(0.0e0)).
!          On input, the right hand side vector f.
!          Overwritten on output.
!          Significant only if IJOB is not equal to 1 or 101.
!
!  X       (input/output) REAL (kind(0.0e0)).
!          On input and if IJOB==10, 12, 212: initial guess
!          On output, the computed solution
!              (IJOB==3: result of the application of the preconditioner)
!
! IJOB     (input) INTEGER. Tells AGMG what has to be done.
!          0: performs setup + solve + memory release, no initial guess
!         10: performs setup + solve + memory release, initial guess in x(1:n)
!          1: performs setup only
!             (preprocessing: prepares all parameters for subsequent solves)
!          2: solves only with the preconditioner based on previous setup,
!             using the system matrix provided in A, JA, IA; no initial guess
!         12: solves only with the preconditioner based on previous setup,
!             using the system matrix provided in A, JA, IA;
!             initial guess in x(1:n)
!        202: solves only with the preconditioner based on previous setup and
!             the same matrix as the one provided for set up;  no initial guess
!             using the system matrix provided during previous setup
!        212: solves only with the preconditioner based on previous setup and
!             the same matrix as the one provided for set up;
!             initial guess in x(1:n)
!          3: the vector returned in x(1:n) is not the solution of the linear
!                 system, but the result of the action of the multigrid
!                 preconditioner to the right hand side in f(1:n)
!         -1: erases the setup and releases internal memory
!
!   !!! IJOB==2,3,12,202,212 require that one has previously called
!       AGMG with IJOB==1
!
!   !!! (change with respect to versions 2.x) !!!
!       The preconditioner defined when calling AGMG
!         with IJOB==1 is entirely kept in internal memory.
!       Hence the arrays A, JA and IA are not accessed upon subsequent calls
!         with IJOB==3 and IJOB==202,212
!       Upon subsequent calls with IJOB==2, 12:  a matrix needs to
!            be supplied in arrays A, JA, IA. It will be used to
!            perform matrix vector product within the main iterative
!            loop (and only for this).
!            Hence the system is solved with this matrix which
!            may differ from the matrix in A, JA, IA that was supplied
!            upon the previous call with IJOB==1;
!            then AGMG attempts to solve a linear system with the "new"
!            matrix (supplied when IJOB==2,12) using the preconditioner
!            set up for the "old" one (supplied when IJOB==1).
!       These functionalities (set up a preconditioner and use it for another
!            matrix) are provided for the sake of generality but should be
!            used with care; in general, set up is fast with AGMG and hence
!            it is recommended to rerun it even if the matrix changes only
!            slightly.
!       In the more standard case where one needs successive solves with a
!       same matrix, it recommended to use IJOB==202 or 212 (in which case
!       the arrays A, JA, IA are not accessed) instead of IJOB==2 or 12
!       (in which case the input matrix has to be kept equal to the one
!       provided for set up). This  slightly faster and helps to save memory
!       since the arrays that store the input matrix may be deallocated after
!       set up.
!
! IPRINT   (input) INTEGER.
!              Indicates the unit number where information is to be printed
!              (N.B.: 5 is converted to 6). If nonpositive, only error
!              messages are printed on standard output. Warning messages about
!              insufficient convergence (in the allowed  maximum number of
!              iterations) are further suppressed supplying a negative number.
!
! NREST    (input) INTEGER.
!             Restart parameter for GCR (an implementation of GMRES)
!             Nonpositive values are converted to NREST=10 (default)
!
! !!  If NREST==1, Flexible CG is used instead of GCR (when IJOB==0,10,2, 12)
!             and also (IJOB==0,1) performs some simplification during
!             the set up based on the assumption that the matrix
!             supplied in A, JA, IA is symmetric.
!
! !!!  NREST=1 Should be used if and only if the matrix is really SYMMETRIC
! !!!         (and positive definite).
!
!  ITER    (input/output) INTEGER.
!          On input, the maximum number of iterations. Should be positive.
!          On output, actual number of iterations.
!            If this number of iteration was insufficient to meet convergence
!            criterion, ITER will be returned negative and equal to the
!            opposite of the number of iterations performed.
!          Significant only if IJOB==0, 2, 10, 12, 202, 212
!
!  TOL     (input) REAL (kind(0.0e0))
!          The tolerance on residual norm. Iterations are stopped whenever
!               || A*x-f || <= TOL* || f ||
!          Should be positive and less than 1.0
!          Significant only if IJOB==0, 2, 10, 12, 202, 212
!
! MPI_COMM_AGMG (input) INTEGER
!               MPI communicator
!
! LISTRANK(IFIRSTLISTRANK:*) INTEGER (input/output)
!          Contains the rank of the tasks to which rows corresponding to
!          nonlocal variable referenced in (A,JA,IA) are assigned.
!          Let Jmin and Jmax be, respectively, the smallest and the largest
!          index of a nonlocal variable referenced in JA(1:IA(N+1)-1)
!          (Note that N < Jmin <= Jmax).
!          listrank(i) will be referenced if and only if  Jmin <= i <= Jmax,
!          and listrank(i) should then be equal to the rank of the "home"
!          task of i if i is effectively present in JA(1:IA(N+1)-1),
!          and equal to some arbitrary NEGATIVE integer otherwise.
!
!          listrank and ifirstlistrank may be modified on output, according
!          to the possible modification of the indexes of nonlocal variables;
!          that is, on output, listrank and ifirstlistrank still carry the
!          correct information about nonlocal variables, but for the
!          matrix as defined in (A,JA,IA) on output.
!
!
!!!!! Remark !!!! Except insufficient number of iterations to achieve
!                 convergence (characterized by a negative value returned
!                 in ITER), all other detected errors are fatal and lead
!                 to a STOP statement.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER, SAVE :: nza
    INTEGER :: i,j,ijb,k,llev,nzext,itef,ierr,errc
    REAL(kind(0.0e0)) :: cputmp
    REAL(kind(0.0e0)) :: resid,tmp
    REAL(kind(0.0e0)), POINTER, DIMENSION(:) :: ascaled
    REAL(kind(0.0e0)) :: adum(1)
    CHARACTER*40 :: valenv
    wfo=.TRUE.
    iout=iprint
    IF (iprint.LE.0) THEN
       iout=6
       wfo=.FALSE.
    ELSE IF (iprint.EQ.5) THEN
       iout=6
    END IF
    ijb=MOD(ijob,100)
    itef=iter
    CALL sagmgpar_init(MPI_COMM_AGMG)
    IF (ijob.NE.-1 .AND. &
        ijob.NE. 0 .AND. &
        ijob.NE. 1 .AND. &
        ijob.NE. 2 .AND. &
        ijob.NE. 3 .AND. &
        ijob.NE. 4 .AND. &
        ijob.NE. 5 .AND. &
        ijob.NE. 6 .AND. &
        ijob.NE. 7 .AND. &
        ijob.NE. 8 .AND. &
        ijob.NE. 9 .AND. &
        ijob.NE.10 .AND. &
        ijob.NE.12 .AND. &
        ijob.NE.100 .AND. &
        ijob.NE.101 .AND. &
        ijob.NE.102 .AND. &
        ijob.NE.110 .AND. &
        ijob.NE.112 .AND. &
        ijob.NE.202 .AND. &
        ijob.NE.212  &
        ) THEN
       WRITE (6,1000) IRANK,ijob
       CALL MPI_ABORT(MPI_COMM_AGMG,errc,ierr)
       STOP
    END IF
1000 FORMAT(i3,'*',' FATAL ERROR: ijob=',i5, &
         ' is not allowed')
    IF (ijob.GE.100 .AND. ijob.LE.201) THEN
       WRITE (6,1002) IRANK,ijob
       CALL MPI_ABORT(MPI_COMM_AGMG,errc,ierr)
       STOP
    END IF
1002 FORMAT(i3,'*',' FATAL ERROR: ijob=',i4,       &
         ' (i.e. >= 100 & < 200: work with transpose) ',/, &
         '     is not allowed with the parallel (OpenMP or MPI) version')
    IF (MOD(ijob,10).GE.2 .AND. .NOT.preprocessed) THEN
       WRITE (6,1001) IRANK,ijob
       CALL MPI_ABORT(MPI_COMM_AGMG,errc,ierr)
       STOP
    END IF
1001 FORMAT(i3,'*',' FATAL ERROR: setup not done: ijob=',i3, &
         ' is not allowed')
    IF (iprint.EQ.-978 .AND. ijob.GE.0) THEN
       valenv(1:10)='AGMG_dump_'
       valenv(19:19)='_'
       CALL DATE_AND_TIME(valenv(11:18),valenv(20:29))
       valenv(28:31)='.mat'
       valenv(32:32)='.'
       WRITE(valenv(33:38),'(i6.6)') IRANK
       OPEN(UNIT=91825,FILE=valenv,FORM='FORMATTED')
       IF (ijob.LE.200) THEN
          WRITE(91825,'(2i18,1x,A4,i1,2i5,i8,ES18.10)')                       &
               n,ia(n+1)-1,'Dump',0,ijob,nrest,iter,tol
         DO i=1,n
            DO j=ia(i),ia(i+1)-1
               WRITE(91825,'(2i18,ES18.10)') i,ja(j),a(j)
            END DO
         END DO
       ELSE
          WRITE(91825,*) n,0,'DUMP',ijob,nrest,iter,tol
       END IF
       CLOSE (91825)
       IF (MOD(ijob,10).NE.1) THEN
          valenv(29:31)='rhs'
          OPEN (UNIT=91825,FILE=valenv,FORM='FORMATTED')
          DO i=1,n
             WRITE(91825,'(ES18.10)') f(i)
          END DO
          CLOSE(91825)
          valenv(29:31)='ini'
          OPEN (UNIT=91825,FILE=valenv,FORM='FORMATTED')
          DO i=1,n
             WRITE(91825,'(ES18.10)') x(i)
          END DO
          CLOSE(91825)
       END IF
       valenv(29:31)='out'
       OPEN (UNIT=91825,FILE=valenv,FORM='FORMATTED')
       iout=91825
       WRITE(iout,910) IRANK
       wfo=.TRUE.
    END IF
    IF (ijob.LT.0 .AND. .NOT.preprocessed) RETURN
    IF (iprint.GT.0 .AND. ijob.GE.0 .AND. MOD(ijob,10).LE.1) THEN
       WRITE(iout,910) IRANK
    END IF
 910 FORMAT(i3,'*ENTERING AGMG **********************************',&
         '***************************')
    wff=wfo.AND.(IRANK.LE.0)
    woo=IRANK.LE.0
    IF (iprint.LT.0 .AND. iprint.NE.-978) woo=.FALSE.
    IF (ijb.LT.0) GOTO 500
    IF (ijob.LE.201) trans=ijob.GE.100
    IF (MOD(ijb,10).GE.2) THEN
      IF (ijob.LE.201 .AND. MOD(ijob,10).LE.2) THEN
        ALLOCATE(dt(1)%idiag(n+1),iextL1(n+1))
        CALL sagmgpar_partroword(n,a,ja,ia,dt(1)%idiag,iextL1,nzext)
        CALL sagmgpar_setsendrecL1(n,ja,ia,iextL1,nzext                        &
                                ,.FALSE.,ifirstlistrank,listrank)
      CALL sagmgpar_restorelstrank(ilstinL1,listrank,ifirstlistrank)
      ifirstlistrank=n+1
        DEALLOCATE(dt(1)%idiag)
      END IF
    ELSE
      IF (preprocessed) THEN
       CALL sagmgpar_relmem(-1)
       solve=.FALSE.
       eltm=0.0e0
       cputm=0.0e0
       flop=0.0e0
       kstat=0
      END IF
      CALL sagmgpar_mestime(-1,0.0e0,0.0e0)
      spd=nrest.EQ.1
      spdcoarse=spd
      transint=trans.AND.(.NOT.spd)
      nlev=0
      imult=1
      ALLOCATE(dt(1)%idiag(n+1),iextL1(n+1))
      CALL sagmgpar_partroword(n,a,ja,ia,dt(1)%idiag,iextL1,nzext)
      CALL sagmgpar_setsendrecL1(n,ja,ia,iextL1,nzext                          &
                              ,.TRUE.,ifirstlistrank,listrank)
      dt(1)%iext => iextL1
      dt(1)%lstin   => lstinL1
      dt(1)%ilstin  => ilstinL1
      dt(1)%lstout  => lstoutL1
      dt(1)%ilstout => ilstoutL1
      nza=ia(n+1)-ia(1)
      CALL sagmgpar_setupL1(n,a,ja,ia)
      CALL sagmgpar_restorelstrank(ilstinL1,listrank,ifirstlistrank)
      ifirstlistrank=n+1
      CALL sagmgpar_smoothsetup
      preprocessed=.TRUE.
      CALL sagmgpar_mestime(1,cputmp,eltmp)
      IF(wfo)THEN
   !CPU_TIME: next line may be uncommented if implemented
   !          (see subroutine mestime)
   !              WRITE(iout,996) IRANK,cputmp
       WRITE(iout,997) IRANK,eltmp
       WRITE(iout,'()')
      END IF
      IF (MOD(ijb,10).EQ.1) THEN
        IF (nlev.EQ.1) THEN
          ALLOCATE(dt(1)%a(ia(n+1)-ia(1)),dt(1)%ja(ia(n+1)-ia(1)),dt(1)%ia(n+1))
          dt(1)%a=a(ia(1):ia(n+1)-1)
          dt(1)%ja=ja(ia(1):ia(n+1)-1)
          dt(1)%ia=ia(1:n+1)
        END IF
        IF (nlev.GT.1) THEN
          DEALLOCATE(iextL1)
        ELSE
          NULLIFY(iextL1)
        END IF
        NULLIFY(lstinL1, ilstinL1, lstoutL1, ilstoutL1)
        IF (iprint.EQ.-978) CLOSE(91825)
        RETURN
      END IF
    END IF
    CALL sagmgpar_mestime(-2,0.0e0,0.0e0)
    resid=tol
    IF (MOD(ijb,10).GE.3) THEN
       IF (wfo) THEN
          WRITE(iout,901) IRANK
       END IF
       IF(.NOT.solve) kstat(1,MAX(2,nlev))=1
       CALL sagmgpar_applyprec(n,f,x,MOD(ijb,10))
       CALL sagmgpar_mestime(2,cputmp,eltmp)
       cputm=cputm+cputmp
       eltm=eltm+eltmp
       solve=.TRUE.
       IF (iprint.EQ.-978) CLOSE(91825)
       RETURN
    ELSE IF (nrest.NE.1) THEN
       nrst=nrest
       IF (nrst.LE.0) nrst=10
       CALL sagmgpar_GCR(n,f,x,itef,resid,a,ja,ia,ijob)
    ELSE
       CALL sagmgpar_FlexCG(n,f,x,itef,resid,a,ja,ia,ijob)
    END IF
    iter=itef
    solve=.FALSE.
    CALL sagmgpar_mestime(2,cputm,eltm)
    IF (wfo) THEN
       IF (wff .AND. iter.NE.0) THEN
          llev=max(2,nlev)
          IF (allzero) llev=nlev-1
          DO i=2,llev
             WRITE(iout,955) min(i,nlev),kstat(2,i-1),kstat(2,i),       &
                  real(kstat(2,i))/real(kstat(2,i-1)),kstat(1,i)
          END DO
       END IF
       WRITE(iout,'()')
       IF (RESID.GT.0 .AND. iter.NE.0) THEN
         tmp=(flop*log(1.0e0/10)/log(RESID))
         WRITE(iout,952) IRANK,tmp/real(2*nza)
       END IF
   !CPU_TIME: next line may be uncommented if implemented
   !          (see subroutine mestime)
   !            WRITE(iout,998) IRANK,cputm
       WRITE(iout,999) IRANK,eltm
       WRITE(iout,'()')
       IF (RESID.GT.0 .AND. iter.NE.0) THEN
          WRITE (iout,902) IRANK
       END IF
    END IF
    IF (MOD(ijb,10).GT.0) THEN
       solve=.FALSE.
       eltm=0.0e0
       cputm=0.0e0
       flop=0.0e0
       kstat=0
       IF (ijob.LE.201) THEN
         DEALLOCATE(iextL1, lstinL1, ilstinL1, lstoutL1, ilstoutL1)
       END IF
       IF (iprint.EQ.-978) CLOSE(91825)
       RETURN
    END IF
500 CONTINUE
    IF (solve .AND. wfo) THEN
       WRITE(iout,'()')
       WRITE(iout,990) IRANK
       WRITE(iout,'()')
       IF (wff) THEN
          llev=max(2,nlev)
          IF (allzero) llev=nlev-1
          DO i=2,llev
             WRITE(iout,955) min(i,nlev),kstat(2,i-1),kstat(2,i),       &
                  real(kstat(2,i))/real(kstat(2,i-1)),kstat(1,i)
          END DO
       END IF
       WRITE(iout,'()')
       WRITE(iout,953) IRANK,flop/real(2*nza)
       !CPU_TIME: next line may be uncommented if implemented
       !          (see subroutine mestime)
       !            WRITE(iout,998) IRANK,cputm
       WRITE(iout,999) IRANK,eltm
       WRITE(iout,'()')
    END IF
    CALL sagmgpar_relmem(ijob)
    preprocessed=.FALSE.
    solve=.FALSE.
    eltm=0.0e0
    cputm=0.0e0
    flop=0.0e0
    kstat=0
    IF (wfo) THEN
       WRITE (iout,903) IRANK
    END IF
    IF (iprint.EQ.-978) CLOSE(91825)
    RETURN
901 FORMAT(i3,'*ONE APPLICATION OF AGMG PRECONDITIONER')
902 FORMAT(i3,            &
  ' (*) 1 work unit represents the cost of 1 (fine grid) residual evaluation ')
903 FORMAT(i3,'*LEAVING AGMG * (MEMORY RELEASED) ***************',&
         '***************************')
952 FORMAT(i3,'*','       Number of work units:',f9.2,             &
              ' per digit of accuracy (*)')
953 FORMAT(i3,'*',' Total number of work units:',f9.2, '(*)')
955 FORMAT('****     level',i2,'   #call=',i6,'   #cycle=',i6,    &
         '   mean=',f7.2,'    max=',i3)
990 FORMAT(i3,'*GLOBAL STATISTICS for preconditioner application:')
996 FORMAT(i3,'*','           Setup time (CPU):   ',ES10.2,     &
         ' seconds')
997 FORMAT(i3,'*','       Setup time (Elapsed):   ',ES10.2,     &
         ' seconds')
998 FORMAT(i3,'*','        Solution time (CPU):   ',ES10.2,     &
         ' seconds')
999 FORMAT(i3,'*','    Solution time (Elapsed):   ',ES10.2,     &
         ' seconds')
  END SUBROUTINE sagmgpar
!!!!!!!!!!!!!!!!!!! END of MAIN DRIVER
