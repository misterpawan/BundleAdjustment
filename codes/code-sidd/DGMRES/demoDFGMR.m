%%-------------------------------------------------
%% generate a "block " finite difference matrix -- using fd3dB
%%-------------------------------------------------
%%Bx = [3 1 -1; 1 3 -1; 0.5 1 3];
%%Bz = 100*ones(3,3);
%%A = fd3dB(4,4,3,Bx,Bx,Bz) ;  
addpath ../precon;
addpath ../sparse;
A = fd3d(30, 40, 1, 0.2, -0.1, 0.0, 0.2) ;
n = size(A,1);
close 
spy(A) 
title([' conv. diff. matrix of size  n = ', num2str(n)]);
pause(0.1) 
%%------------ params for Dfgmres
%% NOTE KRYLOV SUBSPACE DIMENSION DIFFERENT FROM demoILU
im=50; maxits=200; tolIts=1.e-08; nev = 10;

%%
%%-------------------------------------------------
%% create artificial rhs
  rhs = A * ([1:n]') ; 
%% random initial guess
x0 = randn(n,1); 
%%-------------------------------------------------
% noPreconditioning 
% call fgmres accelerator
 PRE = struct('L',sparse(0,0),'U',sparse(0,0)) ; 
 nz1 = 0;
 [sol,res0,its0] = fgmres(A, PRE,'noPRE', rhs, x0,im,maxits,tolIts) ;
 disp(' ** done ** ') 
%--------------------------------------------------
% noPreconditioning 
% call deflated fgmres accelerator
 [sol,res1,its1] = Dfgmres(A, PRE,'noPRE', rhs, x0,im,maxits,tolIts,nev) ;
 disp(' ** done ** ') 
%--------------------------------------------------
% repeat with ILU(0) 
 disp(' solution with ILU0-GMRES') 
 [L, U] = ilu0(A);
 PRE = struct('L',L,'U',U) ;
 nz3 = nnz(L) + nnz(U) - n; 
 [sol,res2,its2] = fgmres(A, PRE,'precLU', rhs, x0,im,maxits,tolIts) ;
 disp(' ** done ** ') 
%--------------------------------------------------
% ILU0 Preconditioning 
% call deflated fgmres accelerator
 disp(' solution with ilu0-DGMRES') 
 [sol,res3,its3] = Dfgmres(A,PRE,'precLU', rhs, x0,im,maxits,tolIts,nev) ;
 disp(' ** done ** ') 
%--------------------------------------------------
close all 

semilogy([0:its0],res0,'linestyle','-','marker','*','LineWidth',3,'color','m')

hold on 

semilogy([0:its1],res1,'linestyle','--','marker','*','LineWidth',3,'color','r')

semilogy([0:its2],res2,'linestyle','-.','marker','v','LineWidth',3,'color','b') 
semilogy([0:its3],res3,'linestyle','--','marker','+','LineWidth',3,'color','k') 
xlabel('Iterations','fontsize',18) 
ylabel('Residual norm','fontsize',18) 

h1 = legend('Noprec-GMRES','Noprec-DGMRES','ILU(0)-GMRES','ILU(0)-DGMRES') 
set(h1,'fontsize',18,'location','southeast') 
