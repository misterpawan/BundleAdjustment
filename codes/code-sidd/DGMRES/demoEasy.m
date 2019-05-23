%%-------------------------------------------------
%% generate a "block " finite difference matrix -- using fd3dB
%%-------------------------------------------------
%Bx = [3 1 -1; 1 3 -1; 0.5 1 3];
%Bz = 100*ones(3,3);
% A = fd3d(8,8,3,Bx,Bx,Bz) ;  
%% shifts: try 0.1 and then 1.0 
A = fd3d(20, 20, 5, 0.1, -0.2, 0.1, 0.1);      %% <--- 3-D dim = 2000
%%A = fd3d(50, 40, 1, 0.1, -0.2, 0.1, 0.1);    %% <--- 2-D dim = 2000
%% see eigenvalues at the end with eigs(A,20,'SR')
n = size(A,1);
%%------------ params for fgmres
  im=100; maxits=200; tolIts=1.e-08;
spy(A) 
title([' block elliptic matrix of size  n = ', num2str(n)]) 
disp('pausing -- push any key to continue ') 
  pause (1.0);
close all 
%%-------------------------------------------------
%% create artificial rhs
rhs = A * ([1:n]') ; 
%% random initial guess
sol0 = rand(n,1); 
%%    S S O R 
%%
% construct SSOR preconditioner 
 L = tril(A); D = diag(A); U = triu(A);
 SSOR.L = L; 
 SSOR.D = full(D);
 SSOR.U = U;
 SSOR.kiter = 25;
%--------------------------------------------------- 
% ready to call fgmres
 disp(' solution with SSOR prec') 
 tic
 [sol,res1,its] = fgmres(A, SSOR,'precLDU', rhs, sol0,im,maxits,tolIts) ;
 toc
 disp(' ** done ** ') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   ILU 
%%
%% construct I L U (0) preconditioner 
%%%
 [L, U] = iluk(A,0) ;
 ILU = struct('L',L,'U',U) 
 [sol,res2,its] = fgmres(A, ILU,'precLU', rhs, sol0,im,maxits,tolIts) ;
 disp(' ** done ** ') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% construct I L U (1) preconditioner 
%%%
 [L, U] = iluk(A,1) ;
 ILU.L = L; ILU.U = U; 
 [sol,res3,its] = fgmres(A, ILU,'precLU', rhs, sol0,im,maxits,tolIts) ;
 disp(' ** done ** ') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-------------------------------------------------- 
% call ILU factorization from matlab -- store result
% in struct PRE 
%
 disp(' solution with ILUT-GMRES') 
 [L, U] = luinc(A, 0.001);
 ILU.L = L; ILU.U = U;
 tic
[sol,res4,its] = fgmres(A, ILU,'precLU', rhs, sol0,im,maxits,tolIts) ;
 toc 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
semilogy(res1);
semilogy(res1,'m--','marker','+','LineWidth',2,'color','k','linewidth',2)
     hold on
semilogy(res2,'g-.','marker','v','LineWidth',2,'color','g','linewidth',2)

semilogy(res3,'b-','marker','d','LineWidth',2,'color','b','linewidth',2)

semilogy(res4,'r-','marker','o','LineWidth',2,'color','b','linewidth',2)
xlabel('Iterations','fontsize',18) 
ylabel('Residual norm','fontsize',18) 
h1=legend('SSOR(4)','ILU(0)','ILU(1)','ilut')
     hold off 
set(h1,'fontsize',18)
