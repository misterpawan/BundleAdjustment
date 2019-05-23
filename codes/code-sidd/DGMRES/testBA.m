% Test Deflated GMRES with/without preconditioning on BA problem

fprintf('reading the matrix file...\n');
tic
M = load('/home/siddhant.katyan/Research/SSBA-master/TEST/LadyBug/JTJ372.txt');  %load matrix from specified path
toc
fprintf('matrix reading done in %d secs\n', toc);


% loading data from file

A = spconvert(M);     % converting the data in matrix form
A = tril(A.',-1) + A; % Converting full Hessian Matrix from Upper Triangular form
[m,n] = size(A);      % get the dimensions of the Matrix



%%------------ params for Dfgmres
im=1200; maxits=1000; tolIts=1.e-15; nev = 3;

% random initial guess
x_exact = rand(m,1);
% Create artificial R.H.S
b = A*x_exact;
x0 = rand(m,1);


%-----------------------------------------

% noPreconditioning 
% call fgmres accelerator
   %PRE = struct('L',sparse(0,0),'U',sparse(0,0)) ; 
  %nz1 = 0;
%  [x,res0,its0] = fgmres(A, PRE,'noPRE', b, x0,im,maxits,tolIts) ;
%  disp(' ** done ** '); 


% noPreconditioning 
% call deflated fgmres accelerator
% tic
%  [x,res,its1] = Dfgmres(A, PRE,'noPRE', b, x0,im,maxits,tolIts,nev);
%  % disp(' ** done ** ');
% toc

 
%-----------------------------------------

%With Preconditioning
%repeat with ILU(0)
% tic
%disp(' solution with ILU0-GMRES') 
 
[L, U] = ilu(A);
PRE = struct('L',L,'U',U);
nz3 = nnz(L) + nnz(U) - n;

%  [x,res2,its2] = fgmres(A, PRE,'precLU', b, x0,im,maxits,tolIts) ;
%  disp(' ** done ** ')   
%  
%  toc
%ILU Preconditioning 
 %call deflated fgmres accelerator
disp(' solution with ilu-DGMRES')
tic
 [x,res3,its3] = Dfgmres(A,PRE,'precLU', b, x0,im,maxits,tolIts,nev) ;
 disp(' ** done ** ') 
toc

% Compute Error in Solution 

err_sol  = norm(x - x_exact, 2);
rel_err = err_sol / norm(x_exact); 
res = norm(b - A*x);
rel_res = norm(b - A*x) / norm(b);

fprintf('\n\n\n\nerror in solution = %g\nrelative error in solution=%g\nresidual=%g\nrel.res=%g\ntime=%g\n', err_sol, rel_err, res, rel_res, toc);

% Results of DMGRES
%fprintf('\nResidual=%d\niters = %d\ntime=%g\n',res3, its3, toc);



