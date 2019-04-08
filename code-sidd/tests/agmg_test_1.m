% First testing of agmg on Bundle Adjustment Hessian Matrix


clc, clear all

%A is JTJ, b=ones(n,1);
load JTJ.dat

% loading data from file
A = spconvert(JTJ);     % converting the data in matrix form
[m,n] = size(A);        % finding the size of JTJ matrix
b = rand(m,1);          % taking right hand side b to mx1 vector
fprintf('Size: %d\n', m);


%solve with agmg:

%x=agmg(A,b,1);         % A is symmetric positive definite
%x=agmg(A,b);           % Consider A as a general matrix

%x=agmg(A,b,1,1e-10,5000,1); % Verbose output
%Example of tolerance below attainable accuracy:

tic
[x,flag,relres,iter,resvec] = agmg(A,b,1,1e-10,500);
toc

%AGMG report normal convergence:

fprintf('flag = %d\nrelres = %g\niters = %d\ntime=%g\n', flag, relres, iter, toc);

%disp('Convergence flag = '),disp(flag)
% The relative residual is nevertheless larger than TOL
% (and than reported by the verbose output):
%disp('Relative residual = '),disp(relres)
% But the attained accuracy is similar to the one obtained with "\":

