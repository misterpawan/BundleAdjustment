

clc, clear all

%A is JTJ, b=ones(n,1);
fprintf('reading the matrix file...\n');
tic
M = load('/home/siddhant.katyan/Research/SSBA-master/TEST/LadyBug/JTJ49.txt');
toc
fprintf('matrix reading done in %d secs\n', toc);

% loading data from file
A = spconvert(M);     % converting the data in matrix form
A = tril(A.',-1) + A;
[m,n] = size(A);        % finding the size of JTJ matrix

fprintf('cond no of A = %g\n', condest(A));
A = -A;
fprintf('doing equilibration\n');
[~, ~, A, ~] = equil(A);
fprintf('cond no of equilibrated A = %g\n', condest(A));


%[P,R,C] = equilibrate(A);
%A = R*P*A*C;

%dd = 1000000; A = A + dd*speye(m,m);
%b = rand(m,1);          % taking right hand side b to mx1 vector
fprintf('Size: %d\n', m);

% symm check
%issymm2 = norm(A - A', inf);

mysymm = issymmetric(A);
IsDiagDom = @(A) all( 2 * abs(diag(A)) > sum(abs(A),2) );

x_exact = rand(m,1);
b = A*x_exact;

%solve with agmg:

%x=agmg(A,b,1);         % A is symmetric positive definite
%x=agmg(A,b);           % Consider A as a general matrix

%x=agmg(A,b,1,1e-20,500,1); % Verbose output
%Example of tolerance below attainable accuracy:

fprintf('doing solve with agmg...\n');
tic
[x,flag,relres,iter,resvec] = agmg(A,b,0,1e-14,500,1);
toc

err_sol  = norm(x - x_exact, 2);
rel_err = err_sol / norm(x_exact); 
fprintf('\n\n\n\nerror in solution = %g\nrelative error in solution=%g\n', err_sol, rel_err);

%AGMG report normal convergence:
if (IsDiagDom(A) == 0)
	fprintf('Warning: matrix is not diag dominant!!!\n')
end

if (mysymm == 0)
   symmnrm = norm(A - A', inf);
   fprintf('norm A-A = %d\nWarning: the matrix is not symmetric!!!!!\n', symmnrm);
end
fprintf('flag = %d\nrelres = %g\niters = %d\ntime=%g\n', flag, relres, iter, toc);

%disp('Convergence flag = '),disp(flag)
% The relative residual is nevertheless larger than TOL
% (and than reported by the verbose output):
%disp('Relative residual = '),disp(relres)
% But the attained accuracy is similar to the one obtained with "\":

