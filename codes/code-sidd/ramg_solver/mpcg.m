% Two-level Grid as Preconditioner for Conjugate Gradient

 M = load('/home/siddhant.katyan/Research/SSBA-master/TEST/LadyBug/JTJ318.txt');
 A = spconvert(M);     % converting the data in matrix form
 A = tril(A.',-1) + A; % constructing full hessian matrix 
 [m,~] = size(A); 
 b = rand(m,1);        % random right hand side mx1 vector
 
 tol = 1e-6;           % tolerance
 maxit = 500;          % maximum number of iterations
 x0  = zeros(n, 1);    % Initial solution vector
 % call to matlab's pcg
 
 fprintf('doing solve with pcg...\n');
 tic
 [x,flag,relres,iter,resvec] = pcg(A, b, tol, maxit ,@pre_solve, x0);
 toc
 
 x_exact = rand(m,1);  % taking random exact solution
 b = A*x_exact;

 err_sol  = norm(x - x_exact, 2);     % 2-norm error in solution
 rel_err = err_sol / norm(x_exact);   % relative error
 
 fprintf('\n\n\n\nerror in solution = %g\nrelative error in solution=%g\n', err_sol, rel_err);
 
 fprintf('flag = %d\nrelres = %g\niters = %d\ntime=%g\n', flag, relres, iter, toc);
 