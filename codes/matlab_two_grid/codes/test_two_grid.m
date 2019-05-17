% test two grid for bundle adjustment problem

addpath ''

Problem = load('JTJ49.mat');
A = Problem.A;

%keyboard;

Problem = load('Parts.txt');

P = Problem.P;

Ac = P' * A * P;
% define smoother
setup.type = 'nofill';
[L, U] = ilu(A, setup);

tic 
[x, its, relres, res] = pcg(A, b, tol, maxit, @tg_solve);
t_pcg = toc;

fprintf('its = %d\nrelres = %g\nflag=%gtime = %g\n', its, relres, flag, t_pcg);

function y = tg_solve(x)

   t = U \ (L \ x);
   g = P * (Ac \ (P' * x));
   q = A*g;
   y = t + g - U \ (L \ q);

end

