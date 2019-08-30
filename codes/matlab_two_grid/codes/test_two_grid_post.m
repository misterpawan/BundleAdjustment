% test two grid for bundle adjustment problem

function test_two_grid
clc, clear

%filename = 'JTJ49';
filename = 'JTJ372_1';
%filename = 'JTJ646';


addpath(genpath('/home/iiit/siddhant.katyan/test1/LadyBug/')); % from IIIT

A = load(strcat(filename, '.txt'));
A = spconvert(A);
%A = Problem.A; 
A = A + triu(A,1)';
nc = '100';
[m,n] = size(A)

fileID = fopen('/home/siddhant.katyan/test1/LadyBug/JTe372_1.txt','r');
formatSpec = '%f';
b = fscanf(fileID,formatSpec);
% 
% [P2,R2,C2] = equilibrate(A);
% A = R2*P2*A*C2;

%matrix_file = strcat(filename, '.mtx');
%matrix_file = strcat('/home/pawan/work/datasets/', matrix_file);

%[status,cmdout] = 

%system(sprintf('./print_nd_perm2 %s %s', matrix_file, nc));

%pause(2)

%system(['kill ' cmdout]);

 % nc = str2num(nc);
%keyboard;

% P = sparse(n,nc);
% parts = load('parts.txt');
% parts = parts + 1;
 
 %P(1:m, parts(1:m)) = 1;
 
%  for i = 1:m
%      P(i, parts(i)) = 1;
%  end
 
% spy(P)
 %keyboard;
 
% rank(full(P))

nc2 = 5;

% tic
% [Vs, Ds] = eigs(A, nc2, 'smallestabs');
% t_smallest = toc
% 
% Ds

tic
[Vl, Dl] = eigs(A, nc2);
t_largest = toc

%Dl

P = Vl;

%P = [P V];
%P = Vl(:,2);

%P = rand(m,100);

%P = zeros(m,1); P(1) = 1;

tic
Ac = P' * A * P;
toc

% 
% size(Ac)
% 
% condest(Ac)
% rank(full(Ac))

addpath('/home/siddhant.katyan/BundleAdjustment/codes/ddm');

nparts = 10;

[sizes_parts, B, p] = domain_decomposition(A, nparts);

nJ1 = sum(sizes_parts(1:nparts));

nJ2 = n - nJ1;

J = blkdiag(A(1:nJ1, 1:nJ1),  A(nJ1+1 : n, nJ1+1: n));


% sz = 369;
% nJ = 614
% J = [];
% r1 = 1; 
% for i = 1 : nJ   
%     r2 = r1 + sz - 1;
%     T = A(r1:r2, r1:r2);
%     J = blkdiag(J,T);
%     r1 = r1 + sz;
% end


[L, U] = lu(J);


%keyboard;


% 
% % define smoother
% setup.type = 'ilutp';
% setup.droptol = 1e-7;
% setup.udiag = 1;
% [L, U] = ilu(A, setup);
 
%  nJacobi = 10;
%  r = rem(m, nJacobi);
 



%b = rand(m,1); 
tol = 1e-2; maxit = 100; 
im = 10;

x_known = rand(m,1);
%b = A*x_known;

norm(x_known)

tic 

[x, flag, relres, its, res] = gmres(A, b, im, tol, maxit, @tg_solve);

%[x, flag, relres, its, res] = gmres(A, b, im, tol, maxit);

%[x, flag, relres, its, res] = bicgstab(A, b, tol, maxit, @tg_solve);

%[x, flag, relres, its, res] = pcg(A, b, tol, maxit, @tg_solve);


%[x, flag, relres, its, res] = pcg(A, b, tol, maxit);


%[x, flag, relres, its, res] = minres(A, b, tol, maxit);

t_pcg = toc;


norm(b - A*x)

norm(b)
norm(x)

norm(b - A*x) / norm(b)

%norm(x - xg)
norm(x - x_known)


iter = (its(1)-1)*im + its(2)
its
%its
relres
%res
flag
% 
% fprintf('its = %d\nrelres = %g\nflag=%gtime = %g\n', its, relres, flag, t_pcg);
% 
function y = tg_solve(x)

    t = U \ (L \ x);           % Pre-Smoothing
    r = x - A*t;
    g = P * (Ac \ (P' * r));
    q = A*g;
    r = r - q;
    g = g + U \ (L \ r);
    y = t+g;
    
% %   %  y =  U \ (L \ x);
%    
%   % y  = P * (Ac \ (P' * x));
%   
%   t = U \ (L \ x);
%     
%   g = P * (Ac \ (P' * x));
%   q = A*g;
%   
%   y = t + g - U \ (L \ q);

end


end

