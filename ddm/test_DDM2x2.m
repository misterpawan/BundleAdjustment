% test DDM 2x2 case
function [Beta1,T] = test_DDM2x2
% clc
%clear all;
p = load('perm_array.txt');
A = mmread('Sky3d16.mtx');
Ap = A(p,p); % permute the matrix
A = Ap; ncols = size(A,1);

nparts = 8;
sizes_parts = load('sizes.txt');
sum = 1;
offset = 1;

for i = 1 : nparts+1
    sum = sum + sizes_parts(i);
    offset = [offset sum];
end

A11 = Ap(1:offset(nparts+1)-1, 1:offset(nparts+1)-1);
A12 = Ap(1:3375, 3376:4096);
A21 = Ap(3376:4096, 1:3375);
A22 = Ap(3376:4096, 3376:4096);
UtO = A12*ones(size(A12,2),1);
tmp = A11 \ UtO; %keyboard;
indNotZeros = (abs(UtO) > 1e-14);
Beta1 = zeros(3375,1);


Beta1(indNotZeros) = tmp(indNotZeros) ./ UtO(indNotZeros);
T = A22 - A21 * diag(Beta1) * A12;

precon = @pre_solve;
% rhs = ones(ncols,1);
rhs = load('rhs.txt');

tol = 1e-6; maxit = 100;

[x,flag,relres,iter] = pcg(Ap, rhs, tol, maxit, precon);

fprintf('\niterations = %d,\n relres = %g,\n flag = %d\n', iter, relres, flag);

function y = pre_solve(x)
    %keyboard;
    t1 = A11 \ x(1:3375);
    t2 = T \ (x(3376:4096) - A21*t1);
    
    y2 = t2;
    y1 = t1 - A11 \ (A12 * y2);
    
    y = [y1;y2];

end

end