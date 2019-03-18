
%A is JTJ, b=ones(n,1);
load JTJ.dat            % loading data from file
A = spconvert(JTJ);     % converting the data in matrix form
[m,n] = size(A);        % finding the size of JTJ matrix
b = ones(m,1);          % taking right hand side b to mx1 vector



%solve with agmg:

x=agmg(A,b,1);         % A is symmetric positive definite
x=agmg(A,b);           % Consider A as a general matrix

x=agmg(A,b,1,1e-10,5000,1); % Verbose output
%Example of tolerance below attainable accuracy:
[x,flag,relres,iter,resvec]=agmg(A,b,1,1e-20);
%AGMG report normal convergence:

disp('Convergence flag = '),disp(flag)
% The relative residual is nevertheless larger than TOL
% (and than reported by the verbose output):
disp('Relative residual = '),disp(relres)
% But the attained accuracy is similar to the one obtained with "\":
y=A\b;
disp('Relative residual with "\" = '),disp(norm(A*y-b)/norm(b))  



