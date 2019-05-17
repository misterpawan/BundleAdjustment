% test non-zeros
P = load('JTJ49.txt'); 
B = spconvert(P);
B = B + triu(B,1)'; % since matrix is stored in symm form, need to copy to lower

nparts = 3;
[sizes_parts, A, p] = domain_decomposition(B, nparts);
[m,n] = size(B); 
%b = rand(m,1);

sizeG = m - sum(sizes_parts(1:nparts));
sizeD = m - sizeG; 

D = A(1 : sizeD, 1 : sizeD);
L = A(sizeD + 1: m, 1 : sizeD);
U = L';
G = A(sizeD + 1:m, sizeD + 1:m);

nnz(A)

% compute Schur compl
S = G - L * (D \ U);


nnz(S)

[LL, UU] = lu(A);

nnz(LL) + nnz(UU)