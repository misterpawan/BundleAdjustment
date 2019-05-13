% test forward and backard solve
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

ncols = size(A,1);
for i = 1 : nparts
    L{i} = A(offset(nparts+1) : ncols, offset(i): offset(i+1)-1);
    U{i} = A(offset(i): offset(i+1)-1, offset(nparts+1) : ncols);
    D{i} = A(offset(i): offset(i+1)-1, offset(i): offset(i+1)-1);
end
D{nparts+1} = A(offset(nparts+1) : ncols, offset(nparts+1): ncols);

T = rand(721,721);
rhs = rand(4096,1);

y = forward_solve(offset, D, L, T, rhs); 
yy = [];
for i=1:nparts+1
    yy = [yy; y{i}];
end

t1 = A11 \ rhs(1:3375);
t2 = T \ (rhs(3376:4096) - A21*t1);

tt = [t1;t2];

norm(yy(3376:4096) - tt(3376:4096))
