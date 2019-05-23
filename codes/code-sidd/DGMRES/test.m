addpath ../sparse;

A = fd3d(30, 40, 1, 0.2, -0.3, 0.0, 0.45);
opts.thresh   = 0.0;
opts.droptol  = 0.05;
[L, U] = luinc(A, opts);
x = ones(1200,1);
b = L*(U*x);

norm(x - U \ (L\b)) 

norm(L\x) 
norm(U\x);
 norm(U\(L\b),'inf') 


z = U \ (L\b) ;
plot(z) 


%norm(full( inv(L*U)*b - U\(L\b) ))
%AA = L * U;
%norm(full( inv(AA)*b - AA\b ))

