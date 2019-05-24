 function [L,U] = ilu0(A) 
%---------------------------------------------
% function [L,U] = ilu0 (A) 
% ILU factorization of A. Uses ikj variant of GE
%---------------------------------------------
n = size(A,1) ;
for i=1:n 
    [ii,jj,rr] = find(A(i,:));
    nzr = length(jj);
    p = 1; 
    k = jj(p) ;
    while (k < i)    
        A(i,k) = A(i,k)/A(k,k); 
        piv = A(i,k);
       for j=p+1:nzr 
          A(i,jj(j)) = A(i,jj(j)) - piv*A(k,jj(j));
       end
	p = p+1;
        k = jj(p) ;
     end  %% while
 end 
 L = tril(A,-1) + spdiags(ones(n,1),0,n,n) ;
 U = triu(A);
