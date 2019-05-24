  function [x, its] = pcg (A, PRE, precfun, rhs, x0, m, tol) 
%% solves A x = rhs by pcg 
%%------------------------------------------ 
 n = size(A,1); 
 x = x0; 
 r = rhs - A * x;
 z = feval(precfun,PRE,r);
 p = z ;
 ro1 = z' * r; 
 tol1 = tol*tol*ro1; 
%%
 its = 0 ;
 while (its < m && ro1 > tol1) 
	 its = its+1; 
	 ro = ro1; 
	 ap = A * p; 
	 alp = ro / ( ap'* p ) ;
	 x = x + alp * p ;
	 r = r - alp * ap; 
         z = feval(precfun,PRE,r);
         ro1= z'*r;
	 bet = ro1 / ro ;
	 p = z + bet * p; 
  end
 end 
         
