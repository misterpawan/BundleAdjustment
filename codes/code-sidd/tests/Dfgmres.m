 function [sol,res,its]=Dfgmres(A,PRE,precfun,rhs,sol,im,maxits,tolIts,NEV) 
%-----------------------------------------------------------------------
%% deflated gmres. 
% restarted gmres with Krylov subspace of dim = im.  
% NOTE: this is actually fllexibe (FGMRES) -- allows
% variations in preconditioner 
%
%-----------------------------------------------------------------------
%% 
 outputG = 1;
 tolmac = 2.2e-16;
 n = size(A,1)    ;
 its = 0    ;
 Y = [];
% 
% main loop 
%
  while (its < maxits) 
      vv(:,1) = rhs - A*sol  ;
      ro = norm(vv(:,1),2)  ; 
      res(its+1) = ro;
      if (its  == 0) 
	tol1=tolIts*ro  ;
      end  ; 
      if (ro <= tol1 | its >= maxits)  
	return
      end
      t = 1.0/ ro   ;
      vv(1:n,1) = vv(1:n,1) * t  ;
%------------initialize 1-st term  of rhs of hessenberg system
      rs(1) = ro  ;
%%------------ print its/residual info
      if (outputG) 
         fprintf(1,' its %d  res %e \n',its,ro)
      end 
%%------------ arnoldi loop
   i = 0  ;
   while (i < im  &  (ro  >  tol1)  &  its < maxits)
      i=i+1  ;
      its = its + 1  ;
      i1 = i + 1 ; 
%%------------ if at end use approx e-vector 
      if ((i>im-NEV)  & size(Y,2)>0) 
          z = Y(:,i+NEV-im);
%%------------ else use preconditioned column of V
      else 
          z = feval(precfun,PRE,vv(:,i));
      end
      Z(:,i) = z; 
%%------------  modified GS  ;
      vv(1:n,i1) = A*z; 
      for j=1:i
           t = vv(1:n,j)'*vv(1:n,i1)  ;
           hh(j,i) = t  ;
	   vv(1:n,i1) = vv(1:n,i1) - t*vv(1:n,j)  ;
      end  ;
%%------------normalize 
      t = norm(vv(1:n,i1),2)  ;
      hh(i1,i) = t  ;
      if (t  ~= 0.0)  
        t = 1.0 / t  ;
        vv(1:n,i1) = vv(1:n,i1)*t  ;
      end  ; %% IF 
%% 
    if (i ~= 1) 
%------------apply previous rot.s on i-th column of H
       for k=2:i 
          k1 = k-1  ;
          t = hh(k1,i)  ;
          hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)  ;
          hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)  ;
       end;  %% FOR 
    end  ; %% IF 
%%
      gam = sqrt(hh(i,i)^2 + hh(i1,i)^2)  ;
      if (gam  == 0.0) 
	gam = tolmac  ;
      end  ; 
%------------determine plane rotation + update rhs of ls pb
      c(i) = hh(i,i)/gam  ;
      s(i) = hh(i1,i)/gam  ;
      rs(i1) = -s(i)*rs(i)  ;
      rs(i) =  c(i)*rs(i)  ;
%------------ determine res. norm. and test convergence
      hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)  ;
      ro = abs(rs(i1))  ;
      res(its+1) = ro   ; 
      if (outputG) 
         fprintf(1,' its %d  res %e \n',its,ro)
      end 
   end   ;  %% end of while (im) arnoldi loop 
%------------compute solnu -1st solve upper triangular system
      rs(i) = rs(i)/hh(i,i)  ;
      for  k=i-1:-1:1 
         t=rs(k);
         for j=k+1:i  
            t = t-hh(k,j)*rs(j)  ;
         end  
         rs(k) = t/hh(k,k)  ;
      end  
%------------done with back substitution..  
%            now form linear combination to get soln
      for j=1:i 
         sol = sol+rs(j)*Z(:,j)   ;
      end  ;
      if ((ro  <=  tol1) | (its >= maxits)) return; end; 
%%------------ compute eigenvectors for deflation
      R = triu(hh(1:im,1:im));
      S = vv'*Z;   
%%------------ apply rotations to S
    for k=2:i
        k1 = k-1  ;
        t = S(k1,:)  ;
        S(k1,:) = c(k1)*t + s(k1)*S(k,:)  ;
        S(k,:) = -s(k1)*t + c(k1)*S(k,:)  ;
    end 
%%------------ solve generalized ev Pb/ 
    [Q, D] = eig(R,S(1:im,1:im));
    d = diag(D);
    [dummy,indx] = sort(abs(d),'ascend') ;
    Q = cmplx2real(d, Q, indx, NEV); 
    Y = Z*Q(:,1:NEV);
end;  %% end while -- restart 
