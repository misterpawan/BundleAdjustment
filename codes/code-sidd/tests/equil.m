%function [R, eqA] = equil2 (A)
 
function          [R,C,eqA,flag] = equil (A)
%
% R contains the row scale factors for A.
% C contains the column scale factors for A
% amax - absolute value of largest matrix element.  If AMAX is very   
%        close to overflow or very close to underflow, the matrix   
%        should be scaled.
% R - contains the ratio of the smallest R(i) to the largest R(i).  
%         If ROWCND >= 0.1 and AMAX is neither too large nor too small, 
%         it is not worth scaling by R.
% C -  contains the ratio of the smallest   
%         C(i) to the largest C(i).  If COLCND >= 0.1, it is not   
%         worth scaling by C.
% flag = -1 : one row or one column is zero
% flag =  1 : scaling in at least one direction
% flag =  0 : no scaling

 thresh = 0.1;
  
 [m, n] = size(A);
  flag=-1 ;
         
% Compute row scale factors
  %Find the maximum element in each row
  R = max (abs(A'));
  rcmax = max (R);
  rcmin = min (R);
  
  amax = rcmax;

  if rcmin == 0.
      rcmin;
      return;
  end
  rowcnd = rcmin / rcmax;
  % Invert the scale factors
  R = 1./ R;
  eqA = A;
  nm = norm(A-A',1);
  
  if rowcnd < thresh
      % scale by R;
      RowScaling = 1;
      
    %  if nm == 0          
      %% for symmetric matrix

    %      for i = 1:n
    %          eqA(i,i:n) = R(i) * eqA(i,i:n);
    %      end
      
   %   else          
      %% for unsymmetric matrix
      
          eqA = spdiags(R',0,m,m)*eqA ;
      
    %  end
      
      flag=1;
  else
      RowScaling = 0;
  end

% Compute column scale factors, assuming the row scalings computed
  % Find the maximum element in each factor

%%
  if nm ==0
      C = max (abs(A));    % for symmetric matrix
  else
      C = max (abs(eqA));  % for unsymmetric matrix
  end
  
%%
  
  rcmax = max (C);
  rcmin = min (C);
  if rcmin == 0.
      rcmin;
      return;
  end
  C = 1./ C;

  colcnd = rcmin / rcmax;
% if colcnd >= thresh, it is not worth scaling by C
  if colcnd < thresh
      %scale by C
      ColScaling = 1;
      
      if nm == 0
      %% for symmetric matrix
          for i = 1:n
              eqA(i+1:n,i) = C(i) * eqA (i+1:n,i);
          end
          
      else
      %% for unsymmetric matrix
      
          eqA=eqA*spdiags(C',0,n,n) ;
          
      end

      flag=1 ;
  else
      ColScaling = 0;
  end
  
  if flag == -1, flag=0 ; end ;

%   rowcnd
%   colcnd
%   amax
%   
  return;
