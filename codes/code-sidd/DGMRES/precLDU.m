 function x = precLDU(PRE, rhs) 
%% function rhs = precLDU(PRE, rhs) 
%% arms preconditioning operation
%% PRE = struct for preconditioner
%%-------------------------------------------------

L = PRE.L;
D = PRE.D;
U = PRE.U; 
kiter = PRE.kiter;

x = U \ ( D .* (L \ rhs)) ;

for k = 2:kiter 
  x = x +  U \( D .* (L \ ( rhs - U*x - L*x +  D .*x ) ));
end
