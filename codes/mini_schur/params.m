%---------------------------------------------------- 
% This file  is read by  many of the codes.  It sets
% all  the  parameters  needed  for setting  up  the
% preconditioner and for  the iterations (global and
% inner).  This avoids  having to  search  for these
% parameters  in the  various files  to  change them
% when doing tests.
%----------------------------------------------------
%%
%% PART 1: parameters for the arms reduction
%% 
%%
tolmac   = 1.e-16  ;    % machine eps [used only in gmres.] can be removed
maxits   = 50;
tol      = 1.e-10;
outputG  = 0;
im       = 50;
                         % have the meaning of level of fill when 
                         % level strategies will be used later. 

%-----------------------------------------------------------------------



