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
tolInd   = 0.1;         % tolerance for diagonal dominance filtration.
ilutolB  = 0.00 ;         % tolerance for dropping in B
bsize    = 50;           % block size in reduction
lfil     = 1000;         % lfil in L, U 
droptolE = 0.000  ;      % drop tol for E*inv(U) and inv(L)*F
lfilE    = 1000;         % lfil for E*inv(U) and inv(L)*F  may 
                         % have the meaning of level of fill when 
                         % level strategies will be used later. 
%%
%% Part 2 : parameters for the last schur complement 
%%
droptolS = 0.00 ;       % droptol used to prune S
ilutolS  = 0.00  ;       % droptol for S = LS  US 
lfilS    = 1000;        % fill paramerer for S = LS US 
imS      = 0;           % krylov dimension and its # for last level
                        % set to zero if no iterations wanted
                        % at last level..
tolInner = 0.00 ;       % tolerance for last-level GMRES (pgmres). 
outputI  = 0;           % whether or not to print residual vs iterations
                        % in inner gmres (pgmres)
%% 
%% Part 3 : parameters for the outer gmres iteration
%% 
im       = 40     ;     % krylov subspace dimension 
tolIts   = 0.00000000001 ;      % stopping criterion
maxits   = 100     ;    % max number of outer steos      --
outputG  = 1       ;    % whether or not to print residual vs iterations
                        % in global gmres. 
%% 
%% Part 4 : general paramters 
%%
tolmac   = 1.e-16  ;    % machine eps [used only in gmres.] can be removed
%
%-----------------------------------------------------------------------



