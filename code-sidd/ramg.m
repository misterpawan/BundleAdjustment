% Algebraic Recursive Multilevel Solver

function [y,res_err] = ramg(A)

    [P,A_c,A,b,b_c] = prep();     % call to Preprocessing function
    
    M_amg = smoother(A,A_c);      % Pre-smoothing fine grid
    z = pcg(M_amg,b);             % using preconditioned conjugate gradient on smoother data   
    
    % Recursive call to Coarse grid solver missing, what should be the base
    % case??
    
    t = cogs(A_c,P,b_c);        % coarse-grid correction 
    t = A*t;                      % Interpolation  
    
    y = z + pcg(M_amg,t);               % Post-smoothing fine grid 

    res_err = norm(A*y-b)/norm(b);    % relative residual error
    
 end   
    
    
