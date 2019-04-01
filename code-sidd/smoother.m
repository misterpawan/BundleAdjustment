% Function to compute the smoothing matrix 
% Smoother can be gauss-seidel, SSOR to reduce the "rough" error

function [M_amg] = smoother(A,A_c)
   
    % Pre-smoothing fine grid
    setup.milu = 'row';          % Choosing two setup options
    setup.type = 'crout';

    
    [L,U] = ilu(A,setup);   % incomplete of cholesky of sparse matrix A
    M_s = inv(L*U);         % smoother matrix reduces "rough error"
    M_c = inv(A_c);         % Coarsening matrix = P* (A)^-1 * P^T 
    M_amg = M_s + M_c - (M_s*A*M_c);   % smoothed matrix 
   
end
