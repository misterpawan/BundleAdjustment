% Function to generate coarse grid Matrix 

% A  : [in]  Fine grid problem matrix - J^TJ (Hessian Matrix)
% CA : [out] Coarse grid problem coefficient matrix
% P  : [out] Interpolation operator, CA = P^T * A * P (Galerkin Projection)

function [CA, P] = cogs(A)
    
    % Computes aggregation and associated prolongation matrix(P)
    [P, Ind] = agtwolev(A);
    
    % Construct coarse grid matrix - Galerkin Projection
    CA = P' * A * P;
    
end 