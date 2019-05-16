% Multigrid as a Preconditioner

function [y] = pre_solve(x)
    
    M = load('/home/siddhant.katyan/Research/SSBA-master/TEST/LadyBug/JTJ318.txt');
    A = spconvert(M);     % converting the data in matrix form
    A = tril(A.',-1) + A; % constructing full hessian matrix 
    [m,~] = size(A); 
    b = rand(m,1);       % random right hand side mx1 vector
    
    % Call to two-level multigrid method
    [y, ~, ~] = ramg(A, b, x, @GS_Iter, 1, 1, 1e-10, 20);
    
end
