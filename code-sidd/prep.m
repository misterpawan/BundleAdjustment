% Function to do the pre-processing

function [P,A_c,A,b,b_c] = prep()
    
    load JTJ.dat            % loading approximate hessian matrix csr format

    A = spconvert(JTJ);     % converting the data in matrix form
    [m,n] = size(A);        % finding the size of JTJ matrix
    b = rand(m,1);          % taking right hand side b to mx1 vector
    
    P = agtwolev(A);        % computes prolongation matrix
    A_c = (P.')*A*P;        % Galerkin Projection, coarse problem
    b_c = (P.')*b;              
end
