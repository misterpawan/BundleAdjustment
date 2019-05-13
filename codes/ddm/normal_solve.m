function  x = normal_solve(A, B)
    B = sparse(B);
    x = A \B;
end