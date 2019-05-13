function x = pre_solve(offset, L, U, T, rhs)
    y = forward_solve(offset, D, L, T, rhs);
    z = backward_solve(offset, D, U, y);
    x = [];
    for i = 1 : length(offset)-1
        x = [x z{i}];  
    end
end