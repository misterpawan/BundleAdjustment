% test for DDM
function [T,BigBeta] = test_DDM(lhs_filename, rhs_filename, nparts, flag)
    A = load(lhs_filename, '-mat');
    disp(size(A.lhs));
    disp(nnz(A.lhs));

    B = A.lhs;
    if flag == 1
        [sizes_parts, A, p] = domain_decomposition(A.lhs, nparts);

    else
        p = load('perm_array.txt');
        Ap = A.lhs(p,p); % permute the matrix
        A = Ap;
        sizes_parts = load('sizes.txt');
    end

    disp(sizes_parts);
    sum = 1;
    offset = 1;

    for i = 1 : nparts+1
        sum = sum + sizes_parts(i);
        offset = [offset sum];
    end


    %% extract D, L, and U
    ncols = size(A,1);
    for i = 1 : nparts
        L{i} = A(offset(nparts+1) : ncols, offset(i): offset(i+1)-1);
        U{i} = A(offset(i): offset(i+1)-1, offset(nparts+1) : ncols);
        D{i} = A(offset(i): offset(i+1)-1, offset(i): offset(i+1)-1);
    end


    D{nparts+1} = A(offset(nparts+1) : ncols, offset(nparts+1): ncols);

    %% Construct the preconditioner
    [T,BigBeta] = buildPreconditioner(offset,L,D,U);

    [TL, TU, TP] = lu(T, 'vector');

    for i = 1 : length(D)
        [DL{i}, DU{i}, DP{i}] = lu(D{i}, 'vector');
    end

    %% Solve using PCG
    precon = @pre_solve;
    rhs = load(rhs_filename, '-mat');
    tol = 1e-12; maxit = 10000;

    format long;

    fprintf('Normal solve:\n ');
    d = @() normal_solve(B, rhs.rhs);
    t0 = timeit(d,1);
    x0 = normal_solve(B, rhs.rhs);
    fprintf('Time required: %f\n', t0);

    fprintf('With preconditioner:\n ');
    f = @() pcg(A, rhs.rhs(p), tol, maxit, precon);
    t1 = timeit(f,4);
    fprintf('Time required: %f\n', t1);
    [x,flag,relres,iter] = pcg(A, rhs.rhs(p), tol, maxit, precon);
    fprintf('\niterations = %d,\n relres = %g,\n flag = %d\n', iter, relres, flag);

    fprintf('Without preconditioner:\n ')
    g = @() pcg(B, rhs.rhs, tol, maxit);
    t1 = timeit(g,4);
    fprintf('Time required: %f\n', t1);
    [x2,flag2,relres2,iter2] = pcg(B, rhs.rhs, tol, maxit);
    fprintf('\niterations = %d,\n relres = %g,\n flag = %d\n', iter2, relres2, flag2);

    P = sparse(size(A,1),size(A,1));
    for i=1:length(p)
        P(i,p(i)) = 1;
    end
    x = P'*x;
    fprintf('\nDifference between normal and pcg solution: %f\n',norm(x-x0));

    function x =  pre_solve(rhs)
        y = forward_solve(offset, DU, DL, DP, L, TU, TL, TP, rhs);    
        z = backward_solve(offset, DU, DL, DP, U, y);
        x = [];
        for i = 1 : length(offset)-1
            x = [x; z{i}];  
        end
    end

end














