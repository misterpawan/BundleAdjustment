%This file tests the domain_decomposition.m function
function test_BA
    clc;
    filepath = '../../../Test'
    addpath(filepath);
    addpath('../ddm');
    lhs_filename = 'JTJ49.mat';
    rhs_filename = 'b49.mat';

    %P = load(strcat(filepath,lhs_filename),'-mat')
    P = load(lhs_filename)
    %size(P.A)
    B = P.lhs;
    nparts = 5;
    %flag = 1; %1 for DDM
    %spy(A);
    %[T,BigBeta] = test_DDM(lhs_filename,rhs_filename,nparts,flag);

    [sizes_parts, A, p] = domain_decomposition(P.lhs, nparts);

    %sizes_parts
     disp(sizes_parts);
        sum = 1;
        offset = 1;

        for i = 1 : nparts+1
            sum = sum + sizes_parts(i);
            offset = [offset sum];
        end

        %offset

        %% extract D,F,E,G
        %%      | D   |  U |   
        %%  A = |----------|
        %%      | L   |  G |
        %%
        ncols = size(A,1);
        for i = 1 : nparts
            L{i} = A(offset(nparts+1) : ncols, offset(i): offset(i+1)-1);
            U{i} = A(offset(i): offset(i+1)-1, offset(nparts+1) : ncols);
            D{i} = A(offset(i): offset(i+1)-1, offset(i): offset(i+1)-1);
        end
        size_L = size(L{1})
        size_U = size(U{1})
        size_D = size(D{1})

        G = A(offset(nparts+1) : ncols, offset(nparts+1): ncols);
        D{nparts+1} = G;
        size_G = size(D{nparts+1})

         %% Construct the preconditioner
        [T,BigBeta] = buildPreconditioner(offset,L,D,U);

        [TL, TU, TP] = lu(T, 'vector');

        for i = 1 : length(D)
            [DL{i}, DU{i}, DP{i}] = lu(D{i}, 'vector');
        end

         %% Solve using PCG
        precon = @pre_solve;
        rhs = load(rhs_filename);
        tol = 1e-12; maxit = 10000;

        format long;

        fprintf('Normal solve:\n ');
        d = @() normal_solve(B, rhs.rhs);
        t0 = timeit(d,1);
        x0 = normal_solve(B, rhs.rhs);
        fprintf('Time required: %f\n', t0);

        fprintf('With preconditioner:\n ');
        %f = @() pcg(A, rhs.rhs(p), tol, maxit, precon);
        f = @() pcg(A, rhs.rhs, tol, maxit, precon);
        t1 = timeit(f,4);
        fprintf('Time required: %f\n', t1);
        %[x,flag,relres,iter] = pcg(A, rhs.rhs(p), tol, maxit, precon);
        [x,flag,relres,iter] = pcg(A, rhs.rhs(p), tol, maxit, precon);
        fprintf('\niterations = %d,\n relres = %g,\n flag = %d\n', iter, relres, flag);

        fprintf('Without preconditioner:\n ')
        g = @() pcg(B, rhs.rhs, tol, maxit);
        t1 = timeit(g,4);
        fprintf('Time required: %f\n', t1);
        [x2,flag2,relres2,iter2] = pcg(B, rhs.rhs, tol, maxit);
        fprintf('\niterations = %d,\n relres = %g,\n flag = %d\n', iter2, relres2, flag2);

        R = sparse(size(A,1),size(A,1));
        for i=1:length(p)
            R(i,p(i)) = 1;
        end
        x = R'*x;
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
    
    