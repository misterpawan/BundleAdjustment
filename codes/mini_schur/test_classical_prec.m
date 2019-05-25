function test_classical_prec
    clc;
    warning('off','all');
    filepath = '../../../Test/49'
    addpath(filepath);
    addpath('../ddm');
    lhs_filename = 'JTJ49_1.mat';
    rhs_filename = 'JTe49_1.mat';

    %P = load(strcat(filepath,lhs_filename),'-mat')
    P = load(lhs_filename,'-mat')
    %size(P.A)
    B = P.lhs; %condest(B)
    [m,n] = size(B);
    %load rhs
    b = load(rhs_filename);
    b = b.rhs;
    %b = rand(m,1);


    fprintf('\nFile read complete...\n');
    
    nparts = 5;
    [sizes_parts, A, p] = domain_decomposition(B, nparts);
    sizeG = m - sum(sizes_parts(1:nparts));
    sizeD = m - sizeG; 
    
    D = A(1 : sizeD, 1 : sizeD);
    L = A(sizeD + 1: m, 1 : sizeD);
    U = L';
    G = A(sizeD + 1:m, sizeD + 1:m);
    
    [LD,UD] = lu(D);
    [UG,LG] = lu(G);
    
    tol = 1e-5; maxit = 20;restart = 20;
    tic, [x,flag,relres,iter] = gmres(B,b,restart,tol,maxit,@prec_jacobi); t_gmres = toc;
    
    its = (iter(1)-1)*restart+iter(2);
    fprintf('flag: %d\nits: %d\nrelres: %d\n', flag, its, relres);
    
    function xx = prec_jacobi(x)
        z1 = UD\LD\x(1:sizeD); 
        z2 = UG\LG\x(sizeD+1:m);
        xx = [z1;z2];
    end
    


end