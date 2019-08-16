%This file tests the domain_decomposition.m function
function test_BA
    clc;
    warning('off','all');
    filepath = '~/Test/138';
    addpath(filepath);
    addpath('../ddm');
    lhs_filename = 'JTJ138_1.mat';
    rhs_filename = 'JTe138_1.mat';

    %P = load(strcat(filepath,lhs_filename),'-mat')
    P = load(lhs_filename,'-mat')
    %size(P.A)
    B = P.lhs; %condest(B)
    [m,n] = size(B);
%     
%     rand_rhs = ones(m,1);
%     fp = fopen("~/rand_rhs_A.txt","w");
%     fprintf(fp,"%12.12f\n",rand_rhs);
%     fclose(fp);
%     
%     [LB,UB] = lu(B);
%     rand_sol = UB\(LB\rand_rhs);
%     fp = fopen("~/rand_sol_A.txt","w");
%     fprintf(fp,"%12.12f\n",rand_sol);
%     fclose(fp);
%     keyboard;
    
    %load rhs
    b = load(rhs_filename);
    b = b.rhs;
    %b = rand(m,1);
     %C = full(B(1:10,1));
%      for  i = 1:10
%          fprintf("\nb[%d] = %15.15f\n",i,b(i,1));
%      end
%         keyboard;
%      LDL_sol = load('~/Test/49/JTe49_1_Perm.txt');
    fprintf('\nFile read complete...\n');
    

    tic;
    xx = normal_solve(B, b);
%     [LB,DB] = ldl(full(B));
%     xx = LB'\(LD\(LB\b));
    normal_solve_time = toc;
    fprintf("\nDirect solve done!\n");
%     xx = full(xx);
%     for i = 1 : 10
%         fprintf("\nxx[%d] = %f\n",i,xx(i));
%     end
%      keyboard;
    
    nparts = 10;
    %[sizes_parts, A, p] = domain_decomposition(B, nparts);
    A=B;
    %sizeG = m - sum(sizes_parts(1:nparts))
    sizeG = 1242;
    sizeD = m - sizeG; 
    %keyboard;
%     D = A(1 : sizeD, 1 : sizeD);
%     L = A(sizeD + 1: m, 1 : sizeD);
%     U = L';
%     G = A(sizeG + 1:m, sizeG + 1:m);

    %% Identify MSCs
    % blocks for schur complements
    %nmsc = 3; 
   % nmsc_blocks = [3 4 5 6 7 8 9 10 15 20];
     nmsc_blocks = [10 20 30];
    %nmsc_blocks = [30];
    %nmsc_blocks = [25 30 35 40];
    %nmsc_blocks = [13 14 15 16 17 18];
    iters = zeros(1,length(nmsc_blocks));
    gmres_time = zeros(1,length(nmsc_blocks));
    rel_res = zeros(1,length(nmsc_blocks));
    total_time = zeros(1,length(nmsc_blocks));
    conv_flag = zeros(1,length(nmsc_blocks));
    
    iters_np = zeros(1,length(nmsc_blocks));
    gmres_time_np = zeros(1,length(nmsc_blocks));
    rel_res_np = zeros(1,length(nmsc_blocks));
    total_time_np = zeros(1,length(nmsc_blocks));
    conv_flag_np = zeros(1,length(nmsc_blocks));
    
    
    
    
    for z = 1:length(nmsc_blocks)
        D = A(1 : sizeD, 1 : sizeD);
        L = A(sizeD + 1: m, 1 : sizeD);
        U = L';
        G = A(sizeD + 1:m, sizeD + 1:m);%%change
        
%         Schur = G - L*(D\U);
        
%         [L_block_G,U_block_G] = lu(G);
        %keyboard;
        r = rem(sizeG, nmsc_blocks(z)); sz = (sizeG - r)/nmsc_blocks(z); r1 = 1; r2 = sz;
        GS = []; %ol = 0;

      
        %nmsc
        %sz

        % blocks for D: good to take same number of blocks for D as for G
        nD = nmsc_blocks(z); rD = rem(sizeD, nD); szD = (sizeD - rD)/nD; rD1 = 1; rD2 = szD;
        %sz is the size of each block
        %nD
        %szD

        % overlap size
        %ol = floor(sz);
        ol = 0;%sz; % try with zero overlap

        %keyboard;

        fprintf('Computing Mini Schur complements ...\n');
        for i = 1 : nmsc_blocks(z)-2
           clear PD PU PL PG
           PD = D(rD1 : rD2 + ol, rD1 : rD2 + ol); 
           PU = U(rD1 : rD2 + ol, r1 : r2 + ol); 
           PL = PU'; % for BA L = U'
           PG = G(r1 : r2 + ol, r1 : r2 + ol);

          % printf block sizes
%            fprintf('\nPrinting block sizes for mini Schur complement: %d\n', i);
%            fprintf('block size of G = %d x %d (%d : %d, %d : %d) \n', size(PG,1), size(PG,2), r1,r2, r1, r2);
%            fprintf('block size of U = %d x %d (%d : %d, %d : %d) \n', size(PU,1), size(PU,2), rD1, rD2, r1, r2);
%            fprintf('block size of L = %d x %d (%d : %d, %d : %d) \n', size(PL,1), size(PL,2), r1, r2, rD1, rD2);
%            fprintf('block size of D = %d x %d (%d : %d, %d : %d) \n\n', size(PD,1), size(PD,2), rD1, rD2, rD1, rD2);

           S{i} = PG - PL * (PD \ PU); 
           GS(r1:r2 + ol, r1 : r2 + ol) = S{i}; 
           GS = sparse(GS); 
           r1 = r2 + 1; 
           r2 = r2 + sz;
           rD1 = rD2 + 1;
           rD2 = rD2 + szD; 
           %nnz(GS)
        end
        %keyboard;
        oll = 0; % overlap size for 2nd last msc
        i = i+1;

        if (i == nmsc_blocks(z)-1)
           clear PD PU PL PG
           PD = A(rD1:rD2+oll, rD1:rD2+oll); 
           PU = U(rD1:rD2+oll, r1:r2+oll);
           %PL = L(r1:r2+oll, r1:r2+oll); 
           PL = PU';
           PG = G(r1:r2+oll, r1:r2+oll);

          % printf block sizes
%            fprintf('Printing block sizes for mini Schur complement: %d\n', i);
%            fprintf('block size of G = %d x %d (%d : %d, %d : %d) \n', size(PG,1), size(PG,2), r1,r2, r1, r2);
%            fprintf('block size of U = %d x %d (%d : %d, %d : %d) \n', size(PU,1), size(PU,2), rD1, rD2, r1, r2);
%            fprintf('block size of L = %d x %d (%d : %d, %d : %d) \n', size(PL,1), size(PL,2), r1, r2, rD1, rD2);
%            fprintf('block size of D = %d x %d (%d : %d, %d : %d) \n\n', size(PD,1), size(PD,2), rD1, rD2, rD1, rD2);

           S{i} = PG - PL * (PD \ PU); 
           GS(r1:r2+oll, r1:r2+oll) = S{i};
           GS = sparse(GS); 
           r1 = r2+1;
           rD1 = rD2+1;
        end

        r2 = sizeG; i = i + 1;
        rD2 = sizeD; 

        if (i == nmsc_blocks(z))
          clear PD PU PL PG
          PD = A(rD1:rD2, rD1:rD2); 
          PU = U(rD1:rD2, r1:r2);
          %PL = L(r1:r2, r1:JTJ138r2); 
          PL = PU';
          PG = G(r1:r2, r1:r2);

           % printf block sizes
%            fprintf('Printing block sizes for mini Schur complement: %d\n', i);
%            fprintf('block size of G = %d x %d (%d : %d, %d : %d) \n', size(PG,1), size(PG,2), r1,r2, r1, r2);
%            fprintf('block size of U = %d x %d (%d : %d, %d : %d) \n', size(PU,1), size(PU,2), rD1, rD2, r1, r2);
%            fprintf('block size of L = %d x %d (%d : %d, %d : %d) \n', size(PL,1), size(PL,2), r1, r2, rD1, rD2);
%            fprintf('block size of D = %d x %d (%d : %d, %d : %d) \n\n', size(PD,1), size(PD,2), rD1, rD2, rD1, rD2);


          S{i} = PG - PL * (PD \ PU); 
          GS(r1:r2, r1:r2) = S{i};
          GS = sparse(GS); 
        end 
        %det_GS = det(GS)
        fprintf('done with mini Schur complements!\n'); 
        fprintf('Computing LU of MSC ...');
        tic,[LG,UG] = lu(GS);t_factor_msc = toc; fprintf('done!!!\n');
%          fprintf('Computing Cholesky of MSC ...');
%          tic,UG = chol(GS);t_factor_msc = toc; LG = UG';
         fprintf('done!!!\n'); % matlab chol factorizes into upper triangular 
        %keyboard;
%         rand_rhs = rand(sizeG,1);
%         fp = fopen("~/rand_rhs_MSC_20.txt","w");
%         fprintf(fp,"%10.9f\n",rand_rhs);
%         fclose(fp);
        
%         rand_sol = UG\(LG\rand_rhs);
%         fp = fopen("~/rand_sol_MSC_20.txt","w");
%         fprintf(fp,"%10.9f\n",rand_sol);
%         fclose(fp);
%           GS = full(GS);
%           for i = 62:124
%               fprintf("GS(%d) = %f\n",i,GS(99,i));
%           end
          
%            keyboard;
        clear S PD PU PL PG
        clear GS
 
%         setup.type='ilutp'; 
%         setup.droptol = 1e-03; %ntol(ii); 
%         setup.udiag = 1;
%         fprintf('Computing ilu(%g) of D...', setup.droptol);
        %tic;

        tic; 
        %[LD,UD] = ilu(D,setup); 
%         [LD,UD] = lu(D);
        [UD,pos_def] = chol(D);
        if(pos_def == 0)
            LD = UD';
        else
           [LD,UD] = lu(D);
        end
        t_factor_d = toc;

        fprintf('done!!!\n');      

       

  J = blkdiag(D,G);
%    keyboard;
tic;
  [LJ,UJ] = lu(J);

%     [L_block_G,U_block_G] = lu(G);
  t_jacobi = toc
 
%    prec_rhs = ones(m,1);
%     prec_sol = nssolve2(prec_rhs);
%     prec_sol = UG\(LG\prec_rhs);
%     pp = J * prec_sol;
%     pp = GS*prec_sol;
%     fp = fopen("~/rhs_MSC_20.txt","w");
%     fprintf(fp,"%d\n",prec_rhs);
%     fclose(fp);
%     keyboard;
%     prec_sol_C = load("~/MSC_sol_20.txt");
%     qq = J * prec_sol_C;
%     keyboard;
 clear J D G ; 
            %% Solve with PCG
            sol=zeros(n,1);
            tol = 1e-2; maxit = 3;restart = 40;
            max_pcg = 200;
            try
                fprintf('Enter GMRES...\n');
                %tic, [x,flag,relres,iter] = pcg(B,b,tol,maxit,@nssolve2); t_pcg = toc;
                %[x,flag,relres,iter,resvec] = pcg(B,b,tol,max_pcg,@jacobi_solve);t_gmres = toc;
                %tic, [x_np,flag_np,relres_np,iter_np,resvec_np] = gmres(B,b,restart,tol,maxit); t_gmres_np = toc;
                tic, [x,flag,relres,iter,resvec] = gmres(B,b,restart,tol,maxit,@nssolve2); t_gmres = toc;
%                 tic, [x,flag,relres,iter,resvec] = gmres(B,b,restart,tol,maxit,@jacobi_solve); t_gmres = toc;
                %tic, [x,flag,relres,iter] = gmres(B,b,restart,tol,maxit); t_gmres = toc;
                %% Display output
                %its_np = (iter_np(1)-1)*restart+iter_np(2);
                its = (iter(1)-1)*restart+iter(2);
                %its = iter;
                %fprintf('flag_np: %d\nits_np: %d\nrelres_np: %d\n', flag_np, its_np, relres_np);
                fprintf('flag: %d\nits: %d\nrelres: %d\n', flag, its, relres);
                %fprintf('time lu MSC: %g\ntime ilu D: %g\ntime PCG: %g\n', t_factor_msc, t_factor_d, t_pcg);
                %t_total_np = t_factor_msc + t_factor_d + t_gmres_np;
%                 t_total = t_factor_msc + t_factor_d + t_gmres;
                t_total = t_gmres + t_jacobi ;
                fprintf('total time: %g\n', t_total);
                fprintf('time lu MSC: %g\ntime ilu D: %g\ntime GMRES: %g\n', t_factor_msc, t_factor_d, t_gmres);
            catch
                fprintf('ERROR in GMRES aborting\n');break;
            end
            %fprintf('\nDifference in residual : %g\n',norm(xx-x));
%             q_np = B*x_np -b;
%             nrm_q = norm(q_np);
%             q = B*x -b;
%             prec_q = nssolve2(q);
%              prec_b = nssolve2(b);
%              for kk = 1: 10
%                fprintf("\nx[%d] = %f\n",kk,x(kk,:));
%              end
% % %             
%               keyboard;
            clear LD UD LG UG x  
            
%             iters_np(z) = its_np;
%             gmres_time_np(z) = t_gmres_np;
%             rel_res_np(z) = relres_np;
%             total_time_np(z) = t_total_np;
%             conv_flag_np(z) = flag_np;
            
            iters(z) = its;
            gmres_time(z) = t_gmres;
            rel_res(z) = relres;
            total_time(z) = t_total;
            conv_flag(z) = flag;

%             fprintf('\n\nNormal solve time : %g\n\n',normal_solve_time);
%             fprintf('\n\nPCG solve time : %g\n\n',t_total);
    end
    clear b
    
%     fprintf('tests done!!!\n');
    
%     fprintf('\n===========================================GMRES NO PRECONDITIONING==============================================\n');
%     fprintf('\nMSC block size | Restarts | Iter | normal_solve_time(s) | gmres_time(s) | Rel Res(gmres) | Total Time(s)  |  Flag \n');
%     fprintf('\n=================================================================================================================\n');
%     for i = 1 : length(nmsc_blocks)
%         fprintf('      %d        |    %d   |   %d  |      %g         |  %g   |     %g   |     %g  |     %d \n',nmsc_blocks(i),restart,iters_np(i),normal_solve_time,...
%                                                             gmres_time_np(i),rel_res_np(i),total_time_np(i),conv_flag_np(i));
%     end
    

    fprintf('\n===========================================GMRES MINI SCHUR PRECONDITIONING======================================\n');
    fprintf('\nMSC block size | Restarts | Iter | normal_solve_time(s) | gmres_time(s) | Rel Res(gmres) | Total Time(s)  |  Flag \n');
    fprintf('\n=================================================================================================================\n');
    for i = 1 : length(nmsc_blocks)
        fprintf('      %d        |    %d   |   %d  |      %g         |  %g   |     %g   |     %g  |     %d \n',nmsc_blocks(i),restart,iters(i),normal_solve_time,...
                                                            gmres_time(i),rel_res(i),total_time(i),conv_flag(i));
    end
    
    %% Functions follow ...
    function xx = nssolve2(y)
      y1 = y(1:sizeD); y2 = y(sizeD+1:sizeD+sizeG);z1 = UD\(LD\y1); 
      z2 = UG\(LG\(y2 - L*z1)); 
      %x2 = z2; x1 = z1 - UD\(LD\(U * x2)); 
      %xx = [x1; x2];
      xx = [z1; z2];
      %sum(isnan(x))
    end

    function xx = jacobi_solve(y)
        xx = UJ\(LJ\y);
%          y1 = y(1:sizeD); y2 = y(sizeD+1:sizeD+sizeG);
%          z1 = UD\(LD\y1);
%          z2 = U_block_G\(L_block_G\y2);
%          xx = [z1;z2]; 
    end

%     function xx = jacobi_split_solve(y)
%          y1 = y(1:sizeD); y2 = y(sizeD+1:sizeD+sizeG);
%          z1 = UD\(LD\y1);
%          z2 = U_block_G\(L_block_G\y2);
%          xx = [z1;z2];
%     end
end
