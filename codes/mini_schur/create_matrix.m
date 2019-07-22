clc
% A1 = rand(3,3);
% A2 = rand(3,3);
% A3 = rand(3,3);
% 
% G1 = rand(3,3);
% G2 = rand(3,3);
% G3 = rand(3,3);
% 
% A = blkdiag(A1,A2,A3,G1,G2,G3);
% [m,n]  = size(A);
% 
% A = A + A';

% [row,col,val] = find(A);
% fp = fopen("")

filepath = '~/Test/49';
addpath(filepath);

lhs_filename = 'JTJ49_1.mat';
P = load(lhs_filename,'-mat');
A = P.lhs;
[m,n] = size(A);

% A = full(A);
% for i = 10:18
%     for j = 1:9
%         A(i,j) = rand(1);
%     end
% end

% A(1:9,10:18) = A(10:18,1:9)';
% 
% fp= fopen("~/mini_schur/test/test_col.txt","w");
% for i = 1:18
%     if i <= 10
%         fprintf(fp,"%d\n",(i-1)*12);
%     else
%         fprintf(fp,"%d\n", 120+((i-11)*12));
%     end
%         
% end
% 
% fclose(fp);
% 
% fp_row = fopen("~/mini_schur/test/test_row.txt","w");
% fp_val = fopen("~/mini_schur/test/test_val.txt","w");
% for i = 1 : 18
%     for j = 1 :18
%         if A(i,j) ~= 0
%             fprintf(fp_row,"%d\n",j-1);
%             fprintf(fp_val,"%15.15f\n",A(i,j));
%         end
%     end
% end
% 
% fclose(fp_row);
% fclose(fp_val);

fprintf("\nDone..!\n");

   %nparts = 10;
    %[sizes_parts, A, p] = domain_decomposition(B, nparts);
%     A=B;
    %sizeG = m - sum(sizes_parts(1:nparts))
    sizeG = 441;
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
%     nmsc_blocks = [3 5 10 20 30];
    nmsc_blocks = [3];
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
           invDU = PD\PU;
           nnz(invDU)
           invDU = spconvert(invDU);
           keyboard;
           S{i} = PG - PL * invDU; 
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
           invDU = PD\PU;
           nnz(invDU)
           S{i} = PG - PL * invDU; 
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

          invDU = PD\PU;
          nnz(invDU)
          S{i} = PG - PL * (PD \ PU); 
          GS(r1:r2, r1:r2) = S{i};
          GS = sparse(GS); 
        end 
        %det_GS = det(GS)
        fprintf('done with mini Schur complements!\n'); 

    end



