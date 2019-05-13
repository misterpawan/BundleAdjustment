%% test driver to solve NS problem using overlapping Schur compl. method
function test_nssolve_metis
clear all; clc;
fout = fopen('metis_results64_leaky_uniform_q1q2_re1000_mxit_20_tri.txt', 'a+');
P = load('/home/pawan/Work/Soft/ifiss3.2/datafiles/matlab.mat');
addpath /home/pawan/Work/Soft/MATLAB_DEMOS/armsC
addpath /home/pawan/Work/Matrices
addpath /home/pawan/Work/Codes/reordering
%% assemble matrix and rhs
A = [P.Anst P.Bst'; P.Bst zeros(P.np,P.np)]; nC = size(A,1); b = [P.fst; P.gst];   % keyboard     
[rscal, A] = equil2(A); b = rscal' .* b; clear rscal;
tic, fprintf('writing matrix...\n'); B = A + A'; 
mmwrite('/home/pawan/Work/Codes/nssolve/matlab.mtx', B);  clear B;
system('/home/pawan/Work/Soft/sparse_matrix_converter/sparse_matrix_converter matlab.mtx "MM" matlab.mtx "MM" -e');
matrix = '/home/pawan/Work/Codes/nssolve/matlab.mtx'; 
fprintf('done in %g secs\n', toc);

%% Set parameters
ntol = [ 0.1 0.01 0.001 0.0001 0.00001 0.000001]; ndoms = [8]; nblocksS = [4, 8]; 
m = 200; tol = 1e-8; maxit = 20;

%% Do tests
for jj = 1 : length(ndoms)                             
   fprintf(fout, 'Reordering ...'); nparts = ndoms(jj);
   tic, fprintf('doing metis reordering...'); 
   [rperm, sizes] = dd_ordering(matrix, nparts); 
   fprintf('done in %g secs\n', toc);
   A = A(rperm, rperm); nD = sum(sizes(1:nparts)); t_ind = toc;
   fprintf(fout, 'done!!!\n'); b = b(rperm); clear rperm; 
   nG = nC-nD; D = A(1:nD,1:nD); G = A(nD+1:nC, nD+1:nC); 
   E = A(1:nD, nD+1:nC); F = A(nD+1:nC, 1:nD);
   fprintf(fout, 'nD = %d, nG = %d\n', nD, nG);
     
   for ii = 1 : length(ntol)
     setup.type='ilutp'; setup.droptol = ntol(ii); 
     tic, 
     [LD,UD] = ilu(D,setup); 
     %[LD,UD] = lu(D); 
     t_factor_d = toc; fprintf(fout, 'done!!!\n');            
     for kk = 1 : length(nblocksS)      
       fprintf(fout, '\n\n================== %s ====================\n\n\n', date);
       fprintf(fout, 'pde: %d, domain: %d, grid_type: %d, qmethod: %d\n', P.pde, P.domain, P.grid_type, P.qmethod);
       fprintf(fout, 'viscosity: %g, nlmethod: %d, maxit_p: %d, maxit_n: %d, tol_nl: %g\n', P.viscosity, P.nlmethod, P.maxit_p, P.maxit_n, P.tol_nl);
       fprintf(fout, 'ntol: %d, ndoms: %d, nblocksS: %d \n', ntol(ii), ndoms(jj), nblocksS(kk));      
       fprintf('\n\nntol: %d, ndoms: %d, nblocksS: %d \n', ntol(ii), ndoms(jj), nblocksS(kk));      
       %% Identify MSCs
       nmsc = nblocksS(kk); r = rem(nG, nmsc); sz = (nG - r)/nmsc; r1 = 1; r2 = sz;
       ol = floor(sz); GS = []; %ol = 0;
       if (ol > sz), fprintf(fout, 'overlap size > msc size, taking ol=sz\n'); ol = floor(sz); end 
         fprintf(fout, 'Computing Schur complements ...');
         tic, fprintf('Computing Schur complements ...');
         for i = 1 : nmsc-2
           clear PD PE PF PG
           %PD = A(r1:r2+ol, r1:r2+ol); 
           PE = E(:, r1:r2+ol); 
           PF = F(r1:r2+ol, :); 
           PG = G(r1:r2+ol, r1:r2+ol);
           S{i} = PG - (PF * (UD \ (LD \ PE))); 
           GS(r1:r2+ol, r1:r2+ol) = S{i}; 
           %GS = sparse(GS); 
           r1 = r2+1; r2 = r2+sz;
         end
           oll = 0; % overlap size for 2nd last msc
           i = i+1;
           if (i == nmsc-1)
             clear PD PE PF PG
             %PD = A(r1:r2+oll, r1:r2+oll); 
             PE = E(:, r1:r2+oll);
             PF = F(r1:r2+oll, :); 
             PG = G(r1:r2+oll, r1:r2+oll);
             S{i} = PG - (PF * (UD \ (LD \ PE))); 
             GS(r1:r2+oll, r1:r2+oll) = S{i};
             %GS = sparse(GS); 
             r1 = r2+1; 
           end
             r2 = nG; i = i + 1;
             if (i == nmsc)
               clear PD PE PF PG
               %PD = A(r1:r2, r1:r2); 
               PE = E(:, r1:r2);
               PF = F(r1:r2, :); 
               PG = G(r1:r2, r1:r2);
               S{i} = PG - (PF * (UD \ (LD \ PE))); 
               GS(r1:r2, r1:r2) = S{i};
               %GS = sparse(GS); 
            end
            fprintf('done in %g secs\n', toc);
            fprintf(fout, 'done!!!\n');
            
            %% Factor MSC
            fprintf(fout, 'Computing LU of MSC ...');tic,[LG,UG] = lu(GS);t_factor_msc = toc;fprintf(fout, 'done!!!\n');
            clear GS S PD PE PF PG
            
            %% Solve with GMRES
            try
                %clear D G; 
                fprintf('Enter GMRES ...');
                %[sol,res,its] = fgmres2(A,nD,nG,precfun,b,sol,E,F,LD,UD,LG,UG); 
                fprintf(fout, 'Enter GMRES ...');
                tic, [~,flag,relres,iter] = gmres(A,b,m,tol,maxit,@nssolve2); t_gmres = toc;
                fprintf(fout, 'done!!!\n');
                fprintf('done!!!\n');
                %% Display output
                its = (iter(1)-1)*m + iter(2);
                fprintf(fout, 'flag: %d\nits: %d\nrelres: %d\n', flag, its, relres);
                fprintf('flag: %d\nits: %d\nrelres: %d\n', flag, its, relres);
                fprintf(fout, 'time lu MSC: %g\ntime ilu D: %g\ntime GMRES: %g\ntime ind. ord.: %g\n', t_factor_msc, t_factor_d, t_gmres, t_ind);
                t_total = t_factor_msc + t_factor_d + t_gmres;
                fprintf(fout, 'total time: %g\ntotal all: %g\n', t_total, t_total+t_ind);
                %fprintf('time lu MSC: %g\ntime ilu D: %g\ntime GMRES: %g\n', t_factor_msc, t_factor_d, t_gmres);
            catch
                fprintf(fout, 'ERROR in ilu(%g) aborting\n', setup.droptol);
            end
            clear LG UG  
      end
    end
 clear D E F G rperm sizes
end
fprintf('tests done!!!\n');

%% Functions follow ...
function xx = nssolve2(y)
    y1 = y(1:nD); y2 = y(nD+1:nD+nG);z1 = UD\(LD\y1); 
    z2 = UG\(LG\(y2 - F*z1)); 
    %x2 = z2; x1 = z1 - UD\(LD\(E * x2)); 
    %xx = [x1; x2];
    xx = [z1; z2];
end
end
