%% test driver to solve NS problem using overlapping Schur compl. method
function xx = test_nssolve3
clear all,clc

ftex = fopen('tables.tex', 'a+');
fout = fopen('rresults64_leaky_stretched_q1q2_re100_mxit_20_tri.txt', 'a+');
P = load('/home/pawan/Work/Soft/ifiss3.2/datafiles/matlab.mat');
addpath /home/pawan/Work/Soft/MATLAB_DEMOS/armsC
 
ntol = [ 0.1 0.01 0.001 0.0001 0.00001 0.000001]; 
%ntol = [0.00001 0.000001 0.0000001];
nblocksA = [50 100 200 400 800 1600 3200 6400]; 
%nblocksA = [6400];
nblocksS = [3, 4, 8, 16, 32, 64]; 
%nblocksS = [32, 64];

for ii = 1 : length(ntol)
    for jj = 1 : length(nblocksA)
        for kk = 1 : length(nblocksS)           
            %% read matrix and rhs
            A = [P.Anst P.Bst'; P.Bst zeros(P.np,P.np)]; nC = size(A,1);
            b = [P.fst; P.gst];
            
            fprintf(fout, '\n\n================== %s ====================\n\n\n', date);
            fprintf(fout, 'pde: %d, domain: %d, grid_type: %d, qmethod: %d\n', P.pde, P.domain, P.grid_type, P.qmethod);
            fprintf(fout, 'viscosity: %g, nlmethod: %d, maxit_p: %d, maxit_n: %d, tol_nl: %g\n', P.viscosity, P.nlmethod, P.maxit_p, P.maxit_n, P.tol_nl);
            fprintf(fout, 'ntol: %d, nblocksA: %d, nblocksS: %d \n', ntol(ii), nblocksA(jj), nblocksS(kk));
   
            fprintf(fout, 'Reordering ...'); sA = nblocksA(jj);
            if(sA > nC/2) 
                fprintf(fout, 'Warning only one ind. block: sA = %d > %d = nC/2\n',sA,nC/2);
            end
            if (kk == 1)
                tic,[rperm,nD]=indset(A,0,sA);t_ind = toc;save('permut.mat','rperm','nD');
            else
                PP = load('permut.mat'); rperm = PP.rperm; nD = PP.nD;
            end
            A = A(rperm,rperm); fprintf(fout, 'done!!!\n'); 
            b = b(rperm); clear rperm; [rscal, A] = equil2(A); b = rscal' .* b; clear rscal;
            nG = nC-nD; D = A(1:nD,1:nD); G = A(nD+1:nC, nD+1:nC); 
            E = A(1:nD, nD+1:nC); F = A(nD+1:nC, 1:nD);
            fprintf(fout, 'nD = %d, nG = %d\n', nD, nG);
            
            %% Identify MSCs
            nmsc = nblocksS(kk); r = rem(nG, nmsc); sz = (nG - r)/nmsc; r1 = 1; r2 = sz;
            ol = floor(sz); GS = []; %ol = 0;
            if (ol > sz), fprintf(fout, 'overlap size > msc size, taking ol=sz\n'); ol = floor(sz); end 
            fprintf(fout, 'Computing Schur complements ...');
            for i = 1 : nmsc-2
                clear PD PE PF PG
                PD = A(r1:r2+ol, r1:r2+ol); PE = E(r1:r2+ol, r1:r2+ol); 
                PF = F(r1:r2+ol, r1:r2+ol); PG = G(r1:r2+ol, r1:r2+ol);
                S{i} = PG - PF * (PD \ PE); GS(r1:r2+ol, r1:r2+ol) = S{i}; 
                GS = sparse(GS); r1 = r2+1; r2 = r2+sz;
            end
            oll = 0; % overlap size for 2nd last msc
            i = i+1;
            if (i == nmsc-1)
                clear PD PE PF PG
                PD = A(r1:r2+oll, r1:r2+oll); PE = E(r1:r2+oll, r1:r2+oll);
                PF = F(r1:r2+oll, r1:r2+oll); PG = G(r1:r2+oll, r1:r2+oll);
                S{i} = PG - PF * (PD \ PE); GS(r1:r2+oll, r1:r2+oll) = S{i};
                GS = sparse(GS); r1 = r2+1; 
            end
            r2 = nG; i = i + 1;
            if (i == nmsc)
                clear PD PE PF PG
                PD = A(r1:r2, r1:r2); PE = E(r1:r2, r1:r2);
                PF = F(r1:r2, r1:r2); PG = G(r1:r2, r1:r2);
                S{i} = PG - PF * (PD \ PE); GS(r1:r2, r1:r2) = S{i};
                GS = sparse(GS); 
            end
            fprintf(fout, 'done!!!\n');
            
            %% Factor MSC
            fprintf(fout, 'Computing LU of MSC ...');tic,[LG,UG] = lu(GS);t_factor_msc = toc;fprintf(fout, 'done!!!\n');
            clear GS S PD PE PF PG
            
            %% Solve with GMRES
            precfun=@nssolve; sol=zeros(nC,1);setup.type='ilutp'; m = 500;tol = 1e-7; maxit = 20;
            setup.droptol = ntol(ii); 
            fprintf(fout, 'Computing ilu(%g) of D...', setup.droptol);
            try
                tic, [LD,UD] = ilu(D,setup); t_factor_d = toc; fprintf(fout, 'done!!!\n');            
                clear D G; %fprintf('Enter GMRES ...\n');
                %[sol,res,its] = fgmres2(A,nD,nG,precfun,b,sol,E,F,LD,UD,LG,UG); 
                fprintf(fout, 'Enter GMRES ...');
                tic, [x,flag,relres,iter] = gmres(A,b,m,tol,maxit,@nssolve2); t_gmres = toc;
                fprintf(fout, 'done!!!\n');
                %% Display output
                its = (iter(1)-1)*m + iter(2);
                fprintf(fout, 'flag: %d\nits: %d\nrelres: %d\n', flag, its, relres);
                fprintf(fout, 'time lu MSC: %g\ntime ilu D: %g\ntime GMRES: %g\ntime ind. ord.: %g\n', t_factor_msc, t_factor_d, t_gmres, t_ind);
                t_total = t_factor_msc + t_factor_d + t_gmres;
                fprintf(fout, 'total time: %g\ntotal all: %g\n', t_total, t_total+t_ind);
                %fprintf('time lu MSC: %g\ntime ilu D: %g\ntime GMRES: %g\n', t_factor_msc, t_factor_d, t_gmres);
            catch
                fprintf(fout, 'ERROR in ilu(%g) aborting\n', setup.droptol);
            end
            clear LD UD LG UG x b 
        end
    end
end
fprintf('tests done!!!\n');

%% Functions follow ...
function xx = nssolve2(y)
    y1 = y(1:nD); y2 = y(nD+1:nD+nG);z1 = UD\(LD\y1); 
    z2 = UG\(LG\(y2 - F*z1)); 
    %x2 = z2; x1 = z1 - UD\(LD\(E * x2)); 
    %xx = [x1; x2];
    xx = [z1; z2];
%sum(isnan(x))
end
end
