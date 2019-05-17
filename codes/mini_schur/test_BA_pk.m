%This file tests the domain_decomposition.m function
%function test_BA

function test

clc; warning('off','all'); addpath('../ddm');
addpath(genpath('/media/pawan/Data/datasets'));

lhs_filename = 'JTJ49.mat';
P = load(lhs_filename); B = P.A;
[m,n] = size(B); b = rand(m,1);
 
tic; xx = normal_solve(B, b); normal_solve_time = toc;
nparts = 5;
  
[sizes_parts, A, p] = domain_decomposition(B, nparts);
sizeG = m - sum(sizes_parts(1:nparts));
sizeD = m - sizeG; 
    
D = A(1 : sizeD, 1 : sizeD);
L = A(sizeD + 1: m, 1 : sizeD);
U = L';
G = A(sizeG + 1:m, sizeG + 1:m);
    
%% Identify MSCs
% blocks for schur complements
nmsc = 3; r = rem(sizeG, nmsc); sz = (sizeG - r)/nmsc; r1 = 1; r2 = sz;
GS = []; %ol = 0;
   
%nmsc
%sz

% blocks for D: good to take same number of blocks for D as for G
nD = nmsc; rD = rem(sizeD, nD); szD = (sizeD - rD)/nD; rD1 = 1; rD2 = szD;

%nD
%szD
    
% overlap size
%ol = floor(sz);
ol = sz; % try with zero overlap

%keyboard;
    
fprintf('Computing Mini Schur complements ...\n');
for i = 1 : nmsc-2
   clear PD PU PL PG
   PD = D(rD1 : rD2 + ol, rD1 : rD2 + ol); 
   PU = U(rD1 : rD2 + ol, r1 : r2 + ol); 
   PL = PU'; % for BA L = U'
   PG = G(r1 : r2 + ol, r1 : r2 + ol);
   
  % printf block sizes
   fprintf('\nPrinting block sizes for mini Schur complement: %d\n', i);
   fprintf('block size of G = %d x %d (%d : %d, %d : %d) \n', size(PG,1), size(PG,2), r1,r2, r1, r2);
   fprintf('block size of U = %d x %d (%d : %d, %d : %d) \n', size(PU,1), size(PU,2), rD1, rD2, r1, r2);
   fprintf('block size of L = %d x %d (%d : %d, %d : %d) \n', size(PL,1), size(PL,2), r1, r2, rD1, rD2);
   fprintf('block size of D = %d x %d (%d : %d, %d : %d) \n\n', size(PD,1), size(PD,2), rD1, rD2, rD1, rD2);
   
   S{i} = PG - PL * (PD \ PU); 
   GS(r1:r2 + ol, r1 : r2 + ol) = S{i}; 
   GS = sparse(GS); 
   r1 = r2 + 1; 
   r2 = r2 + sz;
   rD1 = rD2 + 1;
   rD2 = rD2 + szD;
end
    
oll = 0; % overlap size for 2nd last msc
i = i+1;
    
if (i == nmsc-1)
   clear PD PU PL PG
   PD = A(rD1:rD2+oll, rD1:rD2+oll); 
   PU = U(rD1:rD2+oll, r1:r2+oll);
   %PL = L(r1:r2+oll, r1:r2+oll); 
   PL = PU';
   PG = G(r1:r2+oll, r1:r2+oll);
   
  % printf block sizes
   fprintf('Printing block sizes for mini Schur complement: %d\n', i);
   fprintf('block size of G = %d x %d (%d : %d, %d : %d) \n', size(PG,1), size(PG,2), r1,r2, r1, r2);
   fprintf('block size of U = %d x %d (%d : %d, %d : %d) \n', size(PU,1), size(PU,2), rD1, rD2, r1, r2);
   fprintf('block size of L = %d x %d (%d : %d, %d : %d) \n', size(PL,1), size(PL,2), r1, r2, rD1, rD2);
   fprintf('block size of D = %d x %d (%d : %d, %d : %d) \n\n', size(PD,1), size(PD,2), rD1, rD2, rD1, rD2);

   S{i} = PG - PL * (PD \ PU); 
   GS(r1:r2+oll, r1:r2+oll) = S{i};
   GS = sparse(GS); 
   r1 = r2+1;
   rD1 = rD2+1;
end
    
r2 = sizeG; i = i + 1;
rD2 = sizeD; 

if (i == nmsc)
  clear PD PU PL PG
  PD = A(rD1:rD2, rD1:rD2); 
  PU = U(rD1:rD2, r1:r2);
  %PL = L(r1:r2, r1:r2); 
  PL = PU';
  PG = G(r1:r2, r1:r2);
  
   % printf block sizes
   fprintf('Printing block sizes for mini Schur complement: %d\n', i);
   fprintf('block size of G = %d x %d (%d : %d, %d : %d) \n', size(PG,1), size(PG,2), r1,r2, r1, r2);
   fprintf('block size of U = %d x %d (%d : %d, %d : %d) \n', size(PU,1), size(PU,2), rD1, rD2, r1, r2);
   fprintf('block size of L = %d x %d (%d : %d, %d : %d) \n', size(PL,1), size(PL,2), r1, r2, rD1, rD2);
   fprintf('block size of D = %d x %d (%d : %d, %d : %d) \n\n', size(PD,1), size(PD,2), rD1, rD2, rD1, rD2);

  
  S{i} = PG - PL * (PD \ PU); 
  GS(r1:r2, r1:r2) = S{i};
  GS = sparse(GS); 
end 
       
fprintf('done with mini Schur complements!\n'); 
fprintf('Computing LU of MSC ...');
tic,[LG,UG] = lu(GS);t_factor_msc = toc; fprintf('done!!!\n');
clear GS S PD PU PL PG
        
setup.type='ilutp'; 
setup.droptol = 1e-03; %ntol(ii); 
fprintf('Computing ilu(%g) of D...', setup.droptol);
%tic;
        
tic; 
%[LD,UD] = ilu(D,setup); 
[LD,UD] = lu(D);
t_factor_d = toc;
      
fprintf('done!!!\n'); clear D G ;
       
%condest(B)
        
%% Solve with PCG
%precfun=@nssolve; 
sol=zeros(m,1); tol = 1e-3; maxit = 100;
restart = 60;
%try
  fprintf('Enter GMRES...\n');
  %tic, [x,flag,relres,iter] = pcg(B,b,tol,maxit,@nssolve2); t_pcg = toc;
  tic, [x,flag,relres,iter] = gmres(B,b,restart,tol,maxit,@nssolve2); t_pcg = toc;

%% Display output
  %its = iter;
  its = (iter(1)-1)*restart + iter(2);
  fprintf('flag: %d\nits: %d\nrelres: %d\n', flag, its, relres);
%fprintf('time lu MSC: %g\ntime ilu D: %g\ntime PCG: %g\n', t_factor_msc, t_factor_d, t_pcg);
  t_total = t_factor_msc + t_factor_d + t_pcg;
  fprintf('total time: %g\n', t_total);
  fprintf('time lu MSC: %g\ntime ilu D: %g\n\ntime PCG: %g\n', t_factor_msc, t_factor_d, t_pcg);
%catch
%  fprintf('ERROR in pcg aborting\n');
%end
%fprintf('\nDifference in residual : %g\n',norm(xx-x));
  clear LD UD LG UG x b     
  fprintf('Normal solve time : %g\n\n',normal_solve_time);
%  fprintf('\n\nPCG solve time : %g\n\n',t_total);   
  fprintf('tests done!!!\n');
   
%% Functions follow ...
function xx = nssolve2(y)
  y1 = y(1:sizeD); y2 = y(sizeD+1:sizeD+sizeG);z1 = UD\(LD\y1); 
  z2 = UG\(LG\(y2 - L*z1)); 
  %x2 = z2; x1 = z1 - UD\(LD\(U * x2)); 
  %xx = [x1; x2];
  xx = [z1; z2];
  %sum(isnan(x))
end

end

