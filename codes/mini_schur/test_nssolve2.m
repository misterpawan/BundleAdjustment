%% test driver to solve NS problem using overlapping Schur
%% complement method
function x1 = test_nssolve2
clear all, clc
P = load('/home/pawan/Work/Soft/ifiss3.2/datafiles/matlab.mat');
addpath /home/pawan/Work/Soft/MATLAB_DEMOS/armsC
A = [P.Anst P.Bst'; P.Bst zeros(P.np,P.np)]; nC = size(A,1);
b = [P.fst; P.gst]; [rperm, nD] = indset(A,0,400); A = A(rperm,rperm); 
b = b(rperm); [rscal, A] = equil2(A); b = rscal' .* b;
nG = nC-nD; D = A(1:nD,1:nD); G = A(nD+1:nC, nD+1:nC); 
E = A(1:nD, nD+1:nC); F = A(nD+1:nC, 1:nD);

%% Identify MSCs
nmsc = 10; r = rem(nG, nmsc); sz = (nG - r)/nmsc; r1 = 1; r2 = sz;
ol = floor(sz/2);

for i = 1 : nmsc-1
    PD{i} = A(r1:r2+ol, r1:r2+ol); PE{i} = E(r1:r2+ol, r1:r2+ol);
    PF{i} = F(r1:r2+ol, r1:r2+ol); PG{i} = G(r1:r2+ol, r1:r2+ol);
    r1 = r2+1; r2 = r2+sz;
end

oll = (nG - (r2 + 1))/2; % overlap size for 2nd last msc

if (i == nmsc)
    PD{i} = A(r1:r2+oll, r1:r2+oll); PE{i} = E(r1:r2+oll, r1:r2+oll);
    PF{i} = F(r1:r2+oll, r1:r2+oll); PG{i} = G(r1:r2+oll, r1:r2+oll);
    r1 = r2+1; r2 = r2+sz;
end
%% Compute MSC
for i = 1 : nmsc-1
    %det(PD{i})
    %S{i} = PG{i} - PF{i} * inv(PD{i})*(PE{i});
    S{i} = PG{i} - PF{i} * (PD{i} \ PE{i});
    GS(r11:r22+ol) = S{i};
    r11 = r22+1; r22 = r22+sz;
    %det(S{i})
    %  [LG{i}, UG{i}] = lu(S{i});
end

i = nmsc; 
S{i} = PG{i} - PF{i} * (PD{i} \ PE{i});
GS(r11:r22+oll) = S{i};
r11 = r22+1; r22 = r22+sz;

r2 = nG; last = nmsc+1; r22 = nG;
if (r ~= 0)
    PD{last} = A(r1:r2, r1:r2); PE{last} = E(r1:r2, r1:r2);
    PF{last} = F(r1:r2, r1:r2); PG{last} = G(r1:r2, r1:r2);
    r1 = r2+1; r2 = r2+sz;
    %S{last} = PG{last} - PF{last}*inv(PD{last})*(PE{last});
    S{last} = PG{last} - PF{last} * (PD{last} \ PE{last});
    GS(r11:r22) = S{i};    
    % [LG{last}, UG{last}] = lu(S{last});
end

%% Factor MSC
GS = sparse(GS); [LG,UG] = lu(GS); clear GS S
clear PD PE PF PG
%% Solve with GMRES
precfun = @nssolve; sol = zeros(nC,1); setup.type = 'ilutp'; ...
          setup.droptol = 0.0001;
      fprintf('ilu ...\n');
[LD,UD] = ilu(D,setup);
clear D G
fprintf('Enter GMRES ...\n')
%[sol,res,its] = fgmres2(A,nD,nG,precfun,b,sol,E,F,LD,UD,LG,UG); 
[x,flag,relres,iter] = gmres(A,b,100,1e-9,1,@nssolve2);

iter
flag
relres

function x = nssolve2(y)
    y1 = y(1:nD); y2 = y(nD+1:nD+nG);  
    z1 = UD\(LD\y1); z2 = UG\(LG\(y2 - F*z1)); 
    x2 = z2; x1 = z1 - UD\LD\(E * x2);
    x = [x1; x2];
%sum(isnan(x))
end

end
