%% test driver to solve NS problem using overlapping Schur
%% complement method

function x1 = test_nssolve


clear all, clc
P = load('../extern/ifiss3.5/datafiles/square_stokes_nobc.mat');
addpath ../extern/MATLAB_DEMOS/armsC
np = size(P.B, 1);
A = [P.A P.B'; P.B zeros(np,np)]; nC = size(A,1);
%b = [P.f; P.g]; 
b = rand(nC,1); 

[rperm, nD] = indset(A,0,100); A = A(rperm,rperm); 
%nD = size(P.A,1);

%spy(A),

b = b(rperm); 

[rscal, A] = equil2(A); b = rscal' .* b;

%sum(rscal)


nG = nC-nD; D = A(1:nD,1:nD); G = A(nD+1:nC, nD+1:nC); 
E = A(1:nD, nD+1:nC); F = A(nD+1:nC, 1:nD);

nmsc = 4; r = rem(nG, nmsc); sz = (nG - r)/nmsc; r1 = 1; r2 = sz;

% find neighbors like this: neighbors_g{i} =  find(A(i,:))
% construct G_i
% for j = 1 : length(neighbors_g{i})
%   for k = 1 : length(neighbors_g{i})
%      G(i,j) = A(neighbors_g{i}(j), neighbors_g{i}(k));
%   end
% end


for i = 1 : nmsc
    PD{i} = A(r1:r2, r1:r2); PE{i} = E(r1:r2, r1:r2);
    PF{i} = F(r1:r2, r1:r2); PG{i} = G(r1:r2, r1:r2);
    r1 = r2+1; r2 = r2+sz;
end




for i = 1 : nmsc
    %det(PD{i})
    %S{i} = PG{i} - PF{i} * inv(PD{i})*(PE{i});
    S{i} = PG{i} - PF{i} * (PD{i} \ PE{i});
    %det(S{i})
    %  [LG{i}, UG{i}] = lu(S{i});
end


r2 = nG; last = nmsc+1;
if (r ~= 0)
    PD{last} = A(r1:r2, r1:r2); PE{last} = E(r1:r2, r1:r2);
    PF{last} = F(r1:r2, r1:r2); PG{last} = G(r1:r2, r1:r2);
    r1 = r2+1; r2 = r2+sz;
    %S{last} = PG{last} - PF{last}*inv(PD{last})*(PE{last});
    S{last} = PG{last} - PF{last} * (PD{last} \ PE{last});
    % [LG{last}, UG{last}] = lu(S{last});
end

GS = [];
for i = 1: nmsc
    GS = blkdiag(GS, S{i}); GS = sparse(GS);
end
GS = blkdiag(GS, S{last}); GS = sparse(GS); [LG,UG] = lu(GS); clear GS S
clear PD PE PF PG
precfun = @nssolve; sol = zeros(nC,1); setup.type = 'ilutp'; ...
          setup.droptol = 0.0001;
      fprintf('ilu ...\n');
[LD,UD] = ilu(D,setup);
clear D G
fprintf('Enter GMRES ...\n');

%[sol,res,its] = fgmres2(A,nD,nG,precfun,b,sol,E,F,LD,UD,LG,UG); 
 [x,flag,relres,iter] = gmres(A,b,100,1e-9,60,@nssolve2);
% 

%fprintf('test\n')

fprintf('outer_iter = %d\ninner_iter=%d\nflag = %d\nrelres=%g\n', iter(1), iter(2), flag, relres)




function x = nssolve2(y)
    y1 = y(1:nD); 
    y2 = y(nD+1:nD+nG);  
    z1 = UD\(LD\y1); 
    z2 = UG\(LG\(y2 - F*z1)); % z2 = P'*(UG\(LG\(P*(y2 - F*z1))));
    x2 = z2; 
    x1 = z1 - UD\(LD\(E * x2));
    x = [x1; x2];
%sum(isnan(x))
end

end
