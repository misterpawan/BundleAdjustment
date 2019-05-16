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
    
    %load rhs
    b = load(rhs_filename);
    b = b.rhs;
    
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
        %dd = []; uu = []; ll = []; 
        for i = 1 : nparts
            L{i} = A(offset(nparts+1) : ncols, offset(i): offset(i+1)-1); %ll = [ll L{i}];
            U{i} = A(offset(i): offset(i+1)-1, offset(nparts+1) : ncols); %uu = [uu;U{i}];
            D{i} = A(offset(i): offset(i+1)-1, offset(i): offset(i+1)-1); %dd = [dd;D{i}];
        end
        %clear L D U
        %L = ll; U = uu; D = dd; clear ll uu dd
        %size(L{1})
        %size(U{1})
        %size(D{1})

        G = A(offset(nparts+1) : ncols, offset(nparts+1): ncols);
        %D{nparts+1} = G;
        nG = size(G,1);
        
        %identify MSCs
        nmsc = 3;r = rem(nG, nmsc);
        sz = (nG - r)/nmsc;
        r1_gs = 1; 
        r2_gs = sz;
        ol = 0;%floor(sz);
        GS = [];
        fprintf('Computing Schur complements ...');
        for kk = 1 : nparts
            %fprintf('\nkk = %d\n',kk);
            r1 = 1; 
            r2 = sz;
            
            for i = 1 : nmsc-2
                clear PD PU PL PG
                PD = D{kk}(r1:r2+ol, r1:r2+ol); PU = U{kk}(r1:r2+ol, r1:r2+ol); 
                PL = L{kk}(r1:r2+ol, r1:r2+ol); PG = G(r1:r2+ol, r1:r2+ol);
                S{i} = PG - PL * (PD \ PU); 
                GS(r1_gs:r2_gs+ol, r1_gs:r2_gs+ol) = S{i}; 
                %GS(r1:r2+ol, r1:r2+ol) = S{i}; 
                GS = sparse(GS); 
                r1 = r2+1; r2 = r2+sz;
                r1_gs = r2_gs+1;
                r2_gs = r2_gs+sz;
            end
            
            oll = 0; % overlap size for 2nd last msc
            i = i+1;
            
            if (i == nmsc-1)
                %fprintf('\nr1 : %d\n',r1);
                clear PD PU PL PG
                PD = D{kk}(r1:r2+oll, r1:r2+oll); 
                PU = U{kk}(r1:r2+oll, r1:r2+oll);
                PL = L{kk}(r1:r2+oll, r1:r2+oll); PG = G(r1:r2+oll, r1:r2+oll);
                S{i} = PG - PL * (PD \ PU); 
                GS(r1_gs:r2_gs+oll, r1_gs:r2_gs+oll) = S{i};
                %GS(r1:r2+ol, r1:r2+ol) = S{i}; 
                GS = sparse(GS); 
                r1 = r2+1; 
                r1_gs = r2_gs + 1;
            end
            
            r2 = nG; r2_gs = kk*nG;
            i = i + 1;
            if (i == nmsc)
                clear PD PU PL PG
                PD = D{kk}(r1:r2, r1:r2); PU = U{kk}(r1:r2, r1:r2);
                PL = L{kk}(r1:r2, r1:r2); PG = G(r1:r2, r1:r2);
                S{i} = PG - PL * (PD \ PU);
                GS(r1_gs:r2_gs, r1_gs:r2_gs) = S{i};
                %GS(r1:r2, r1:r2) = S{i};
                GS = sparse(GS); 
            end
            
            r1_gs = r2_gs+1;
            r2_gs = r2_gs+sz;
        end
        
        
        
        fprintf('done!!!\n');
        
        spy(GS);
        %mgs = size(GS,1)
        %ngs = size(GS,2)
        
        fprintf('Computing LU of MSC ...');tic,[LG,UG] = lu(GS);t_factor_msc = toc;
        fprintf('done!!!\n');
        clear GS S PD PU PL PG
        
        setup.type='ilutp'; m = 500;
        setup.droptol = 1e-03; %ntol(ii); 
        fprintf('Computing ilu(%g) of D...', setup.droptol);
        %tic;
        dd = []; %for storing the D block
        ll = []; %for storing the L block
        row1 = 1;
        row2 = 0;
        for i = 1: nparts
            row2 = row2 + size(D{i},1);
            dd(row1:row2,row1:row2) = D{i};
            ll(1:nG,row1:row2) = L{i};
            row1 = row2+1;
        end
        clear D L
        D = sparse(dd);%spy(D)
        L = sparse(ll);
        tic;
        [LD,UD] = ilu(D,setup); 
        t_factor_d = toc;
        nD = size(D,1);
        fprintf('done!!!\n');            
        clear D G ll;
        
        %% Solve with PCG
        precfun=@nssolve; sol=zeros(ncols,1);
        tol = 1e-7; maxit = 20;
        %try
            fprintf('Enter PCG...');
            tic, [x,flag,relres,iter] = pcg(A,b,tol,maxit,@nssolve2); t_pcg = toc;
            %% Display output
            its = (iter(1)-1)*m + iter(2);
            fprintf(fout, 'flag: %d\nits: %d\nrelres: %d\n', flag, its, relres);
            fprintf(fout, 'time lu MSC: %g\ntime ilu D: %g\ntime GMRES: %g\ntime ind. ord.: %g\n', t_factor_msc, t_factor_d, t_gmres, t_ind);
            t_total = t_factor_msc + t_factor_d + t_pcg;
            fprintf('total time: %g\ntotal all: %g\n', t_total, t_total);
            fprintf('time lu MSC: %g\ntime ilu D: %g\ntime PCG: %g\n', t_factor_msc, t_factor_d, t_pcg);
        %catch
         %   fprintf('ERROR in pcg aborting\n');
        %end
        clear LD UD LG UG x b 
        
        fprintf('tests done!!!\n');

        %% Functions follow ...
        function xx = nssolve2(y)
            y1 = y(1:nD); y2 = y(nD+1:nD+nG);z1 = UD\(LD\y1); 
            z2 = UG\(LG\(y2 - L*z1)); 
            %x2 = z2; x1 = z1 - UD\(LD\(E * x2)); 
            %xx = [x1; x2];
            xx = [z1; z2];
        %sum(isnan(x))
        end
end
    
    
    
    