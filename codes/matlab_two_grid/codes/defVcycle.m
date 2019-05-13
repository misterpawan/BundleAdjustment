function x = defVcycle(A, level, k, x, r, itmax, tol, levmax, matrix)
    if level == levmax 
        x = A \ r;
    else
        %[x, relres, W] = srcg3b(A, r, x, k, itmax, tol); % relax with SRCG
        x = (diag(diag(A)) + tril(A)) \ x; W = [];
        fprintf('level: %d, pre-smoothing done \n', level);
%         if (level == 1) % first aggregation by metis
%             cf = 3.5; n = (size(A,1))^(1/3); nc = int2str((floor(n/cf))^3), metis_aggregation_fun(matrix, nc); 
%             fprintf('level: %d, metis aggregation done\n', level);
%             
%             cd /home/pawan/Work/Codes/gmg/gmg_c_new/
%             ind = load('part_array.txt'); fprintf('level: %d, size of ind: %d\n', level, length(ind));
%             cd /home/pawan/Work/Codes/amg/
%             max(ind)
%             ind2 = find(ind); I = sparse(ind2, ind(ind2), sign(ind(ind2)), size(A, 1), max(ind)); clear ind ind2;
%             %ind2 = find(ind); P = sparse(ind2, ind(ind2), sign(ind(ind2)), size(A, 1), max(ind)); clear ind ind2; % aggregation to prolongation
%          
%             fprintf('level: %d, size of metis I %d x %d\n', level, size(I,1), size(I,2));
%         else
            [I ind] = agtwolev(A, 2, 1, 2); clear ind;
            fprintf('level: %d, size of I %d x %d\n', level, size(I,1), size(I,2));
%        end
        cd /home/pawan/Work/Codes/amg ; 
        It = [I W]; % enrich interpolation
        Ac = It' * A * It; % construct coarse grid
        fprintf('level: %d, coarse grid constructed\n', level);
        rc = It' * ( r - A * x); % restrict residual
        xc = zeros(length(rc), 1); % coarse grid x
        fprintf('level: %d, size of Ac: (%d, %d), xc: %d, rc: %d\n', level, size(Ac,1), size(Ac,2), length(xc), length(rc));
        xc = defVcycle(Ac, level+1, k, xc, rc, itmax, tol, levmax, matrix); % recursively apply Vcycle
        fprintf('level: %d, Vcycle done\n', level);
        x = x + It * xc;
        %[x, relres, W] = srcg3b(A, r, x, k, itmax, tol); % relax with SRCG
        x = (diag(diag(A)) + triu(A)) \ x;
        %x = (diag(diag(A)) + tril(A)) \ x; W = [];
        fprintf('level: %d, post-smoothing done \n', level); clear W It xc rc A r
    end
end