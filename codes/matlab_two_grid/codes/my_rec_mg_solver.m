function [u, ucache, nu, err_vec] = my_rec_mg_solver(A, b, u0, PP, maxnu, tol, varargin)

m = length(PP)+1;
AA = cell(1,m);
AA{m} = A;
for k=m:-1:2
    AA{k-1} = PP{k-1}'*AA{k}*PP{k-1};
end

use_cache = (length(varargin)>0);
if (use_cache)
    ucache = zeros(size(A,1), maxnu+1);
    ucache(:,1) = u0;
end

u = u0;
for nu=1:maxnu
    r = b - A*u;
    %u = u + mg_step(r, AA, PP);
    %keyboard;
    u = u + recursive_Vcycle(m, AA, PP, r); 
    if (use_cache)
        ucache(:,nu+1) = u;
    end
    err_vec(nu) = norm(r)/norm(b);
    %fprintf('Error: %g\n', err_vec(nu));
    if err_vec(nu) <= tol
        break;
    end
end
