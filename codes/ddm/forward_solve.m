function x = forward_solve(offset, DU, DL, DP, L, TU, TL, TP, rhs)

    ncols = length(rhs);
    nparts = length(offset)-2;
      
    Lx = zeros(ncols - offset(nparts+1) + 1, 1);
        
    for i = 1 : nparts
       p = DP{i};
       yy = rhs(offset(i):offset(i+1)-1);
       x{i} = DU{i} \ (DL{i} \ yy(p));
    end
    
    for i = 1 : nparts
        Lx = Lx + L{i}*x{i};
    end

    zz = rhs(offset(nparts+1):ncols) - Lx;
    x{nparts+1} = TU \ (TL \ zz(TP));
end
