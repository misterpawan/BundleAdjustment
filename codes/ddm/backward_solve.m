function x = backward_solve(offset, DU, DL, DP, U, rhs)
   nparts = length(offset)-2;
   x{nparts+1} = rhs{nparts+1};
     % keyboard;

   for i = nparts : -1 : 1 
        p = DP{i};
        yy = (U{i}*x{nparts+1});
        x{i} = rhs{i} - (DU{i} \ (DL{i} \ yy(p)));
   end
end