function x = backward_solve(offset, D, U, rhs)
   nparts = length(offset)-2;
   x{nparts+1} = rhs(nparts+1);
   for i = nparts : 1 : -1 
        x{i} = rhs{i} - D{i} \ (U{i}*x{nparts+1});
   end
end