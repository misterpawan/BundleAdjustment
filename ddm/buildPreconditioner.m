function [T, BigBeta] = buildPreconditioner(offset,L,D,U)

   nparts = length(offset)-2;
   T = D{nparts+1};
   BigBeta = [];

   for i = 1 : nparts
      Beta1{i} = sparse(diag(betaApprox(U,i,D)));
      T = T - (L{i} * (Beta1{i} * U{i})); 
      
      BigBeta = blkdiag(BigBeta, Beta1{i});
      
   end
   BigBeta = diag(BigBeta);

end