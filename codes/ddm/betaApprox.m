function beta1 = betaApprox(U,i,D)
   UtOnes = U{i}*ones(size(U{i}, 2), 1);
   tmp = D{i} \ UtOnes;
   beta1 = tmp; %zeros(size(UtOnes,1),1);
   indNotZeros = (abs(UtOnes) > 1e-14);
   beta1(indNotZeros) = tmp(indNotZeros) ./ UtOnes(indNotZeros);
end