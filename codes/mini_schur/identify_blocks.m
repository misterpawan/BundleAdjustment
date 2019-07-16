matfile = 'JTJ372.mat';

P = load(matfile); A = P.lhs;

[m,n] = size(A);

p = 372*9;

min_row = m;

for i=1:m-p
    
    i
   S = find(A(i,i+3:n) > 0);
   %idx = S > i+2;
   min_row =  min(min_row, min(S));
  
end

%D = A(1:min_row, 1:min_row);
%G = A(min_row+1:end, min_row+1:end);

