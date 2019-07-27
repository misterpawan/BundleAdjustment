clc;
matfile = '~/Test/1778/JTJ1778_1.txt';

% P = load(matfile); A = P.lhs;
A = load(matfile);
A = spconvert(A);

[m,n] = size(A)
keyboard;
p = 1778*9;

min_row = m;
% 
% for i = 4:n
%     if(abs(A(1,i)) > 0 && i < min_row)
%         min_row = i;
%         break;
%     end
% end

% for i = 1 :(m-p)
%     i
%     %find(abs(A(i,i+3:n))>0,1)
%    for j = i+3 : min_row
%        if(abs(A(i,j)) > 0)
%            if (j < min_row)
%                min_row = j;
%                break;
%            else
%                break;
%            end
%        end
%    end
% end

for i= 1 :(m-p)
     i
   if(A(i,i+3) ~= 0) 
       break;
   end
   S = find(abs(A(i,:))>0);
%     min_tmp = min(S)
%     A(i,S)
%     find(A(i,S) ~= 0)
%     R = (min_tmp-1) + find(A(i,S));
   idx = S > (i+2);
   %S(idx(1))
   min_tmp = min(S(idx));
   min_row =  min(min_row, min_tmp)
  
end

min_row 
D = A(1:min_row, 1:min_row);
G = A(min_row+1:end, min_row+1:end);

% for 372 dataset, sizeG = 3399
% for 539 dataset, sizeG = 5739
% for 885 dataset, sizeG = 12948
