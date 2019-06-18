clc;
filepath = '~/Test/' 
addpath(strcat(filepath,'49'))
filename = 'JTJ49_2';
A = load(strcat(filename,'.txt'));

lhs = spconvert(A);
%lhs = lhs + triu(lhs,1)';
%size(A)
[m,n] = size(lhs)
%save mat file
save(strcat(strcat(filepath,'49/'),strcat(filename,'.mat')),'lhs');


%rhs_filename = strcat('JTe',extractAfter(filename,'JTJ'))
rhs_filename = 'JTe49_2'
rhs = load(strcat(rhs_filename,'.txt'));
save(strcat(strcat(filepath,'49/'),strcat(rhs_filename,'.mat')),'rhs');