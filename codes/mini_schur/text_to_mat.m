
filepath = '../../../Test/' 
addpath(strcat(filepath,'356'))
filename = 'JTJ356_5'
A = load(strcat(filename,'.txt'));

lhs = spconvert(A);
lhs = lhs + triu(lhs,1)';
%size(A)
[m,n] = size(lhs)
%save mat file
save(strcat(strcat(filepath,'356/'),strcat(filename,'.mat')),'lhs');


%rhs_filename = strcat('JTe',extractAfter(filename,'JTJ'))
rhs_filename = 'JTe356_5'
rhs = load(strcat(rhs_filename,'.txt'));
save(strcat(strcat(filepath,'356/'),strcat(rhs_filename,'.mat')),'rhs');