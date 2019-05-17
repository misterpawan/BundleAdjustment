filepath = '../../../Test/' 
addpath(strcat(filepath,'LadyBug'))
filename = 'JTJ356'
A = load(strcat(filename,'.txt'));

lhs = spconvert(A);
lhs = lhs + triu(lhs,1)';
%size(A)
[m,n] = size(lhs)
%save mat file
save(strcat(filepath,strcat(filename,'.mat')),'lhs');

%generating random RHS
rhs = rand(m,1);
rhs_filename = strcat('b',extractAfter(filename,'JTJ'))
save(strcat(filepath,strcat(rhs_filename,'.mat')),'rhs');