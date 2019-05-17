clc, clear 
filename = 'JTJ49';
addpath(genpath('/home/pawan/work/datasets')); % from IIIT
P = load(strcat('JTJ49', '.txt'));
A = spconvert(P);
A = A + triu(A,1)';
[m,n] = size(A);
[i,j,val] = find(A);     % find to get index & value vectors
data_dump = [i,j,val];     % save this in data_dump
fid = fopen(strcat('/home/pawan/work/datasets/', filename, '.mtx'),'w');    % open file in write mode
fprintf(fid, '%d %d %d\n', m, n, nnz(A));
fprintf( fid,'%d %d %f\n', transpose(data_dump));     % print the JTJ to the file in COO 
fclose(fid);     % close file 