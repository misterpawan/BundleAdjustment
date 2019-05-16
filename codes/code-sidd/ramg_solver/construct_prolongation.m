clear all
% Generate Prolongation Matrix from parts.txt

fileID = fopen('parts.txt','r');  % loading parts.txt generated from METIS graph partitioner
formatSpec = '%d';                % type of values in the file - parts.txt
p = fscanf(fileID,formatSpec);    % p is a column vector containing data from file
fclose(fileID);                   % close the file

% number of parts of aggregation : nc
nc = max(p)+1;
fprintf('No of parts of aggregation is %d\n ', nc);

% number of rows of prolongation matrix
[n,~] = size(p);
fprintf('Number of Rows of Prolongation matrix is %d\n', n);

% Construct Prolonagtion Matrix

P = sparse(n,nc);                 % generates an n-by-nc all zero sparse matrix
p = bsxfun(@plus,p,1);            % indexing 1
P(1:n, p(1:n)) = 1;               % Put only single 1 in each row
[a,b] = size(P);
fprintf('Size of Prolongation Matrix is %d by %d\n', a,b);
fprintf('Number of Non Zero elements in Prolongation Matrix %d\n:', nnz(P));










