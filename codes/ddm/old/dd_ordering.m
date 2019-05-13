function [perm, sizes] = dd_ordering(matrix, nparts)

[perm, sizes] = PartitionAndOrdUNF_new(matrix, num2str(nparts)); 

nlev = floor(log2(nparts))+1;
sizes_tree = get_tree_sizes(sizes, nlev);
node = 1; offsets = []; sz = 0;
[offsets, sz] = locate_offsets(sizes_tree, nlev, node, offsets, sz);

perm = perm'; rows = [];

for i = 1 : 2: length(offsets)
    rows = [rows offsets(i):offsets(i+1)];  
end
perm = [ perm perm(rows) ];

perm(rows) = []; perm = perm';

end

