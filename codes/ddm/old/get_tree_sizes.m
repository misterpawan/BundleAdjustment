function sizes_tree = get_tree_sizes(sizes, nlev)
% make tree data structure for sizes
% sizes: the sizes array returned from metis that stores the sizes of nodes
% ans separators
% level: level of the ND tree
% node: node number at the current level
count = 0;
for j = 1 : nlev 
    for i = 1 : 2^(nlev-j)
        count = count+1;
        sizes_tree{j}{i} = sizes(count);
    end
end
end
