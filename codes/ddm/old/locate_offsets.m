function [offsets, sz] = locate_offsets(size_tree, level, node, offsets, sz)
%This function locates the the row indices of the separators of 2-way
%nested dissection
%
%size_tree(input): the sizes array returned by metis for 2-way nested dissection
%level(input): number of levels on input
%node(input): the root node
%sz(dummy) is a dummy output useful for recursive implementation
%offsets(output): array containing start and end row indices of separators

if (level == 1)
    sz = size_tree{level}{node};
else
        [offsets, sz_left] = locate_offsets(size_tree, level-1, 2*node - 1, offsets, sz);
        [offsets, sz_right] = locate_offsets(size_tree, level-1, 2*node, offsets, sz);
                
        if isempty(offsets)
            start_offset = sz_left + sz_right + 1;
            end_offset = start_offset + size_tree{level}{node} - 1;
        else
            start_offset = offsets(end) + sz_left + sz_right + 1;
            end_offset = start_offset + size_tree{level}{node} - 1;
        end
        
        offsets = [offsets start_offset end_offset];
        
end

end
