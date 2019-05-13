function [sizes_parts, B] = domain_decomposition(A, nparts)
    %p = size(A,1):-1:1;
    B = A;
    end_index = size(B,1);
    for i=size(B,1):-1:1
        if find(B(i,1:i),1) == find(B(1:i,i),1) && find(B(1:i,i),1) == i 
            end_index = i;
            break;
        end
    end
    for i=end_index:size(B,2)    
        if nnz(B(i,1:end_index-1)) ~= 0
            end_index = i-1;
            break;
        end
    end
    sizes_parts = [];
    ind = floor(end_index/nparts);
    prev = 1;
    for i=1:nparts-1
        while ~(find(B(ind,1:ind),1) == find(B(1:ind,ind),1) && find(B(1:ind,ind),1) == ind)
            ind = ind + 1;
        end
        sizes_parts = [sizes_parts, (ind-prev)];
        prev = ind;
        ind = ind + floor(end_index/nparts);
    end
    sizes_parts = [sizes_parts, (end_index-prev+1), (size(B,1) - end_index)];
end