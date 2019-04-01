% Function to do coarse grid correction

function [t] = cogs(A_c,P,b_c)
    
    y = pcg(A_c,b_c);
    t = P*y;
    
end 
