function y = recursive_Vcycle(nlev, AA, PP, r)

    t_smooth = apply_smoother(AA{nlev}, r); 
    
    tmp = r - AA{nlev} * t_smooth;
    
    if (nlev == 2) % if two grid
        t_corr = PP{nlev-1} * ( AA{nlev-1} \ (PP{nlev-1}' * tmp) );
        %y = t_smooth + t_corr; 
    else
        t_corr = PP{nlev-1} * recursive_Vcycle(nlev-1, AA, PP, PP{nlev-1}' * tmp);
        %y = t_smooth + recursive_Vcycle(nlev-1, AA, PP, PP{nlev-1}' * tmp);
    end
    y = t_smooth + t_corr; 

end