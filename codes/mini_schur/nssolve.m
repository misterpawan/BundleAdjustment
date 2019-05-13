function [x] = nssolve(nD, nG, E, F, LD, UD, LG, UG, y)
    y1 = y(1:nD); y2 = y(nD+1:nD+nG);  
    z1 = UD\(LD\y1); z2 = UG\(LG\(y2 - F*z1)); 
    x2 = z2; x1 = z1 - UD\LD\(E * x2);
    x = [x1; x2];
end