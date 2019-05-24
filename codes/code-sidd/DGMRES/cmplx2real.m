 function X = cmplx3real(d, Y, indx, nev)
%% builds a real eigenbasis from a complex one
%% indx is a sorting array -- select first
%% nev eigenvectors in order given by indx.
 nev = length(d);
 i = 1 ; 
 X = [] ;
 while (i<=nev) 
      j = indx(i);
      if (imag(d(j))==0) 
         X = [X, Y(:,j)] ;
         i = i+1;
      else 
	 X = [X, real(Y(:,j)),imag(Y(:,j+1))] ;
        i = i+2;
      end
end
