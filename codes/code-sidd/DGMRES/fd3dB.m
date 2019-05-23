 function A = fd3dB(nx, ny, nz, Bx, By, Bz) 
%------------------------------
% 
% A = fd3dB(nx, ny, nz, Bx, By, Bz) 
%
% 5- or 7-point block-Diffusion/conv. matrix. 
% with constant blocks in each direction.. 
% discretizes: 
% d / dx (Bx d . /dx) +  d / dy (By d . /dy) + d / dz (Bz d . /dz)
% on a regular nx x ny x nz grid
% 
% Bx, By, Bz are constant small matrices of the same dimension.
% 
% NOTES: nx must be > 1. 
%------------------------------
%% x-direction only 
if (nx > 1) 
   tx = tridiag(2, -1, -1, nx) ;
   A  = kron(tx,Bx); 
end
% y-direction
if (ny > 1) 
  ty = tridiag(2, -1, -1, ny) ;
  A  = kron(speye(ny,ny),A) + kron(ty,kron(speye(nx,nx),By)) ;  
end
% z-direction
if (nz > 1) 
  tz = tridiag(2, -1, -1, nz) ;
  A = kron(speye(nz,nz),A)+kron(tz,kron(speye(nx*ny,nx*ny),Bz)); 
end

