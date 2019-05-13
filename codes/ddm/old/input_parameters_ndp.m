%% input parameters
%% for gmg
%matrix   = '/media/Brahmaputra/Matrices/mat_j4/matvf3dCS120.mtx';
method   = 'def';   % possible options: (a) agmg (b) hamg
smoother = 'ilu';   % possible options: (a) gs (b) jacobi (c) ilutp_sm
tol_s    = 0.1; % tol for ilutp when ilutp is the smoother
tol_c    = 0.001;  % tol for ilutp when ilutp is used for coarse solve
tol      = 1e-9;
%cf       = 3.0;     % coarsening factor
doagmg   = 0;
def      = 1;
ndef     = 77;
mgdef    = 0; %% multigrid plus deflation
skip     = 1; %% Apply deflation every 5 steps
doprecon = 1;
%% for laplace
%addpath /home/pawan/Work/Codes/ch
%nn       = 50;
%l1       = 1;
%l2       = 1;
%l3       = 1;
%% for metis
wgt          = '3'; % 3: do metis ndp, for k way part use input_parameters.m file
opt_ptype    = '0'; % 0: recursive 1: kway
opt_objtype  = '0'; 
opt_ctype    = '0'; 
opt_iptype   = '0'; 
opt_rtype    = '0'; 
opt_ncuts    = '1';
opt_nseps    = '1';
opt_num      = '0';
opt_niter    = '10';
opt_seed     = '0';
opt_minconn  = '0';
opt_contig   = '1';
opt_compress = '0';
opt_ccorder  = '0';
opt_pfactor  = '0';
opt_ufactor  = '0';
opt_dbglvl   = '0';
opt_default  = '1';
