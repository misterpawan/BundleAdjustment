% To plot the memory requirements

%for 30 MSC blocks
% b_jacobi_mem = [19217160,49199616,117670344,141427968,162024696,242722080,381057120,442847160,605884608,1030508640];
% 
% msc_mem=[19255424,49375312,118528576,141719208,163715552,248045072,381651232,443591544,606821072,1033169544];


b_jacobi_mem = [7606080,19480320,46597440,56006400,64163520,96122880,188172480,218686080,299198400,508889280];

msc_mem=[7644380,19656320,47455440,56297640,65854376,101445872,188766592,219430464,300134864,511550184];


direct_mem = [7763760,23448024,72165216,53012496,106024536,182655504,253902216,308714520,420342792,749947512];

num_col = [23769,60876,145617,175020,200511,300384,588039,688394,934995,1590279];

%semilogy(b_jacobi_mem,num_col);
plot(num_col,b_jacobi_mem,'r-^',num_col,msc_mem,'b-^',num_col,direct_mem,'g-^');
%hold on;
title('Memory requirement for direct solve, and block jacobi and MSC preconditioned GMRES');
xlabel('No. of columns');
ylabel('Size(in bytes)');
legend('Block Jacobi','MSC(30)','Direct');


% G_mem = [31800,89400,794000,819720,5764536,22103424,2930688,2577312,1209384,4115448];
% MSC_30_mem = [38300,176000,858000,291240,1690856,5322992,594112,744384,936464,2660904];
% num_col = [23769,60876,145617,175020,200511,300384,588039,688394,934995,1590279];
% 
% 
%  plot(num_col,G_mem,'r-o',num_col,MSC_30_mem,'b-^');
% % %hold on;
% title('Memory requirement for G block and MSC(30) block');
%  xlabel('No. of columns');
%  ylabel('Size(in bytes)');
%  legend('Block G','MSC(30)');



% G_iters = [46,19,76,89,75,86,60,82,52,92];
% MSC_30_lower_iters = [39,32,55,87,38,62,40,37,43,36];
% MSC_30_full_iters = [55,56,88,86,49,94,68,77,80,93];
% num_col = [23769,60876,145617,175020,200511,300384,588039,688394,934995,1590279];
% 
% 
%  plot(num_col,G_iters,'r-o',num_col,MSC_30_lower_iters,'b-^',num_col,MSC_30_full_iters,'g-x');
% % %hold on;
% title('#GMRES iters per LM iteration');
%  xlabel('No. of columns');
%  ylabel('#GMRES iterations');
%  legend('Block G','MSC\_half(30)','MSC\_full(30)');

% direct_time = [0.11,0.44,2.94,0.9,5.39,11.45,8.17,10.97,16.69,49.11];
% G_time = [0.58,1.12,6.8,10.43,9.52,25.71,29.38,41.46,43.12,131.65];
% MSC_30_lower_time = [0.63,1.59,5.86,9.49,6.66,17,26.19,30.6,46.63,87.16];
% MSC_30_full_time = [1.18,3.29,13.07,15.46,13.12,34.61,53.52,69.41,106,250.52];
% num_col = [23769,60876,145617,175020,200511,300384,588039,688394,934995,1590279];
% 
% 
%  plot(num_col,direct_time,'m-[]',num_col,G_time,'r-^',num_col,MSC_30_lower_time,'b-^',num_col,MSC_30_full_time,'g-^');
% % %hold on;
% title('Time(in secs) per LM iteration');
%  xlabel('No. of columns');
%  ylabel('Time(in secs)');
%  legend('Direct(LDL)','Block Jacobi','MSC\_half(30)','MSC\_full(30)');


% direct_error = [0.7337,0.8938,0.7043,0.8754,0.7444,0.693,0.7132,0.8028,0.7523,0.9951];
% G_error = [0.9661,0.8289,0.7311,0.7963,0,0.6932,0.7113,0.7706,0.7699,0];
% MSC_30_lower_error = [0,0.885,0.7421,0.6773,0.7676,0.7118,0.7202,0.798,1.12,0.8527];
% MSC_30_full_error = [0.8462,0.8938,0.7057,0.6788,0.73,0.7143,0.7178,0.7893,0.755,0.9566];
% num_col = [23769,60876,145617,175020,200511,300384,588039,688394,934995,1590279];
% 
% 
%  plot(num_col,direct_error,'m-^',num_col,G_error,'r-^',num_col,MSC_30_lower_error,'b-^',num_col,MSC_30_full_error,'g-^');
% % %hold on;
% title('Mean Reprojection Error for Lifted Schur Cost Function(Zero errors are erroneous)');
%  xlabel('No. of columns');
%  ylabel('Error');
%  legend('Direct(LDL)','Block Jacobi','MSC\_half(30)','MSC\_full(30)');
