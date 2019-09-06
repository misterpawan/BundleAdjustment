% To plot the memory requirements

%for 30 MSC blocks
% b_jacobi_mem = [19217160,49199616,117670344,141427968,162024696,242722080,381057120,442847160,605884608,1030508640];
% 
% msc_mem=[19255424,49375312,118528576,141719208,163715552,248045072,381651232,443591544,606821072,1033169544];


% b_jacobi_mem = [7641360,19579680,61359184,213299152,241762456,80121048,217654944,938095616];
% 
% msc_mem=[7686136,19845992,61984752,189434216,220246056,48475048,67909064,109743320];
% 
% 
% direct_mem = [7763760,23448024,53012496,253902216,308714520,72165216,106024536,182655504];
% 
% b_jacobi_mem = [31170774,7988283,12953810,23277506];
% 
% msc_mem=[29622119,6172150,8535290,13079639];
% 
% 
% direct_mem = [39262193,9166163,13453379,23131987];
% 
% 
% % % num_cam = [49,138,225,308,356,372,539,885];
% num_cam = [356,372,539,885];
% %num_col = [23769,60876,145617,200511,300384,476238,588039,688394];
% 
% %semilogy(b_jacobi_mem,num_col);
% plot(num_cam,b_jacobi_mem,'r-o',num_cam,msc_mem,'b-^',num_cam,direct_mem,'g-*');
% %hold on;
% %title('Memory requirement for direct solve, and block jacobi and MSC preconditioned GMRES');
% xlabel('Number of cameras','fontweight','bold','fontsize',16);
% ylabel('Size(in bytes)','fontweight','bold','fontsize',16);
% legend('Block Jacobi','MSC(30)','Direct','location','NorthEast','fontsize',30);
% set(gca,'fontsize',30);
% %saveas(gcf,'mem_compare_all.eps','epsc')
% saveas(gcf,'mem_compare_all.fig')


% G_mem = [35280,99360,25126672,23076376,33523608,153491424,841972736,336859072];
% MSC_30_mem = [80056,365672,1261736,1559976,1877608,3745544,13620440,15783112];
% %num_col = [23769,60876,145617,175020,200511,300384,588039,688394,934995,1590279];
% num_cam = [49,138,308,356,372,539,885,1642];
% 
%  plot(num_cam,G_mem,'r-^',num_cam,MSC_30_mem,'b-^');
% % %hold on;
% % title('Memory requirement for G block and 30 MSC blocks');
%  xlabel('Number of cameras','fontweight','bold');
%  ylabel('Size(in bytes)','fontweight','bold');
%  legend('Block G','MSC(30)');
% saveas(gcf,'~/mem_compare_G_MSC.eps','epsc')
% 
% for factor time
% G_mem = [35280,99360,25126672,23076376,33523608,153491424,841972736,336859072];
% MSC_30_mem = [80056,365672,1261736,1559976,1877608,3745544,13620440,15783112];
% num_col = [23769,60876,145617,175020,200511,300384,588039,688394,934995,1590279];
% num_cam = [49,138,225,308,356,372,539,885];
% 
%  plot(num_cam,G_mem,'r-^',num_cam,MSC_30_mem,'b-^');
% % %hold on;
% % title('Memory requirement for G block and 30 MSC blocks');
%  xlabel('Number of cameras','fontweight','bold');
%  ylabel('Size(in bytes)','fontweight','bold');
%  legend('Block G','MSC(30)');
%  set(gca,'fontsize',12);
% saveas(gcf,'mem_compare_G_MSC.eps','epsc')


% G_time = [0.018,0.048,0.215,1.59,1.98,0.135,0.403,2.848];
% MSC_30_lower_time = [0.016,0.039,0.173,1.553,1.736,0.129,0.204,0.401];
% %MSC_30_full_iters = [55,56,88,86,49,94,68,77,80,93];
% %num_col = [23769,60876,145617,175020,200511,300384,588039,688394,934995,1590279];
% num_cam = [49,138,225,308,356,372,539,885];
% 
% 
%  plot(num_cam,G_time,'r-o',num_cam,MSC_30_lower_time,'b-^');
% % %hold on;
% % title('#GMRES iters per LM iteration');
%  xlabel('Number of cameras','fontweight','bold','fontsize',16);
%  ylabel('Time(in secs)','fontweight','bold','fontsize',16);
%  legend('Block G','MSC(30)','location','NorthWest');
%  set(gca,'fontsize',12);
%  saveas(gcf,'time_compare_factor.eps','epsc')

direct_time = [11,44,90,817,1097,294,539,1145];
G_time = [58,112,1043,2938,4146,680,952,2571];
MSC_30_lower_time = [63,159,949,2619,3060,586,664,1700];
% MSC_30_full_time = [1.18,3.29,13.07,15.46,13.12,34.61,53.52,69.41,106,250.52];
% num_col = [23769,60876,145617,175020,200511,300384,588039,688394,934995,1590279];
num_cam = [49,138,225,308,356,372,539,885];

 plot(num_cam,direct_time,'g-*',num_cam,G_time,'r-o',num_cam,MSC_30_lower_time,'b-^');
% %hold on;
% title('Time(in secs) per LM iteration');
 xlabel('Number of cameras','fontweight','bold','fontsize',16);
 ylabel('Time(in secs)','fontweight','bold','fontsize',16);
 legend('Direct(LDL)','Block Jacobi','MSC(30)','location','NorthEast','fontsize',30);
 set(gca,'fontsize',30);
 %saveas(gcf,'time_compare_all.eps','epsc')
 saveas(gcf,'time_compare_all.fig')


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
