% function misfit plot
%plotting  FLAC result versus semi-analytic result
function []=SDR_misfit_plot(particle_x,particle_y,particle_y_semi_analytic,i_time,num_trial,time,Te_guess,misfit)
%subplot(length_TIME,max_num_trial,i_time*num_trial)
particle_x = particle_x/1000;
particle_y = particle_y/1000;
particle_y_semi_analytic = particle_y_semi_analytic/1000;
fig = figure('visible','off');
% For setting the size of the plot
x0=10;
y0=10;
width=1600;
height=400
set(gcf,'units','points','position',[x0,y0,width,height])

plot(particle_x,particle_y,'ko')   %plot result from FLAC
hold on;
plot(particle_x,particle_y_semi_analytic,'m+') % plot semi-analytic result
axis([0 200 -30 1]); 
%title([''  'time=' num2str((time-1000)) 'kyrs; '...
%title([''  'time=' num2str((time-2000)) 'kyrs; '...
%title([''  'time=' num2str((time-3000)) 'kyrs; '...
title([''  'time=' num2str((time)) 'kyrs; '...
    'Te=' num2str(Te_guess) '; ' 'misfit=' num2str(misfit) ';'],'Fontsize',20');
xlabel('distance from the axis [km]','Fontsize',22');
ylabel('depth [km]','Fontsize',22');
set(gca,'Fontsize',22','Linewidth',2)
legend('FLAC','semi-analytic','location','southeast');
grid on;

%get the directory of your input files:
dirname = 'PDF_6_SDR3_20170502_dTe5_alpha2_tolerance20000_particle_4';
if exist(dirname) ~= 7
    mkdir(pwd,dirname);
end
%use that when you save
%fig = fullfile(pathname, strcat(num2str([i_time num_trial]),'.pdf'))
%saveas(figfile, ...');
saveas(fig,[pwd strcat('/',dirname,'/','T',num2str(i_time),'_','Trial',num2str(num_trial),'.epsc')])
%ref:https://www.mathworks.com/help/matlab/ref/saveas.html
%saveas(fig,...')