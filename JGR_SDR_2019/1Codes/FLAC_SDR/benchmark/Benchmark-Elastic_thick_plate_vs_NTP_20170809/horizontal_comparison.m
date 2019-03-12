clc
clear all
close all
% plot_particles
if exist('platform.s','file') == 0
  plat = 'native';
elseif exist('platform.s','file') == 2
  plat = 'ieee-le';
end
dirname = 'SDR_particle_video';
if exist(dirname) ~= 7
    %mkdir(pwd,dirname);
end
fid_t = fopen('time.0','r',plat);

load nxnz.0
nx=nxnz(1);
ny=nxnz(2);

t = fread(fid_t,'single');
%nt=length(t);
nt=length(t);

num_particles = 600;

if exist('particles_01.0','file')>0
    fid_part = fopen('particles_01.0','r',plat);
    junk = fread(fid_part,'single'); 
%    particles = reshape(junk,length(junk)/2/nt,2,nt);
    particles_1 = reshape(junk,num_particles,2,floor(length(junk)...
        /(num_particles*2)));
end
if exist('particles_02.0','file')>0
    fid_part = fopen('particles_02.0','r',plat);
    junk = fread(fid_part,'single'); 
%    particles = reshape(junk,length(junk)/2/nt,2,nt);
    particles_2 = reshape(junk,num_particles,2,floor(length(junk)...
        /(num_particles*2)));
end
if exist('particles_03.0','file')>0
    fid_part = fopen('particles_03.0','r',plat);
    junk = fread(fid_part,'single'); 
%    particles = reshape(junk,length(junk)/2/nt,2,nt);
    particles_3 = reshape(junk,num_particles,2,floor(length(junk)...
        /(num_particles*2)));
end
if exist('particles_04.0','file')>0
    fid_part = fopen('particles_04.0','r',plat);
    junk = fread(fid_part,'single'); 
%    particles = reshape(junk,length(junk)/2/nt,2,nt);
    particles_4 = reshape(junk,num_particles,2,floor(length(junk)...
        /(num_particles*2)));
end
%num_particles = 800;
if exist('particles_05.0','file')>0
    fid_part = fopen('particles_05.0','r',plat);
    junk = fread(fid_part,'single'); 
%    particles = reshape(junk,length(junk)/2/nt,2,nt);
    particles_5 = reshape(junk,num_particles,2,floor(length(junk)...
        /(num_particles*2)));
end
if exist('particles_06.0','file')>0
    fid_part = fopen('particles_06.0','r',plat);
    junk = fread(fid_part,'single'); 
%    particles = reshape(junk,length(junk)/2/nt,2,nt);
    particles_6 = reshape(junk,num_particles,2,floor(length(junk)...
        /(num_particles*2)));
end

%plot particles, how surface move with time
size_particle = size(particles_1);

%n = size_particle(3);
%n=100;  % for SDR_3 from 68~230
%n=35
%figure
%fig = figure('visible','off');

% For setting the size of the plot
x0=10;
y0=10;
width=1600;
height=400
set(gcf,'units','points','position',[x0,y0,width,height])
w_flac_1 = particles_1 * 1000;   %convert from km to meters
w_flac_2 = particles_2 * 1000;   %convert from km to meters
w_flac_3 = particles_3 * 1000;   %convert from km to meters
w_flac_4 = particles_4 * 1000;   %convert from km to meters
w_flac_5 = particles_5 * 1000;   %convert from km to meters
%w_flac_6 = particles_6 * 1000;   %convert from km to meters

i_start_particles = 1;
dN_sampling = 1; 
i = 50;
data_x1 = particles_1((i_start_particles:dN_sampling:size_particle(1)),1,i);
data_y1 = w_flac_1((i_start_particles:dN_sampling:size_particle(1)),2,i);

data_x2 = particles_2((i_start_particles:dN_sampling:size_particle(1)),1,i);
data_y2 = w_flac_2((i_start_particles:dN_sampling:size_particle(1)),2,i);

data_x3 = particles_3((i_start_particles:dN_sampling:size_particle(1)),1,i);
data_y3 = w_flac_3((i_start_particles:dN_sampling:size_particle(1)),2,i);

data_x4 = particles_4((i_start_particles:dN_sampling:size_particle(1)),1,i);
data_y4 = w_flac_4((i_start_particles:dN_sampling:size_particle(1)),2,i);

size_particle = size(particles_5);
data_x5 = particles_5((i_start_particles:dN_sampling:size_particle(1)),1,i);
data_y5 = w_flac_5((i_start_particles:dN_sampling:size_particle(1)),2,i);

%data_x6 = particles_6((i_start_particles:dN_sampling:size_particle(1)),1,i);
%data_y6 = w_flac_6((i_start_particles:dN_sampling:size_particle(1)),2,i);

SDR_1 = plot(data_x1,data_y1,'r-.', 'Linewidth', 3); hold on;
SDR_2 = plot(data_x2,data_y2,'y-.', 'Linewidth', 3); hold on;
SDR_3 = plot(data_x3,data_y3,'g-.', 'Linewidth', 3); hold on;
SDR_4 = plot(data_x4,data_y4,'m-.', 'Linewidth', 3); hold on;
SDR_5 = plot(data_x5,data_y5,'b-.', 'Linewidth', 3); hold on;
%SDR_6 = plot(data_x6,data_y6,'k-.', 'Linewidth', 2); hold on;

%legend('SDR1','SDR2','SDR3','SDR4','SDR5')
legend('3 km','2 km','1 km','0.5 km','0.25 km')
    axis([0 100 -400 20]); 
    grid on
    xlabel('distance from the axis [km]','Fontsize',26');
    ylabel('depth [km]','Fontsize',26');
    set(gca,'Fontsize',26','Linewidth',3)
    

    
%---------------------------------------------------------------------

npoly = 26;   % degree of the fitting polynomial

%[p,S,mu] = polyfit(X(1:dN_sampling:length(X)),Ws(1:dN_sampling:length(X)),npoly);
[p,S] = polyfit(data_x5,data_y5,npoly);

Y = polyval(p,data_x5);  % get the fitting curve Y value
figure
plot(data_x5,Y,'m+'); hold on;
plot(data_x5,data_y5,'k-.', 'Linewidth', 2); hold on;
axis([0 100 -400 20]); 
grid on
xlabel('distance from the axis [km]','Fontsize',26');
ylabel('depth [km]','Fontsize',26');
set(gca,'Fontsize',26','Linewidth',3)
norm_poly=S.normr

Y1 = polyval(p,data_x1);  % get the fitting curve Y value
Y2 = polyval(p,data_x2);  % get the fitting curve Y value
Y3 = polyval(p,data_x3);  % get the fitting curve Y value
Y4 = polyval(p,data_x4);  % get the fitting curve Y value
Y5 = polyval(p,data_x5);  % get the fitting curve Y value
%Y6 = polyval(p,data_x6);  % get the fitting curve Y value



%residual_L2 = sum((data_y(1:length(data_x))-YY(1:length(data_x))).^2)
%L2_norm = sqrt(residual_L2)
%L2_norm_avg = L2_norm / length(data_x)


residual_L1_1 = sum(abs(data_y1-Y1));
L1_norm_1 = residual_L1_1
L1_norm_avg_1 = L1_norm_1 / length(data_x1)

residual_L1_2 = sum(abs(data_y2-Y2));
L1_norm_2 = residual_L1_2
L1_norm_avg_2 = L1_norm_2 / length(data_x2)

residual_L1_3 = sum(abs(data_y3-Y3));
L1_norm_3 = residual_L1_3
L1_norm_avg_3 = L1_norm_3 / length(data_x3)

residual_L1_4 = sum(abs(data_y4-Y4));
L1_norm_4 = residual_L1_4
L1_norm_avg_4 = L1_norm_4 / length(data_x4)

residual_L1_5 = sum(abs(data_y5-Y5));
L1_norm_5 = residual_L1_5
L1_norm_avg_5 = L1_norm_5 / length(data_x5)

%residual_L1_6 = sum(abs(data_y6-Y6));
%L1_norm_6 = residual_L1_6
%L1_norm_avg_6 = L1_norm_6 / length(data_x6)

reso = [3000 2000 1000 500 250];
L1_norm_avg = [L1_norm_avg_1 L1_norm_avg_2 L1_norm_avg_3 L1_norm_avg_4 ...
    L1_norm_avg_5];

figure
npoly2 = 3;
[p2,S2] = polyfit(reso,L1_norm_avg,npoly2);
reso_x = linspace(min(reso)-100,max(reso)+100,100);
Yp2 = polyval(p2,reso_x);
plot(reso_x/1000, Yp2, 'r-.', 'Linewidth', 3)
S2.normr
hold on;

plot(reso/1000, L1_norm_avg,'b.', 'Markersize', 68); hold on;
grid on
xlabel("grid size [km]")
ylabel("error in L1 norm on average [m]")
set(gca,'Fontsize',26','Linewidth',3)
axis([0 3.2 -0.8 12])




figure
plot(data_x1, Y1); hold on;
plot(data_x2, Y2); hold on;
plot(data_x3, Y3); hold on;
plot(data_x4, Y4); hold on;
plot(data_x5, Y5); hold on;
%plot(data_x6, Y6); hold on;
%}