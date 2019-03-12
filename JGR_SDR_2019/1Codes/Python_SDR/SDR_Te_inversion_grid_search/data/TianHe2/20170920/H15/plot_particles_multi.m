clc
clear
close all
% plot_particles
% move to the data directory
%dirname = ['/Users/Tian/Documents/MATLAB/SDR/SDR_Te_Inversion/' ...
%                 'constant_thermal_Te/data/Roche/Te15_w0_20_1km_new']
%dirname = ['/Users/Tian/Documents/MATLAB/SDR/SDR_Te_Inversion/'...
%    'constant_thermal_Te/data/TH2/20170514/H6']
%cd(dirname)

if exist('platform.s','file') == 0
  plat = 'native';
elseif exist('platform.s','file') == 2
  plat = 'ieee-le';
end
%{
dirname = 'SDR_particle_video2';
if exist(dirname) ~= 7
    mkdir(pwd,dirname);
end
%}
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

if exist('particles_07.0','file')>0
    fid_part = fopen('particles_07.0','r',plat);
    junk = fread(fid_part,'single'); 
%    particles = reshape(junk,length(junk)/2/nt,2,nt);
    particles_7 = reshape(junk,num_particles,2,floor(length(junk)...
        /(num_particles*2)));
end
if exist('particles_08.0','file')>0
    fid_part = fopen('particles_08.0','r',plat);
    junk = fread(fid_part,'single'); 
%    particles = reshape(junk,length(junk)/2/nt,2,nt);
    particles_8 = reshape(junk,num_particles,2,floor(length(junk)...
        /(num_particles*2)));
end

if exist('particles_09.0','file')>0
    fid_part = fopen('particles_09.0','r',plat);
    junk = fread(fid_part,'single'); 
%    particles = reshape(junk,length(junk)/2/nt,2,nt);
    particles_9 = reshape(junk,num_particles,2,floor(length(junk)...
        /(num_particles*2)));
end

if exist('particles_10.0','file')>0
    fid_part = fopen('particles_10.0','r',plat);
    junk = fread(fid_part,'single'); 
%    particles = reshape(junk,length(junk)/2/nt,2,nt);
    particles_10= reshape(junk,num_particles,2,floor(length(junk)...
        /(num_particles*2)));
end
if exist('particles_11.0','file')>0
    fid_part = fopen('particles_11.0','r',plat);
    junk = fread(fid_part,'single'); 
%    particles = reshape(junk,length(junk)/2/nt,2,nt);
    particles_11= reshape(junk,num_particles,2,floor(length(junk)...
        /(num_particles*2)));
end
if exist('particles_12.0','file')>0
    fid_part = fopen('particles_12.0','r',plat);
    junk = fread(fid_part,'single'); 
%    particles = reshape(junk,length(junk)/2/nt,2,nt);
    particles_12= reshape(junk,num_particles,2,floor(length(junk)...
        /(num_particles*2)));
end
if exist('particles_13.0','file')>0
    fid_part = fopen('particles_13.0','r',plat);
    junk = fread(fid_part,'single'); 
%    particles = reshape(junk,length(junk)/2/nt,2,nt);
    particles_13= reshape(junk,num_particles,2,floor(length(junk)...
        /(num_particles*2)));
end
if exist('particles_14.0','file')>0
    fid_part = fopen('particles_14.0','r',plat);
    junk = fread(fid_part,'single'); 
%    particles = reshape(junk,length(junk)/2/nt,2,nt);
    particles_14= reshape(junk,num_particles,2,floor(length(junk)...
        /(num_particles*2)));
end
if exist('particles_15.0','file')>0
    fid_part = fopen('particles_15.0','r',plat);
    junk = fread(fid_part,'single'); 
%    particles = reshape(junk,length(junk)/2/nt,2,nt);
    particles_15= reshape(junk,num_particles,2,floor(length(junk)...
        /(num_particles*2)));
end
if exist('particles_16.0','file')>0
    fid_part = fopen('particles_16.0','r',plat);
    junk = fread(fid_part,'single'); 
%    particles = reshape(junk,length(junk)/2/nt,2,nt);
    particles_16= reshape(junk,num_particles,2,floor(length(junk)...
        /(num_particles*2)));
end
%plot particles, how surface move with time
size_particle = size(particles_1);

n = size_particle(3);

%n=35
%n=10000
%figure
%fig = figure('visible','off');

% For setting the size of the plot
x0=10;
y0=10;
width=2000;
height=400
set(gcf,'units','points','position',[x0,y0,width,height])


tip_x = zeros(1,16)
tip_y = zeros(1,16)
%{
for j = 1:16
    %filename = sprintf('%s_%d','particles',j)
    tip_x(j) = sprintf('%s_%d(1,1,n)','particles',j)
    tip_y(j) = sprintf('%s_%d(1,2,n)','particles',j)
end
%}
tip_x(1) = particles_16(1,1,n)
tip_x(2) = particles_15(1,1,n)
tip_x(3) = particles_14(1,1,n)
tip_x(4) = particles_13(1,1,n)
tip_x(5) = particles_12(1,1,n)
tip_x(6) = particles_11(1,1,n)
tip_x(7) = particles_10(1,1,n)
tip_x(8) = particles_9(1,1,n)
tip_x(9) = particles_8(1,1,n)
tip_x(10) = particles_7(1,1,n)
tip_x(11) = particles_6(1,1,n)
tip_x(12) = particles_5(1,1,n)
tip_x(13) = particles_4(1,1,n)
tip_x(14) = particles_3(1,1,n)
tip_x(15) = particles_2(1,1,n)
tip_x(16) = particles_1(1,1,n)


tip_y(1) = particles_16(1,2,n)
tip_y(2) = particles_15(1,2,n)
tip_y(3) = particles_14(1,2,n)
tip_y(4) = particles_13(1,2,n)
tip_y(5) = particles_12(1,2,n)
tip_y(6) = particles_11(1,2,n)
tip_y(7) = particles_10(1,2,n)
tip_y(8) = particles_9(1,2,n)
tip_y(9) = particles_8(1,2,n)
tip_y(10) = particles_7(1,2,n)
tip_y(11) = particles_6(1,2,n)
tip_y(12) = particles_5(1,2,n)
tip_y(13) = particles_4(1,2,n)
tip_y(14) = particles_3(1,2,n)
tip_y(15) = particles_2(1,2,n)
tip_y(16) = particles_1(1,2,n)

dif_Xf = diff(tip_y)./diff(tip_x)

%for i = 1:n
    for i = n:n
       
    SDR_1 = plot(particles_1(:,1,i),particles_1(:,2,i),'r.');
    hold on
    
    %SDR_2 = plot(particles_2(:,1,i),particles_2(:,2,i),'color',rand(1,3));
    SDR_2 = plot(particles_2(:,1,i),particles_2(:,2,i),'b.');
    hold on
  
    SDR_3 = plot(particles_3(:,1,i),particles_3(:,2,i),'g.');
    hold on
 
    SDR_4 = plot(particles_4(:,1,i),particles_4(:,2,i),'k.');
    hold on
    SDR_5 = plot(particles_5(:,1,i),particles_5(:,2,i),'m.');
    hold on
    SDR_6 = plot(particles_6(:,1,i),particles_6(:,2,i),'c.');
    hold on
    
    
    SDR_7 = plot(particles_7(:,1,i),particles_7(:,2,i),'y.');
    hold on
     
    
    SDR_8 = plot(particles_8(:,1,i),particles_8(:,2,i),'r.');
    hold on
    
    SDR_9 = plot(particles_9(:,1,i),particles_9(:,2,i),'b.');
    hold on
    SDR_10 = plot(particles_10(:,1,i),particles_10(:,2,i),'g.');
    hold on
    SDR_11 = plot(particles_11(:,1,i),particles_11(:,2,i),'k.');
    hold on
    SDR_12 = plot(particles_12(:,1,i),particles_12(:,2,i),'m.');
    hold on
    SDR_13 = plot(particles_13(:,1,i),particles_13(:,2,i),'c.');
    hold on
    SDR_14 = plot(particles_14(:,1,i),particles_14(:,2,i),'y.');
    hold on
    SDR_15 = plot(particles_15(:,1,i),particles_15(:,2,i),'r.');
    hold on
 
    SDR_16 = plot(particles_16(:,1,i),particles_16(:,2,i),'b.');
    hold on
    
    plot(tip_x,tip_y,'ro','markersize',22)
    %}
%    axis([0 80 -25*10^3 1000]);
    axis([0 200 -30 1]); 
    %axis([105 112 -0.4 0.01]); 
    title(['' num2str(i*20) 'kyrs'],'Fontsize',26');
    pause(0.02);
    xlabel('distance from the axis [km]','Fontsize',26');
    ylabel('depth [km]','Fontsize',26');
    set(gca,'Fontsize',26','Linewidth',3)
    hold on;
    grid on
    
    
    if (i<=(750))
        
        delete(SDR_1);
     
        delete(SDR_2);
         
        delete(SDR_3);
       
        delete(SDR_4);
        delete(SDR_5);
        delete(SDR_6);
        
        delete(SDR_7);
        delete(SDR_8);
        
        delete(SDR_9);
        delete(SDR_10);
        delete(SDR_11);
        delete(SDR_12);
        delete(SDR_13);
        delete(SDR_14);
        delete(SDR_15);
        delete(SDR_16);
        %}
    end


% plot the boundary between SDRs and dike
%{
    SDR_boundary_1 = plot(particles_1(1:3,1,i),particles_1(1:3,2,i),'r.');
    %hold on

    SDR_boundary_2 = plot(particles_2(1:3,1,i),particles_2(1:3,2,i),'b.');
    %hold on
    SDR_boundary_3 = plot(particles_3(1:3,1,i),particles_3(1:3,2,i),'g.');
    %hold on
    
    SDR_boundary_4 = plot(particles_4(1:3,1,i),particles_4(1:3,2,i),'k.');
    %hold on
    SDR_boundary_5 = plot(particles_5(1:3,1,i),particles_5(1:3,2,i),'m.');
    %hold on
    SDR_boundary_6 = plot(particles_6(1:3,1,i),particles_6(1:3,2,i),'c.');
     %hold on
%}
    %{
    SDR_boundary_7 = plot(particles_7(1:3,1,i),particles_7(1:3,2,i),'y.');
     %hold on
    SDR_boundary_8 = plot(particles_8(1:3,1,i),particles_8(1:3,2,i),'r.');
    %hold on
    SDR_boundary_9 = plot(particles_9(1:3,1,i),particles_9(1:3,2,i),'b.');
    %hold on
    SDR_boundary_10 = plot(particles_10(1:3,1,i),particles_10(1:3,2,i),'b.');
    %hold on
%}
daspect([1 1 1])

%saveas(fig,[pwd strcat('/',dirname,'/','T',num2str(i),'.pdf')])
    end

    %gradient_limit = 0.01
    gradient_limit = tand(1)
    figure 
    plot(dif_Xf); hold on;
    plot(dif_Xf, 'r+'); hold on;
    plot(zeros(length(dif_Xf))); hold on; 
    plot(gradient_limit * ones(length(dif_Xf))); hold on;
    plot(-gradient_limit * ones(length(dif_Xf))); hold on;
    
    for k = 1:length(dif_Xf)
        
        if abs(dif_Xf(k)) <= gradient_limit       
            plot(k,dif_Xf(k),'bo','markersize',16)
            k
            tip_x(k)
            %print("the corresponding Xf is %f",)
            xlabel('number of SDR','Fontsize',26');
            ylabel('slope between tip of SDRs','Fontsize',26');
            set(gca,'Fontsize',26','Linewidth',3)
            grid on
        end
    end
    





