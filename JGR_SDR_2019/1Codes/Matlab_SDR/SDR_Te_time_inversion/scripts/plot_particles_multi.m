clc
clear
close all
% plot_particles
addpath('../data_from_flac')
if exist('platform.s','file') == 0
  plat = 'native';
elseif exist('platform.s','file') == 2
  plat = 'ieee-le';
end
dirname = 'SDR_particle_video2';
if exist(dirname) ~= 7
    mkdir(pwd,dirname);
end
fid_t = fopen('time.0','r',plat);

load nxnz.0
nx=nxnz(1);
ny=nxnz(2);

t = fread(fid_t,'single');
%nt=length(t);
nt=length(t);

num_particles = 360;

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
%plot particles, how surface move with time
size_particle = size(particles_1);

n = size_particle(3);
n=900;  % for SDR_3 from 68~230
%n=35
%figure
%fig = figure('visible','off');

% For setting the size of the plot
x0=10;
y0=10;
width=1600;
height=400
set(gcf,'units','points','position',[x0,y0,width,height])

for i = 1:n
%    for i = n:n
    SDR_1 = plot(particles_1(:,1,i),particles_1(:,2,i),'r.');
    hold on
    
    %SDR_2 = plot(particles_2(:,1,i),particles_2(:,2,i),'color',rand(1,3));
    SDR_2 = plot(particles_2(:,1,i),particles_2(:,2,i),'b.');
    hold on
  
    SDR_3 = plot(particles_3(:,1,i),particles_3(:,2,i),'g.');
    %hold on
 
    SDR_4 = plot(particles_4(:,1,i),particles_4(:,2,i),'k.');
    %hold on
    SDR_5 = plot(particles_5(:,1,i),particles_5(:,2,i),'m.');
    %hold on
    SDR_6 = plot(particles_6(:,1,i),particles_6(:,2,i),'c.');
     %hold on
    SDR_7 = plot(particles_7(:,1,i),particles_7(:,2,i),'y.');
     %hold on
     
    
    SDR_8 = plot(particles_8(:,1,i),particles_8(:,2,i),'r.');
    %hold on
    SDR_9 = plot(particles_9(:,1,i),particles_9(:,2,i),'b.');
    %hold on
    SDR_10 = plot(particles_10(:,1,i),particles_10(:,2,i),'g.');
    %hold on
     SDR_11 = plot(particles_11(:,1,i),particles_11(:,2,i),'r.');
    %hold on
     SDR_12 = plot(particles_12(:,1,i),particles_12(:,2,i),'b.');
    %hold on
%    SDR_10 = plot(particles_10(:,1,i),particles_10(:,2,i),'k*');
    %}
%    axis([0 80 -25*10^3 1000]);
    axis([0 140 -30 1]); 
    %axis([.4 10.4 -2 0]); 
    title(['' num2str(i*30) 'kyrs'],'Fontsize',26');
    pause(0.05);
    xlabel('distance from the axis [km]','Fontsize',26');
    ylabel('depth [km]','Fontsize',26');
    set(gca,'Fontsize',26','Linewidth',3)
    hold on;
    
    
    if (i<(280))
        
        %delete(SDR_1);
     
        %delete(SDR_2);
         %{
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
        %}
    end
    
%    delete(SDR_5);

% plot the boundary between SDRs and dike

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
% for producing a video
%saveas(fig,[pwd strcat('/',dirname,'/','T',num2str(i),'.pdf')])
end






