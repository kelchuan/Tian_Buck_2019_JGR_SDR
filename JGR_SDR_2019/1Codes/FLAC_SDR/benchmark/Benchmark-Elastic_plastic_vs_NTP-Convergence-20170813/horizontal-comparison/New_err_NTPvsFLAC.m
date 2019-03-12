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

%plot particles, how surface move with time
size_particle = size(particles_1);

n = size_particle(3);

x0=10;
y0=10;
width=1600;
height=400
set(gcf,'units','points','position',[x0,y0,width,height])
w_flac = particles_1 * 1000;   %convert from km to meters

i_start_particles = 1;
dN_sampling = 1; 
for i = n:n
%    for i = n:n
data_x = particles_1((i_start_particles:dN_sampling:size_particle(1)),1,i);
data_y = w_flac((i_start_particles:dN_sampling:size_particle(1)),2,i);
if i==n
     SDR_1 = plot(data_x,data_y,'g-.', 'Linewidth', 5);
else
    SDR_1 = plot(data_x,data_y,'b.');
    %hold on
end
    axis([0 100 -80 20]); 
    %axis([105 112 -0.4 0.01]); 
    title(['' num2str(i*2) 'kyrs'],'Fontsize',26');
    pause(0.05);
    xlabel('distance from the axis [km]','Fontsize',26');
    ylabel('depth [km]','Fontsize',26');
    set(gca,'Fontsize',26','Linewidth',3)
    hold on;
    
end

%----------------------------------------------------------------
% Initialization
%----------------------------------------------------------------
scale_factor_dV0 = 10 % this decrease the load due to resolution change
N_SDR = 30;    %total number of SDRs
Length = 200000;             %[m] total length 
N = 2001;               % number of node points
dx = Length / (N - 1);       %[m] distance between node points
x = linspace(0, Length, N);  %[m] initialize the x array
w = zeros(1, N);        %[m] analytic solution to deflection w(x)
dw = zeros(1, N);       %[m] analytic solution to deflection dw(x)
wn = zeros(1, N);       %[m] numerical solution to deflection

Hd = 6000;              %[m] height of the dike
Te = 6000;              %[m] effective plate thickness
g = 10;                 %[m/s^2] gravitational acceleration
rho_d = 3000;           %[kg/m^3] density of the solidified dike
rho_f = 2800;           %[kg/m^3] density of the fluid dike
delta_rho_d = rho_d - rho_f;
rho_i_sedi = 2300;
rho_i = 0;           %[kg/m^3] density of the infill
rho_c = 3000;           %[kg/m^3] density of the underlying lower crust or mantle
sediment = 0;  % whether infill is sedi or lava  1 means sedi
if sediment == 1
    delta_rho_c = rho_c - rho_i_sedi;
else
    delta_rho_c = rho_c - rho_i;
end

E = 7.5 * 10^10;          %[Pa] Young's modulus
%E = 10 * 10^10;          %[Pa] Young's modulus
mu = 0.25;              % Poisson's ratio
D = (E * Te^3) / (12 * (1 - mu^2)); %[N*m] Flexural rigidity
alpha = (4 * D / (delta_rho_c * g))^0.25;%[m] Flexural wavelength

dV_0 = dx * g * Hd * delta_rho_d * scale_factor_dV0; %[N/m] volcanic line load due to denser dike


%Plot the analytic solution using function analytic_plot()
%wa=analytic_plot(Hd, delta_rho_c, delta_rho_d, data_x, alpha,dV_0,g);
dw_0 = 2 * dV_0 / (alpha * g * delta_rho_c);
wa = -dw_0 * exp(-data_x*1000/alpha) .* cos(data_x*1000/alpha); % (eq.1)

%----------------------------------------------------------------
% 2. Use polyfit to approximate the semi-analytic result
%----------------------------------------------------------------
plot(data_x,wa,'r-.','linewidth',3)
legend('FLAC', 'ATP');
grid on;

residual_L2 = sum((data_y(1:length(data_x))-wa(1:length(data_x))).^2)
L2_norm = sqrt(residual_L2)
L2_norm_avg = L2_norm / length(data_x)


residual_L1 = sum(abs(data_y(1:length(data_x))-wa(1:length(data_x))));
L1_norm = residual_L1
L1_norm_avg = L1_norm / length(data_x)