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
w_flac = particles_1 * 1000;   %convert from km to meters

i_start_particles = 1;
dN_sampling = 1; 
for i = 1:n
%    for i = n:n
data_x = particles_1((i_start_particles:dN_sampling:size_particle(1)),1,i);
data_y = w_flac((i_start_particles:dN_sampling:size_particle(1)),2,i);
if i==n
     SDR_1 = plot(data_x,data_y,'g-.', 'Linewidth', 5);
else
    SDR_1 = plot(data_x,data_y,'b.');
    %hold on
end
    axis([0 100 -1500 100]); 
    %axis([105 112 -0.4 0.01]); 
    title(['' num2str(i*2) 'kyrs'],'Fontsize',26');
    pause(0.05);
    xlabel('distance from the axis [km]','Fontsize',26');
    ylabel('depth [km]','Fontsize',26');
    set(gca,'Fontsize',26','Linewidth',3)
    hold on;
    
end

% plot the boundary between SDRs and dike

    %SDR_boundary_1 = plot(particles_1(1:3,1,i),particles_1(1:3,2,i),'r.');
    %hold on

%daspect([1 1 1])

%saveas(fig,[pwd strcat('/',dirname,'/','T',num2str(i),'.pdf')])











%----------------------------------------------------------------
% Iterative Implicit FD with updated distributed load q(x)
% Assuming lava fills in wherever wn flex beneath sea level (initial zero)
% q(x) = (rho_i - rho_c) * g * wn(x)
% when wn(x) >= 0, rho_i = 0; (air density)
% when wn(x) < 0, rho_i = 2800; [kg/m^3] (lava density)
% Tian 2016 Nov. @Lamont
%----------------------------------------------------------------
%clc
%clear all
%close all
%----------------------------------------------------------------
% Initialization
%----------------------------------------------------------------
scale_factor_dV0 = 40 % this decrease the load due to resolution change
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
rho_i = 2800;           %[kg/m^3] density of the infill
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
%----------------------------------------------------------------
% Finite difference solution to D * d4w/dx4 + q(x) = 0
% Implicit solution L * W  = R
% L are the corresponding coefficients of w(i-2 ~ i+2)
% W are all the w
% R are the right hand side which are the extra loads
%----------------------------------------------------------------

%----------------------------------------------------------------
% Lava sea or not
%----------------------------------------------------------------
index_lava = 0;
% index_lava = 1 means lava sea where:
                        % q(x) = (rho_i - rho_c) * g * w(x) for any w(x)
% index_lava = 0 means:
                        % q(x) = (rho_i - rho_c) * g at w(x) < 0
                        % q(x) = - rho_c * g at w(x) >= 0
%----------------------------------------------------------------

%----------------------------------------------------------------
% Setup the R matrix (N-5+1+4, 1)  % THe added 4 is for 4 BCs
%----------------------------------------------------------------
R = zeros(N-5+1+4, 1);


%----------------------------------------------------------------
% setup the L matrix (coefficients of w)
% D/dx^4 * (w(i+2) - 4*w(i+1) + 6*w(i) - 4*w(i-1) + w(i-2)) + ...
% delta_rho_c * g * w(i) = 0
%----------------------------------------------------------------
num_i = N;  
num_j = N - 4;
L = zeros(N, N);

%------------------------------------------
% Boundary conditions
%------------------------------------------
index_broken = 1; % if index_broken = 1, then it is calculating broken plate
                  % if index_broken != 1, it is for continuous plate

% BC1 approximate line load due to denser dike at the center 
% R(1) = -dV_0 / 2 / (D / dx^4);
% Rather than applying on the RHS of the first equation, we use BCs of 
% V = dV_0 at x = 0

if (index_broken == 1)
    % The remaining four BCs can be added to the L and R matrices
    % The row number of the four BCs is random and interchangeable
    %---------------------
    % BC1: wn(N) = 0;
    %---------------------
    L(num_j+1,num_i) = 1;
    %---------------------
    % BC2: w'(inf) = 0;
    % wn(N-1) = wn(N);
    %---------------------
    L(num_j+2,num_i) = 1;
    L(num_j+2,num_i-1) = -1;
    %---------------------
    % BC3: M = 0 at x = 0  --> d2w/dx2|(x=0) = 0
    % d2w/dx2 = (wn(i+1) - 2wn(i) + wn(i-1)) / (dx)^2;
    % M = -D * d2w/dx2 = -D * (wn(i+1) - 2wn(i) + wn(i-1)) / (dx)^2;
    %---------------------
    L(num_j+3, 1) = 1;
    L(num_j+3, 2) = -2;
    L(num_j+3, 3) = 1;
    %---------------------
    % BC4 V = dM/dx = -D * d3w/dx3 = V_0/2 = dV_0 at x = 0
    % d3w/dx3|(i-0.5) = (d2w/dx2|i - d2w/dx2|i-1) / dx
    %               = (wn(i+1) -3*wn(i)+3*wn(i-1)-wn(i-2)) / (dx)^3
    % -D * (-wn(1) + 3 * wn(2) - 3 * wn(3) + wn(4)) / dx^3;
    %---------------------
    L(num_j+4, 1) = -1;
    L(num_j+4, 2) = 3;
    L(num_j+4, 3) = -3;         
    L(num_j+4, 4) = 1;
    R(num_j+4,1) = dV_0 / (-D/dx^3); % value at the RHS
    
else % Continuous plate (index_broken != 1)
    %---------------------
    %BC1
    %wn(N) = 0;
    %---------------------
    L(num_j+1,num_i) = 1;
    %---------------------
    %BC2
    %wn(N-1) = wn(N);
    %---------------------
    L(num_j+2,num_i) = 1;
    L(num_j+2,num_i-1) = -1;
    %---------------------
    % BC3 wn' = 0 at x = 0
    % wn(1) = wn(2);
    %---------------------
    L(num_j+3, 1) = 1;
    L(num_j+3, 2) = -1;
    %---------------------
    % BC4 V = 0 at x = 0
    %wn(4) = wn(1) - 3 * wn(2) + 3 * wn(3);
    %---------------------
    L(num_j+4, 1) = -1;
    L(num_j+4, 2) = 3;
    L(num_j+4, 3) = -3;
    L(num_j+4, 4) = 1;
    R(num_j+4,1) = dV_0 / (-D/dx^3); % value at the RHS
end
%figure
num_iterate = 10;
for i=1:1:num_iterate;
% update L base on wn   --> Assigning distributed load q(x)
if (index_lava == 0)
    L = distributed_q( L,N,wn,g,dx,rho_c,delta_rho_c, D); 
end

% update wn
wn = L\R;
%wn = wn; %convert from meters to km
if (i==num_iterate)
    %plot(x/1000, wn,'k.');hold on;
    %plot(x/1000, wn,'k-.','Linewidth', 5);hold on;
    plot(x/1000, wn,'k-.','Linewidth', 5);hold on;

elseif (i<=3)
    if(i ==1)
        %plot(x/1000, wn,'r','Linewidth', 1);hold on;
    elseif(i==2)
        %plot(x/1000, wn,'r','Linewidth', 1);hold on;
    elseif(i==3)
         %plot(x/1000, wn,'r','Linewidth', 1);hold on;
    end
else
     %plot(x/1000, wn,'color',rand(1,3));hold on;
     %plot(x/1000, wn,'r');hold on;
end
%legend('1st', '2nd', '3rd');
%axis([0 Length/1000/9 -400 100])
%pause(0.0001);
grid on
end


xlabel('distance from the axis [km]','Fontsize',26');
ylabel('deflection [m]','Fontsize',26');
%title('iterated steady state dw(x) with half dike load of width of dx'...
%,'Fontsize',26')
set(gca,'Fontsize',26','Linewidth',3)
%axis([0 L/1000 -8000 1000])

%----------------------------------------------------------------
% 2. Use polyfit to approximate the semi-analytic result
%----------------------------------------------------------------
x = x ./ 1000;
x = x';
npoly = 26;   % degree of the fitting polynomial

%[p,S,mu] = polyfit(X(1:dN_sampling:length(X)),Ws(1:dN_sampling:length(X)),npoly);
[p,S] = polyfit(x(1:length(x)),wn(1:length(x)),npoly);


%Y = polyval(p,x);  % get the fitting curve Y value
Y = polyval(p,x);  % get the fitting curve Y value
%plot(X(1:dN_sampling:length(X))/1000,Y(1:dN_sampling:length(X)),'m+');
norm_poly=S.normr
%legend('boundary','SDR0','SDR1','SDR2','SDR3','SDR4','SDR5','SDR6',...
%    'semi-analytic','semi-analytic(sampled)','polyfit')

%ploting the error between polyfit and semi-analytic data
%figure
%plot(X(1:dN_sampling:length(X))/1000,Ws(1:dN_sampling:length(X))-Y(1:dN_sampling:length(X)),'k+');
%ylabel('absolute error')
%xlabel('x [km]')

residual_poly_L2=sum((wn(1:length(x))-Y(1:length(x))).^2);
residual_poly_L1=sum(abs(wn(1:length(x))-Y(1:length(x))));
L2_norm_poly = sqrt(residual_poly_L2);
L2_norm_poly_avg = sqrt(residual_poly_L2)/length(wn(1:length(x)));
L1_norm_poly = residual_poly_L1;
L1_norm_poly_avg = L1_norm_poly / length(wn(1:length(x)))
%}


%YY = polyval(p,(data_x-0.5)); % get the NTP YY corresponding to the data_x location
YY = polyval(p,(data_x));
% the data_x - 0.5km means the center should be offseted to the center of
% the dike which has a width of 1 km

%plot(data_x-0.5,YY,'r-.','linewidth',5)
%legend('FLAC', 'NTP', 'NTP polyfit');
plot(data_x,YY,'r-.','linewidth',5)

residual_L2 = sum((data_y(1:length(data_x))-YY(1:length(data_x))).^2)
L2_norm = sqrt(residual_L2)
L2_norm_avg = L2_norm / length(data_x)


residual_L1 = sum(abs(data_y(1:length(data_x))-YY(1:length(data_x))));
L1_norm = residual_L1
L1_norm_avg = L1_norm / length(data_x)