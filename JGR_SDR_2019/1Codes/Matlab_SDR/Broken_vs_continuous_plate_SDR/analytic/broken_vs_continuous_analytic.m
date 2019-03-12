%----------------------------------------------------------------
% Broken Plate vs Continuous plate (Analytic solution)
% Tian 2017 Jan. @Lamont
%----------------------------------------------------------------
clc
clear all
close all
%----------------------------------------------------------------
% Initialization
%----------------------------------------------------------------
L = 200000;             %[m] total length 
N = 2001;               % number of node points
%dx = L / (N - 1);       %[m] distance between node points
x = linspace(0, L, N);  %[m] initialize the x array
%w = zeros(1, N);        %[m] analytic solution to deflection w(x)
%dw = zeros(1, N);       %[m] analytic solution to deflection dw(x)
%wn = zeros(1, N);       %[m] numerical solution to deflection

Hd = 5000;              %[m] height of the dike
Te = 5000;              %[m] effective plate thickness
g = 10;                 %[m/s^2] gravitational acceleration
rho_d = 3000;           %[kg/m^3] density of the solidified dike
rho_f = 2800;           %[kg/m^3] density of the fluid dike
delta_rho_d = rho_d - rho_f;
rho_i = 2800;           %[kg/m^3] density of the infill
rho_c = 3000;           %[kg/m^3] density of the underlying lower crust or mantle
delta_rho_c = rho_c - rho_i;

E = 7.5 * 10^10;          %[Pa] Young's modulus
mu = 0.25;              % Poisson's ratio
D = (E * Te^3) / (12 * (1 - mu^2)); %[N*m] Flexural rigidity
alpha = (4 * D / (delta_rho_c * g))^0.25;%[m] Flexural wavelength
%dV_0 = dx * g * Hd * delta_rho_d; %[N/m] volcanic line load due to denser dike

%Plot the analytic solution using function analytic_plot()
brok_or_cont = 1;
analytic_plot(Hd, delta_rho_c, delta_rho_d, x, alpha, L, brok_or_cont);
hold on;
brok_or_cont = 2;
analytic_plot(Hd, delta_rho_c, delta_rho_d, x, alpha, L, brok_or_cont);

title({'Analytical solutions to eq3 and eq4 of Buck 2017 EPSL',...
    'broken plate (red) vs continuous plate (green)'}...
        ,'Fontsize',26)
grid on;
axis([0 200 -7000 1200])
axis tight

x0=10;
y0=10;
width=1050;
height=360
set(gcf,'position',[x0,y0,width,height])
                      
