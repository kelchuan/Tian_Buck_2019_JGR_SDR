% This matlab script plots the analytic solution to Buck's 2016 SDR paper 
% with updated gap and sediment infill after SDRs
% Tian 2017 June 23rd 

clc
clear all
close all

% Initialization
L = 200000;             %[m] total length 
N = 1001;               % number of node points
dx = L / (N - 1);       %[m] distance between node points
x = linspace(0, L, N);  %[m] initialize the x array
w = zeros(1, N);        %[m] analytic solution to deflection w(x)
dw = zeros(1, N);       %[m] analytic solution to deflection dw(x)
%wn = zeros(1, N);       %[m] numerical solution to deflection

% Analytic solution
Hd = 5000;              %[m] height of the dike
Te = 5000;              %[m] effective plate thickness
g = 10;                 %[m/s^2] gravitational acceleration
rho_d = 3000;           %[kg/m^3] density of the solidified dike
rho_f = 2800;           %[kg/m^3] density of the fluid dike
delta_rho_d = rho_d - rho_f;
rho_i_sdr = 2800;           %[kg/m^3] density of the SDR infill
rho_i_gap = 0;           %[kg/m^3] density of the gap(air) infill
rho_i_sedi = 2800;           %[kg/m^3] density of the sediment infill
rho_c = 3000;           %[kg/m^3] density of the underlying lower crust or mantle
delta_rho_c_sdr = rho_c - rho_i_sdr;
delta_rho_c_gap = rho_c - rho_i_gap;
delta_rho_c_sedi = rho_c - rho_i_sedi;

E = 5 * 10^10;          %[Pa] Young's modulus
mu = 0.25;              % Poisson's ratio
D = (E * Te^3) / (12 * (1 - mu^2)); %[N*m] Flexural rigidity
alpha_sdr = (4 * D / (delta_rho_c_sdr * g))^0.25;%[m] Flexural wavelength sdr
alpha_gap = (4 * D / (delta_rho_c_gap * g))^0.25;%[m] Flexural wavelength sdr
alpha_sedi = (4 * D / (delta_rho_c_sedi * g))^0.25;%[m] Flexural wavelength sdr

dV_0 = dx * g * Hd * delta_rho_d; %[N/m] volcanic line load due to denser dike

w_0_sdr = Hd * (delta_rho_d / delta_rho_c_sdr); %
w_0_gap = Hd * (delta_rho_d / delta_rho_c_gap); %
w_0_sedi = Hd * (delta_rho_d / delta_rho_c_sedi); %

dw_0_sdr = 2 * dV_0 / (alpha_sdr * g * delta_rho_c_sdr);
dw_0_gap = 2 * dV_0 / (alpha_gap * g * delta_rho_c_gap);
dw_0_sedi = 2 * dV_0 / (alpha_sedi * g * delta_rho_c_sedi);


x_sdr = 70000;
x_gap = 12000;

x_sedi = 200000;



%--------------------------------------------------------
% Plot sediments reflectors
%--------------------------------------------------
analytic_plot(Hd, delta_rho_c_gap, delta_rho_d, x...
    , alpha_gap, L, x_sedi)





xlabel('distance from the axis [km]','Fontsize',26');
ylabel('deflection [m]','Fontsize',26');
%title('analytical solutions to SDR surfaces (eq3 and eq4 of Buck2016)'...
%    ,'Fontsize',26');
set(gca,'Fontsize',26','Linewidth',2)
%axis([0 L/1000 -8000 1000])