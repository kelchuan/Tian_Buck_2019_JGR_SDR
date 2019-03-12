function [L,N,dx,x,w,dw,wn,Hd,g,delta_rho_d,rho_c,delta_rho_c,D,dV_0,Te] = semi_analytic_init(Te)

%----------------------------------------------------------------
% Initialization
%----------------------------------------------------------------
L = 200000;             %[m] total length 
N = 2001;               % number of node points
dx = L / (N - 1);       %[m] distance between node points
x = linspace(0, L, N);  %[m] initialize the x array
w = zeros(1, N);        %[m] analytic solution to deflection w(x)
dw = zeros(1, N);       %[m] analytic solution to deflection dw(x)
wn = zeros(1, N);       %[m] numerical solution to deflection

%Hd = input('enter the dike height [m]');  %[m] height of the dike
Hd = 20000;
%Te = 10000;              %[m] effective plate thickness
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
mu = 0.25;              % Poisson's ratio
D = (E * Te^3) / (12 * (1 - mu^2)); %[N*m] Flexural rigidity
alpha = (4 * D / (delta_rho_c * g))^0.25;%[m] Flexural wavelength
scale_factor = 1; %for 1km wide element, dx here is 100m, so need scale factor of 10
dV_0 = dx * g * Hd * delta_rho_d * scale_factor; %[N/m] volcanic line load due to denser dike