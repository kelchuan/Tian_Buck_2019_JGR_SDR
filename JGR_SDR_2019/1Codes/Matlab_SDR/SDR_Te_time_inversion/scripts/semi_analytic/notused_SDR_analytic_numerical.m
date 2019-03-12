% This matlab script plots the analytic solution to Buck's 2016 SDR paper 
% with updated gap and sediment infill after SDRs
% Tian 2017 April

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
rho_i_gap = 1000;           %[kg/m^3] density of the gap(water) infill
rho_i_sedi = 2400;           %[kg/m^3] density of the sediment infill
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

%w = w_0 * exp(- x / alpha) .* (sin(x/alpha) - cos(x/alpha)) + w_0; %(eq.3)
%dw = dw_0 * exp(-x/alpha) .* cos(x/alpha); % (eq.1)
x_sdr = 7000;
x_gap = 12000;
%x_gap = 0;
x_sedi = 10000;

w_sdr_end = w_0_sdr * (exp(- x_sdr / alpha_sdr) .* (sin(x_sdr/alpha_sdr)...
    - cos(x_sdr/alpha_sdr)) + 1);
w_gap_end = w_0_gap * (exp(- (x_sdr + x_gap) / alpha_gap) .* ...
    (sin((x_sdr + x_gap)/alpha_gap)- cos((x_sdr + x_gap)/alpha_gap)) -...
    exp(- x_sdr / alpha_gap) .* (sin(x_sdr/alpha_gap)...
    - cos(x_sdr/alpha_gap)));
%{
w_os = w_0_sedi * (exp(- x / alpha_sedi) .* ...
    (sin(x/alpha_sedi)- cos(x/alpha_sedi)) -...
    exp(- (x_sdr + x_gap) / alpha_sedi) .* ...
    (sin((x_sdr + x_gap)/alpha_sedi)- cos((x_sdr + x_gap)/alpha_sedi)));
%}

w_sdr = w_0_sdr * (exp(- x / alpha_sdr) .* (sin(x/alpha_sdr)...
    - cos(x/alpha_sdr)) + 1);
w_gap = w_0_gap * (exp(- (x) / alpha_gap) .* ...
    (sin((x)/alpha_gap)- cos((x)/alpha_gap)) -...
    exp(- x_sdr / alpha_gap) .* (sin(x_sdr/alpha_gap)...
    - cos(x_sdr/alpha_gap)));
w_os = w_0_sedi * (exp(- x / alpha_sedi) .* ...
    (sin(x/alpha_sedi)- cos(x/alpha_sedi)) -...
    exp(- (x_sdr + x_gap) / alpha_sedi) .* ...
    (sin((x_sdr + x_gap)/alpha_sedi)- cos((x_sdr + x_gap)/alpha_sedi)));

%w = w_sdr + w_gap + w_os;
% Plotting the boundary between dike and SDRs/sediments
figure
flag_sdr = (x <= x_sdr);
plot(x(find(flag_sdr,1):find(flag_sdr,1,'last'))/1000,...
    -w_sdr(find(flag_sdr, 1):find(flag_sdr,1,'last')), 'r.');
hold on;
flag_gap = ((x<=(x_sdr+x_gap)) .* (x>x_sdr));
plot(x(find(flag_gap, 1):find(flag_gap, 1, 'last'))/1000,...
   -(w_sdr_end+w_gap(find(flag_gap,1):find(flag_gap, 1, 'last' ))), 'b.');
hold on;
flag_sedi = (x > (x_sdr + x_gap));
plot(x(find(flag_sedi,1):find(flag_sedi, 1, 'last' ))/1000, ...
    -(w_sdr_end+w_gap_end...
    +w_os(find(flag_sedi,1):find(flag_sedi,1,'last'))), 'g.');
hold on;
%axis([0 L/1000/2 -6000 1000])

% Plot sediments reflectors

analytic_plot(Hd, delta_rho_c_sedi, delta_rho_d, x...
    , alpha_sedi, L, x_sedi)

%{
% equation 4 with first SDR, then GAP, last SEDIMENTS
sdr_lava = w_0_sdr * (exp(- (x-x_gap-x_sedi) / alpha_sdr) .* ...
    (sin((x-x_gap-x_sedi)/alpha_sdr) - cos((x-x_gap-x_sedi)/alpha_sdr)) - ...
    exp(- (x-x_gap-x_sedi-x_sdr) / alpha_sdr) .* ...
    (sin((x-x_gap-x_sedi-x_sdr)/alpha_sdr) ...
    - cos((x-x_gap-x_sedi-x_sdr)/alpha_sdr)));
sdr_gap = w_0_gap * (exp(- (x-x_sedi) / alpha_gap) .* ...
    (sin((x-x_sedi)/alpha_gap) - cos((x-x_sedi)/alpha_gap)) - ...
    exp(- (x-x_gap-x_sedi) / alpha_gap) .* ...
    (sin((x-x_gap-x_sedi)/alpha_gap) ...
    - cos((x-x_gap-x_sedi)/alpha_gap)));
sdr_sedi = w_0_sedi * (exp(- (x) / alpha_sedi) .* ...
    (sin((x)/alpha_sedi) - cos((x)/alpha_sedi)) - ...
    exp(- (x-x_sedi) / alpha_sedi) .* ...
    (sin((x-x_sedi)/alpha_sedi) ...
    - cos((x-x_sedi)/alpha_sedi)));

sdrs = sdr_lava + sdr_gap + sdr_sedi;
plot(x(floor((x_sdr+x_sedi+x_gap)/200000*N)+1:N)/1000, ...
    -sdrs(floor((x_sdr+x_sedi+x_gap)/200000*N)+1:N), 'r');
hold on;

sdrs = sdr_lava + sdr_gap;
plot(x(floor((x_sdr+x_gap)/200000*N)+1:N)/1000, ...
    -sdrs(floor((x_sdr+x_gap)/200000*N)+1:N), 'g');
%}



% plotting individual reflector

%SDR+gap+sedi
analytic_plot_sedi_gap(w_0_sdr, w_0_gap, w_0_sedi, ...
    x_sdr, x_gap, x_sedi, alpha_sdr, alpha_gap, alpha_sedi, x, N, 'g');
% Comparison with only continuous lava infill, namely SDR only
x_0 = x_sdr+x_sedi+x_gap;
w = w_0_sdr * (exp(- x / alpha_sdr) .* (sin(x/alpha_sdr) - cos(x/alpha_sdr))...
        - exp(- (x - x_0)/ alpha_sdr) .* (sin((x - x_0)/alpha_sdr) - ...
        cos((x - x_0)/alpha_sdr)));
plot(x(floor(x_0/200000*N)+1:N)/1000, -w(floor(x_0/200000*N)+1:N), 'r+');
hold on;
%SDR+gap
x_sedi = 0;
analytic_plot_sedi_gap(w_0_sdr, w_0_gap, w_0_sedi, ...
    x_sdr, x_gap, x_sedi, alpha_sdr, alpha_gap, alpha_sedi, x, N, 'b');
% Comparison with only continuous lava infill, namely SDR only
x_0 = x_sdr+x_sedi+x_gap;
w = w_0_sdr * (exp(- x / alpha_sdr) .* (sin(x/alpha_sdr) - cos(x/alpha_sdr))...
        - exp(- (x - x_0)/ alpha_sdr) .* (sin((x - x_0)/alpha_sdr) - ...
        cos((x - x_0)/alpha_sdr)));
plot(x(floor(x_0/200000*N)+1:N)/1000, -w(floor(x_0/200000*N)+1:N), 'r+');
hold on;
%SDR only
x_sedi = 0;
x_gap = 0;
analytic_plot_sedi_gap(w_0_sdr, w_0_gap, w_0_sedi, ...
    x_sdr, x_gap, x_sedi, alpha_sdr, alpha_gap, alpha_sedi, x, N, 'r');
% Comparison with only continuous lava infill, namely SDR only
x_0 = x_sdr+x_sedi+x_gap;
w = w_0_sdr * (exp(- x / alpha_sdr) .* (sin(x/alpha_sdr) - cos(x/alpha_sdr))...
        - exp(- (x - x_0)/ alpha_sdr) .* (sin((x - x_0)/alpha_sdr) - ...
        cos((x - x_0)/alpha_sdr)));
%plot(x(floor(x_0/200000*N)+1:N)/1000, -w(floor(x_0/200000*N)+1:N), 'r+');
hold on;



%{
for i = 0:1:10
    x_0 = i * 5000;
    w = w_0_sdr * (exp(- x / alpha_sdr) .* (sin(x/alpha_sdr) - cos(x/alpha_sdr))...
        - exp(- (x - x_0)/ alpha_sdr) .* (sin((x - x_0)/alpha_sdr) - ...
        cos((x - x_0)/alpha_sdr)));
    plot(x(floor(x_0/200000*N)+1:N)/1000, -w(floor(x_0/200000*N)+1:N), 'r+');
    hold on
end
hold on;

analytic_plot(Hd, delta_rho_c_sdr, delta_rho_d, x, alpha_sdr, L);
%}




xlabel('distance from the axis [km]','Fontsize',26');
ylabel('deflection [m]','Fontsize',26');
title('analytical solutions to SDR surfaces (eq3 and eq4 of Buck2016)'...
    ,'Fontsize',26');
set(gca,'Fontsize',26','Linewidth',2)
%axis([0 L/1000 -8000 1000])


%{
% Finite difference solution to D * d4w/dx4 + q(x) = 0
% Implicit solution L * W  = R
% L are the corresponding coefficients of w(i-2 ~ i+2)
% W are all the w
% R are the right hand side which are the extra loads


% Lava sea or not
index_lava = 1;
% index_lava = 1 means lava sea where:
                        % q(x) = (rho_i - rho_c) * g * w(x) for any w(x)
% index_lava = 2 means:
                        % q(x) = (rho_i - rho_c) * g at w(x) < 0
                        % q(x) = - rho_c * g at w(x) >= 0
%if (index_lava = 1)

%broken plate
% setup the L matrix (coefficients of w)
% D/dx^4 * (w(i+2) - 4*w(i+1) + 6*w(i) - 4*w(i-1) + w(i-2)) + ...
% delta_rho_c * g * w(i) = 0
% Setup the L matrix
num_i = N;  
num_j = N - 4;
L = zeros(N, N);
%{
coef_i_minus_2 = D / dx^4;               % coefficient for w(i-2) (j == i)
coef_i_minus_1 = -4 * D / dx^4;          % coefficient for w(i-1) (i = j+1)
coef_i = 6 * D / dx^4 + delta_rho_c * g; % coeeficient for w(i)   (i = j+2)
coef_i_plus_1 = -4 * D / dx^4;           % coefficient for w(i+1) (i = j+3)
coef_i_plus_2 = D / dx^4;                % coefficient for w(i+2) (i = j+4)
%}
coef_i_minus_2 = 1;                         % coefficient for w(i-2) (j == i)
coef_i_minus_1 = -4;                        % coefficient for w(i-1) (i = j+1)
coef_i = 6 + delta_rho_c * g / (D / dx^4);  % coeeficient for w(i)   (i = j+2)
coef_i_plus_1 = -4;                         % coefficient for w(i+1) (i = j+3)
coef_i_plus_2 = 1;                          % coefficient for w(i+2) (i = j+4)

for j = 1:1:num_j
    for i = 1:1:num_i
        if(i == j)
           L(j,i) = coef_i_minus_2;
        elseif(i==j+1)
            L(j,i) = coef_i_minus_1;
        elseif(i==j+2)
            L(j,i) = coef_i;
        elseif(i==j+3)
            L(j,i) = coef_i_plus_1;
        elseif(i==j+4)
            L(j,i) = coef_i_plus_2;
        end
    end
end
%spy(L) % use spy function to show the sparse pattern of matrix L


% Setup the R matrix (N-5+1+4, 1)  % THe added 4 is for 4 BCs

    R = zeros(N-5+1+4, 1);


% Boundary conditions
index_broken = 1; % if index_broken = 1, then it is calculating broken plate
                  % if index_broken != 1, it is for continuous plate

% BC1 approximate line load due to denser dike at the center 
% R(1) = -dV_0 / 2 / (D / dx^4);
% Rather than applying on the RHS of the first equation, we use BCs of 
% V = dV_0 at x = 0

if (index_broken == 1)
    % The remaining four BCs can be added to the L and R matrices
    % The row number of the four BCs is random and interchangeable
    % BC1&2 w(inf) = 0; w'(inf) = 0
    % wn(N) = 0;
    L(num_j+1,num_i) = 1;
    % w'(inf) = 0;
    % wn(N-1) = wn(N);
    L(num_j+2,num_i) = 1;
    L(num_j+2,num_i-1) = -1;
    % BC3 M = 0 at x = 0  --> d2w/dx2|(x=0) = 0
    % d2w/dx2 = (wn(i+1) - 2wn(i) + wn(i-1)) / (dx)^2;
    % M = -D * d2w/dx2 = -D * (wn(i+1) - 2wn(i) + wn(i-1)) / (dx)^2;
    L(num_j+3, 1) = 1;
    L(num_j+3, 2) = -2;
    L(num_j+3, 3) = 1;
    % BC4 V = dM/dx = -D * d3w/dx3 = V_0/2 = dV_0 at x = 0
    % d3w/dx3|(i-0.5) = (d2w/dx2|i - d2w/dx2|i-1) / dx
    %               = (wn(i+1) -3*wn(i)+3*wn(i-1)-wn(i-2)) / (dx)^3
    % -D * (-wn(1) + 3 * wn(2) - 3 * wn(3) + wn(4)) / dx^3;
    L(num_j+4, 1) = -1;
    L(num_j+4, 2) = 3;
    L(num_j+4, 3) = -3;         
    L(num_j+4, 4) = 1;
    R(num_j+4,1) = dV_0 / (-D/dx^3); % value at the RHS
    
else % Continuous plate (index_broken != 1)
    
   % BC2&3 w(inf) = 0; w'(inf) = 0
    %wn(N) = 0;
    L(num_j+1,num_i) = 1;
    %wn(N-1) = wn(N);
    L(num_j+2,num_i) = 1;
    L(num_j+2,num_i-1) = -1;
    % BC4 wn' = 0 at x = 0
    % wn(1) = wn(2);
    L(num_j+3, 1) = 1;
    L(num_j+3, 2) = -1;

    % BC5 V = 0 at x = 0
    %wn(4) = wn(1) - 3 * wn(2) + 3 * wn(3);
    L(num_j+4, 1) = -1;
    L(num_j+4, 2) = 3;
    L(num_j+4, 3) = -3;
    L(num_j+4, 4) = 1;
end

%elseif (index_lava = 2) % q(x) = (rho_i - rho_c) * g at w(x) < 0
                        % q(x) = - rho_c * g at w(x) >= 0



% Use backslash solver "\" to get W
wn = L\R;
%figure

if (index_lava == 2)
    L = distributed_q( L,N,wn,g,dx,rho_c,delta_rho_c, D); 
end

wn = L\R;


%Analytic solution and numerical solution
figure

subplot(2,1,1)
plot(x(1:20:N)/1000, -dw(1:20:N),'bo','Markersize', 12);
hold on;
plot(x/1000, wn,'k');
hold on
plot(x/1000, 0);
xlabel('distance from the axis [km]','Fontsize',26');
ylabel('deflection [m]','Fontsize',26');
title('(equation 1) dw(x) numerical and analytical solutions comparison'...
    ,'Fontsize',26');
legend('analytic', 'implicit FDM');
set(gca,'Fontsize',26','Linewidth',3)

subplot(2,1,2)
plot(x/1000, abs(wn+dw'), 'r');
set(gca,'Fontsize',26','Linewidth',3)
ylabel('absolute error','Fontsize',26');
xlabel('distance from the axis [km]','Fontsize',26');
title('absolute error of wn-dw');

%}

   










