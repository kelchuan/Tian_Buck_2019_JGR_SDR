%----------------------------------------------------------------
% Iterative Implicit FD with updated distributed load q(x)
% Assuming lava fills in wherever wn flex beneath sea level (initial zero)
% q(x) = (rho_i - rho_c) * g * wn(x)
% when wn(x) >= 0, rho_i = 0; (air density)
% when wn(x) < 0, rho_i = 2800; [kg/m^3] (lava density)
% Tian 2016 Nov. @Lamont
%----------------------------------------------------------------
clc
clear all
close all
%----------------------------------------------------------------
% Initialization
%----------------------------------------------------------------
L = 200000;             %[m] total length 
N = 1001;               % number of node points
dx = L / (N - 1);       %[m] distance between node points
x = linspace(0, L, N);  %[m] initialize the x array
w = zeros(1, N);        %[m] analytic solution to deflection w(x)
dw = zeros(1, N);       %[m] analytic solution to deflection dw(x)
wn = zeros(1, N);       %[m] numerical solution to deflection

Hd = 5000;              %[m] height of the dike
Te = 5000;              %[m] effective plate thickness
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

E = 5 * 10^10;          %[Pa] Young's modulus
mu = 0.25;              % Poisson's ratio
D = (E * Te^3) / (12 * (1 - mu^2)); %[N*m] Flexural rigidity
alpha = (4 * D / (delta_rho_c * g))^0.25;%[m] Flexural wavelength

dV_0 = dx * g * Hd * delta_rho_d; %[N/m] volcanic line load due to denser dike
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
num_iterate = 6;
for i=1:1:num_iterate;
% update L base on wn   --> Assigning distributed load q(x)
if (index_lava == 0)
    L = distributed_q( L,N,wn,g,dx,rho_c,delta_rho_c, D); 
end

% update wn
wn = L\R;

if (i==num_iterate)
    %plot(x/1000, wn,'k.');hold on;
    plot(x/1000, wn,'k.');hold on;
elseif (i<=3)
    if(i ==1)
        plot(x/1000, wn,'g','Linewidth', 2);hold on;
    elseif(i==2)
        plot(x/1000, wn,'m','Linewidth', 2);hold on;
    elseif(i==3)
         plot(x/1000, wn,'r','Linewidth', 2);hold on;
    end
else
     plot(x/1000, wn,'color',rand(1,3));hold on;    
end
legend('1st', '2nd', '3rd');
axis([0 200 -35 5])
pause(0.0001);

end


xlabel('distance from the axis [km]','Fontsize',26');
ylabel('deflection [m]','Fontsize',26');
%title('iterated steady state dw(x) with half dike load of width of dx'...
%    ,'Fontsize',26')
set(gca,'Fontsize',26','Linewidth',3)
%axis([0 L/1000 -8000 1000])

%Plot the analytic solution using function analytic_plot()
%analytic_plot(Hd, delta_rho_c, delta_rho_d, x, alpha, L);



% For numerical eq3 wt and eq4 ws
figure
% For setting the size of the plot
x0=10;
y0=10;
width=1600;
height=400
set(gcf,'units','points','position',[x0,y0,width,height])
%--------------
% ploting eq3
%--------------
wt = zeros(1, N);

for i = 1:1:N
    %wt(i) = sum(wn(1:i)) * dx; %[?] why *dx doesn't work?
    wt(i) = sum(wn(1:i));
end
plot(x(3:N)/1000, wt(3:N),'k.-', 'Linewidth', 3.5);hold on;
axis([-10 210 -8000 2000])
daspect([10 1000 1])
%legend();

%--------------
% ploting eq4 ws
%--------------


for i = 0:1:20
    x_0 = i * 10000;
    ws = zeros(1,N);
    for j = (floor(x_0/dx)+1):1:N
        ws(j) = sum(wn((j-floor(x_0/dx)):j));
    end
    %plot(x(30:N)/1000, ws(30:N),'m');hold on;
    plot(x(floor(x_0/200000*N)+1:N)/1000, ws(floor(x_0/200000*N)+1:N),'k', 'Linewidth', 3);hold on;
    
end
axis([-10 L/1000+10 -8000 2000])
daspect([10 1000 1])
%legend('w_t(x) broken','w_s(x) broken');
%legend('w_t(x) ','w_s(x) continuous');
%axis([0 260 -7000 1200])
%axis([0 80 -7000 1200])
xlabel('distance from the axis [km]','Fontsize',26');
ylabel('deflection [m]','Fontsize',26');
set(gca,'Fontsize',26','Linewidth',3)



%{
%--------------------------------
% for subsidence rate of eq3 Vs
%--------------------------------
w_0 = Hd * (delta_rho_d / delta_rho_c); %
%dw_0 = 2 * dV_0 / (alpha * g * delta_rho_c);
sec_in_year = 365 * 24 * 60 * 60; 
sec_in_Myr = 10^6 * sec_in_year;
Vx = 0.01 / sec_in_year;  %half spreading rate [m/s]  assumming 1cm/yr
t = linspace(0, 9 * sec_in_Myr, 1000); % time for 9Myrs
Vs = w_0 * 2 * Vx / alpha * exp(-Vx / alpha * t) .* cos(Vx / alpha * t);
Vs = Vs * sec_in_year * 100; % transfer the unit from m/s to cm/yr

figure
plot(t / sec_in_year / 1000, Vs)
xlabel('time (kyrs)','Fontsize',26');
ylabel('Vs [cm/yr]','Fontsize',26');
set(gca,'Fontsize',26','Linewidth',3)
%axis([0 L/1000 -8000 1000])
axis tight;
title('rate of subsidence of eq3','Fontsize',26');



%--------------
% ploting numerical thin plate dip angle based on ws
%--------------
figure

for i = 0:1:10
    x_0 = i * 10000;
    ws = zeros(1,N);
    for j = (floor(x_0/dx)+1):1:N
        ws(j) = sum(wn((j-floor(x_0/dx)):j));
    end
    %plot(x(30:N)/1000, ws(30:N),'m');hold on;
    %plot(x(floor(x_0/200000*N)+1:N)/1000, ws(floor(x_0/200000*N)+1:N),'k', 'Linewidth', 3);hold on;
    %dip_func_numerical(x_0,N,ws,dx);
    dwsn_dx = diff(ws) / dx;

    phi = atand(abs(dwsn_dx));
    plot(abs(ws(1:(length(ws)-1))/1000), phi, 'o','color',rand(1,3)*.88)
    hold on;
    
    
end
xlabel('km','Fontsize',26')
ylabel('phi','Fontsize',26')

%}
%{
if sediment == 1
    axis([0 2 -2 8])
    title('dip angles of sediments with depth','Fontsize',26')
%    legend('10km', '20km', '30km', '40km', '50km', '60km');
else
    axis([0 8 -2 20])
    title('dip angles of SDRs with depth','Fontsize',26')
%    legend('10km', '20km', '30km', '40km', '50km', '60km');
end
%}

