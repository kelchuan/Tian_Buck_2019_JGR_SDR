%----------------------------------------------------------------
% Iterative Implicit FD with updated distributed load q(x)
% Assuming lava fills in wherever wn flex beneath sea level (initial zero)
% q(x) = (rho_i - rho_c) * g * wn(x)
% when wn(x) >= 0, rho_i = 0; (air density)
% when wn(x) < 0, rho_i = 2800; [kg/m^3] (lava density)
% Tian 2016 Nov. @Lamont
%----------------------------------------------------------------

function [X,Ws] = semi_analytic_results(time, Te, dt)

%initialization
%function [L,N,dx,x,w,dw,wn,Hd,g,delta_rho_d,rho_c,delta_rho_c,D,dV_0,Te] = semi_analytic_init(Te);
[L,N,dx,x,~,~,wn,~,g,~,rho_c,delta_rho_c,D,dV_0,~] = semi_analytic_init(Te);
%function [L] = FDM_setup(N,D,dx,dV_0);
[Lm,R,index_lava] = FDM_setup(N,D,dx,dV_0);

%figure
num_iterate = 6;
for i=1:1:num_iterate
% update L base on wn   --> Assigning distributed load q(x)
if (index_lava == 0)
    Lm = distributed_q( Lm,N,wn,g,dx,rho_c,delta_rho_c, D); 
end

% update wn
wn = Lm\R;

%Plotting of wn
%{
if (i==num_iterate)
    plot(x/1000, wn,'k.');hold on;
else
     plot(x/1000, wn,'color',rand(1,3));hold on;    
end
%axis([0 200 -50 5])
pause(0.0001);

xlabel('distance from the axis [km]','Fontsize',26');
ylabel('deflection [m]','Fontsize',26');
title('iterated steady state dw(x) with half dike load of width of dx'...
    ,'Fontsize',26')
set(gca,'Fontsize',26','Linewidth',3)
axis([0 L/1000 -50 5])
%}
end




%Plot the analytic solution using function analytic_plot()
%analytic_plot(Hd, delta_rho_c, delta_rho_d, x, alpha, L);


vx = 10000;    % [meter/Myr] half spreading rate  10km/Myr = 1cm/yr
%{
% For numerical eq3 wt and eq4 ws
figure  %figure 1
%--------------
% ploting eq3
%--------------
wt = zeros(1, N);

for i = 1:1:N
    %wt(i) = sum(wn(1:i)) * dx; %[?] why *dx doesn't work?
    wt(i) = sum(wn(1:i));
end
plot(x(3:N)/1000, wt(3:N),'r.-', 'Linewidth', 3.5);hold on;
%legend();


%--------------
% ploting eq4 ws
%--------------
for i = 0:2:12
    
    x_0 = i * vx;
    ws = zeros(1,N);
    for j = (floor(x_0/dx)+1):1:N
        ws(j) = sum(wn((j-floor(x_0/dx)):j));
    end
    %plot(x(floor(x_0/200000*N)+1:N)/1000, ws(floor(x_0/200000*N)+1:N),'k', 'Linewidth', 2);hold on;
    plot(x(floor(x_0/200000*N)+1:N)/1000, ws(floor(x_0/200000*N)+1:N),...
        'color',rand(1,3),'Linewidth', 2);hold on;  
end
%}

time_Myr = time/1000;   % time in Myr
x0 = time_Myr * vx;
for j = (floor(x0/dx)+1):1:N
        Ws(j) = sum(wn((j-floor(x0/dx)):j));
end
X = x(floor(x0/L*N)+1:N);
Ws = Ws(floor(x0/L*N)+1:N);
%{
plot(X/1000,Ws,'bo', 'Linewidth', 1, 'markersize', 12);hold on;
%axis([0 L/1000 -26000 1200])
%axis([0 L/1000 -26000 1200])
xlabel('distance from the axis [km]','Fontsize',26');
ylabel('deflection [m]','Fontsize',26');
set(gca,'Fontsize',18','Linewidth',2)
%}

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
%}


%--------------
% ploting numerical thin plate dip angle based on ws
%--------------
%{
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

if sediment == 1
    axis([0 2 -2 8])
    title('dip angles of sediments with depth','Fontsize',26')
    legend('10km', '20km', '30km', '40km', '50km', '60km');
else
    axis([0 8 -2 20])
    title('dip angles of SDRs with depth','Fontsize',26')
    legend('10km', '20km', '30km', '40km', '50km', '60km');
end
%}

%legend('w_t(x) broken','w_s(x) broken');
%legend('w_t(x) ','w_s(x) continuous');
%axis([0 160 -7000 1200])
%axis([0 80 -7000 1200])
%xlabel('distance from the axis [km]','Fontsize',26');
%ylabel('deflection [m]','Fontsize',26');
%set(gca,'Fontsize',26','Linewidth',3)
