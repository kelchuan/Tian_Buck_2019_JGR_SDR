function analytic_plot(Hd, delta_rho_c, delta_rho_d, x, alpha, L, x_sedi)
%This function plot the analytic solution to eq3 and eq4
w_0 = Hd * (delta_rho_d / delta_rho_c); %
%dw_0 = 2 * dV_0 / (alpha * g * delta_rho_c);

sedi_flag = (x<=x_sedi);
x_sedi_reflector = x(find(sedi_flag,1):find(sedi_flag,1,'last'));

%w = w_0 * exp(- x / alpha) .* (sin(x/alpha) - cos(x/alpha)) + w_0; %(eq.3)
w = w_0 * exp(- x_sedi_reflector / alpha) .* (sin(x_sedi_reflector/alpha)...
    - cos(x_sedi_reflector/alpha)) + w_0; %(eq.3)
%dw = dw_0 * exp(-x/alpha) .* cos(x/alpha); % (eq.1)

%flag_sedi = (x<x_sedi);

%figure
plot(x_sedi_reflector/1000, -w, 'g.');  %basal line for lava sea case
hold on;
A = size(x);
N = A(2);
%'plot'


dx_sedi = 5000;
N_sedi = x_sedi/dx_sedi;
% equation 4
for i = 0:1:N_sedi
    x_0 = i * dx_sedi;
    w = w_0 * (exp(- x / alpha) .* (sin(x/alpha) - cos(x/alpha))...
        - exp(- (x - x_0)/ alpha) .* (sin((x - x_0)/alpha) - ...
        cos((x - x_0)/alpha)));
    if i>0
        
        plot(x(floor(x_0/L*N)+1:N)/1000, -w(floor(x_0/L*N)+1:N),'g.');
        %200000 is the total length, not using L because L is a bug with 
        % both totoal length and the Left Matrix
    %plot(x/1000, -w);
    hold on
end
xlabel('distance from the axis [km]','Fontsize',26');
ylabel('deflection [m]','Fontsize',26');
%title('analytical solutions to SDR surfaces (eq3 and eq4 of Buck2016)'...
%    ,'Fontsize',26');
set(gca,'Fontsize',26','Linewidth',3)
%axis([0 L/1000 -8000 1000])
%axis tight;


end

