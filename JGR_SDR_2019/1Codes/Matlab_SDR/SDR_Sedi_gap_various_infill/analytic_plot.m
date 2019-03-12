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
% For setting the size of the plot
x0=10;
y0=10;
width=1600;
height=400
set(gcf,'units','points','position',[x0,y0,width,height])

%plot(x_sedi_reflector/1000, -w, 'g.');  %basal line for lava sea case
light_blue = [0.6 0.7 1];   % for air
dark_blue = [0.1 0.1 1];    % for ocean
dark_green = [0.6 0.8 0.5];    % for sediment
redhot = [1 0.35 0.12];      %for lava
black = [.01 .01 .01];
colorful = black;
plot(x_sedi_reflector/1000, -w, 'color',colorful,'linewidth',5);  %basal line for lava sea case
hold on;
A = size(x);
N = A(2);
%'plot'


dx_sedi = 10000;
N_sedi = x_sedi/dx_sedi;
% equation 4
for i = 0:1:N_sedi
    x_0 = i * dx_sedi;
    w = w_0 * (exp(- x / alpha) .* (sin(x/alpha) - cos(x/alpha))...
        - exp(- (x - x_0)/ alpha) .* (sin((x - x_0)/alpha) - ...
        cos((x - x_0)/alpha)));
    if i>0
        
        plot(x(floor(x_0/L*N)+1:N)/1000, -w(floor(x_0/L*N)+1:N),'color',colorful,'linewidth',1.2);
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
axis([-10 L/1000+10 -8000 2000])
daspect([10 1000 1])

end

