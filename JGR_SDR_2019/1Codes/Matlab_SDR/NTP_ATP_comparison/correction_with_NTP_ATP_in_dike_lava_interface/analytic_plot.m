function analytic_plot(Hd, delta_rho_c, delta_rho_d, x, alpha, L, Length, N_SDR)
%This function plot the analytic solution to eq3 and eq4
w_0 = Hd * (delta_rho_d / delta_rho_c); %
%dw_0 = 2 * dV_0 / (alpha * g * delta_rho_c);

w = w_0 * exp(- x / alpha) .* (sin(x/alpha) - cos(x/alpha)) + w_0; %(eq.3)
%dw = dw_0 * exp(-x/alpha) .* cos(x/alpha); % (eq.1)

figure
plot(x(1:length(x)*0.4)/1000, -w(1:length(x)*0.4), 'b.');  %basal line for lava sea case
hold on;
A = size(x);
N = A(2);
% equation 4
for i = 0:1:N_SDR
    if i==N_SDR-1
        x_0=0;
    else
    x_0 = i * 4000;
    end
    w = w_0 * (exp(- x / alpha) .* (sin(x/alpha) - cos(x/alpha))...
        - exp(- (x - x_0)/ alpha) .* (sin((x - x_0)/alpha) - ...
        cos((x - x_0)/alpha)));
    if i>0
        
        plot(x(floor(x_0/Length*N)+1:N)/1000, -w(floor(x_0/Length*N)+1:N),'b.','Linewidth', 3);
        %L is the total length, not using L because L is a bug with 
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

