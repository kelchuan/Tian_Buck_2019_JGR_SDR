function analytic_plot(Hd, delta_rho_c, delta_rho_d, x, alpha, L, brok_or_cont)
%This function plot the analytic solution to eq3 and eq4
if(brok_or_cont == 1) % brok_or_cont == 1 for broken plate; 2 for cont plate
    w_0 = Hd * (delta_rho_d / delta_rho_c); %
    %dw_0 = 2 * dV_0 / (alpha * g * delta_rho_c);

    w = w_0 * exp(- x / alpha) .* (sin(x/alpha) - cos(x/alpha)) + w_0; %(eq.3)
    %dw = dw_0 * exp(-x/alpha) .* cos(x/alpha); % (eq.1)
else % for continuous plate
    w_0 = Hd * (delta_rho_d / delta_rho_c); % 
    %dw_0 = dV_0 / (alpha * g * delta_rho_c);
    w = w_0 * (1 - exp(- x / alpha) .* cos(x/alpha)); %(eq.3)
    %dw = dw_0 * exp(-x/alpha) .* (cos(x/alpha) + sin(x/alpha)); % (eq.1)
end

%figure
if(brok_or_cont == 1)
    cl = [.6 0 0] ;  %'Color',[.6 0 0]
else
    cl = [0 0.6 0];
end
plot(x/1000, -w, 'Color', cl);  %basal line for lava sea case
hold on;
A = size(x);
N = A(2);

% equation 4
if(brok_or_cont == 1) % brok_or_cont == 1 for broken plate; 2 for cont plate
    for i = 0:1:20
        x_0 = i * 5000;
        w = w_0 * (exp(- x / alpha) .* (sin(x/alpha) - cos(x/alpha))...
            - exp(- (x - x_0)/ alpha) .* (sin((x - x_0)/alpha) - ...
            cos((x - x_0)/alpha)));
        if i>0
            plot(x(floor(x_0/200000*N)+1:N)/1000,...
                -w(floor(x_0/200000*N)+1:N), 'Color', cl);
            %200000 is the total length, not using L because L is a bug with 
            % both totoal length and the Left Matrix
        %plot(x/1000, -w);
            hold on
        end
    end
else
    for i = 0:1:20
        x_0 = i * 5000;
        w = w_0 * (- exp(- x / alpha) .* cos(x/alpha)...
            + exp(- (x - x_0)/ alpha) .* cos((x - x_0)/alpha));
        if i>0
            plot(x(floor(x_0/200000*N)+1:N)/1000,...
                -w(floor(x_0/200000*N)+1:N), 'Color', cl);
            %200000 is the total length, not using L because L is a bug with 
            % both totoal length and the Left Matrix
        %plot(x/1000, -w);
            hold on
        end
    end
end
    
xlabel('Distance from the axis [km]','Fontsize',26');
ylabel('Depth [m]','Fontsize',26');
title('Analytical solutions to eq3 and eq4 of Buck 2017 EPSL'...
    ,'Fontsize',26');
set(gca,'Fontsize',26','Linewidth',3)
%axis([0 L/1000 -8000 1000])
axis tight;

end

