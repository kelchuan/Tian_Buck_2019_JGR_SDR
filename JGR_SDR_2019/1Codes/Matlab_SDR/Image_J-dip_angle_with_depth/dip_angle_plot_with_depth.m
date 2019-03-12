% Plotting the analytic angle
%close all
%figure
hold on
sediment = 1;
rho_sedi = 2400;
del_rho_c_sedi = rho_c - rho_sedi;
w_0_sedi = w_0 * delta_rho_c / del_rho_c_sedi;
alpha_sedi = alpha * (delta_rho_c / del_rho_c_sedi)^0.25;


for i = 1:8
    time = sec_in_Myr * i;  % 2 Myr, x0 = 20km for Vx = 1cm/yr
    X0_flag = (x >  Vx * time);
    x_dip = x .* X0_flag;
    if (sediment == 1)
        dws_dx = 2 * w_0_sedi / alpha_sedi * (exp(-x_dip/alpha_sedi) ...
            .* cos(x_dip/alpha_sedi) - ...
        exp(-(x_dip - Vx * time)/alpha_sedi) ...
        .* cos((x_dip - Vx * time)/alpha_sedi));
        
        WS = w_0_sedi * (exp(- x_dip / alpha) .* (sin(x_dip/alpha) - cos(x_dip/alpha))...
            - exp(- (x_dip - Vx * time)/ alpha) .* (sin((x_dip - Vx * time)/alpha) - ...
            cos((x_dip - Vx * time)/alpha)));
    else
        dws_dx = 2 * w_0 / alpha * (exp(-x_dip/alpha) .* cos(x_dip/alpha) - ...
        exp(-(x_dip - Vx * time)/alpha) .* cos((x_dip - Vx * time)/alpha));
        
        WS = w_0 * (exp(- x_dip / alpha) .* (sin(x_dip/alpha) - cos(x_dip/alpha))...
            - exp(- (x_dip - Vx * time)/ alpha) .* (sin((x_dip - Vx * time)/alpha) - ...
            cos((x_dip - Vx * time)/alpha)));
    end

    WS = WS .* X0_flag;
    phi = atand(abs(dws_dx));
    

%figure
    hold on
    %plot(WS/1000, phi, 'o','color',rand(1,3)*.88)
    plot(WS/1000, phi, 'r+')

end


xlabel('km','Fontsize',26')
ylabel('phi','Fontsize',26')
if (sediment==1)
    axis([0 8 -2 20])
    title('dip angles of sediments with depth','Fontsize',26')
else
    axis([0 8 -2 20])
    title('dip angles of SDRs with depth','Fontsize',26')
end

legend('1 Myr','2 Myrs','3 Myrs','4 Myrs','5 Myrs','6 Myrs','7 Myrs'...
    ,'8 Myrs')
set(gca,'Fontsize',26','Linewidth',3)

%axis tight 