%-----------------------------------------
% SDR_top_corrected_err_estimation.m
%-----------------------------------------
% top corrected error estimation
clear 
clc
close all

%This script estimate based on analytic equation the error of gamma
% when measure along the dike-lava interface and along the SDR start
% from Xf
%-----------------------------------------
% 1. Along dike lava interface:  (x0 and x are non-dimensionalized with
% alpha)
%-----------------------------------------

% first, define the depth_correct, namely, the top envelop of SDRs
% ploting SDRs surfaces following SDR_top_for_loop.m
N = 1001;
t_max = 10;
t = linspace(0,t_max,N);
dx = t_max/(N-1);

N_SDR = 91;
X0m = 9;
t0 = linspace(0,X0m,N_SDR);
dt = X0m / (N_SDR - 1);

tt = zeros(N_SDR, N);
for i=1:N_SDR
    tt(i,:)=t;
end
SDR = zeros(N_SDR, N);
% ploting SDRs

for i = 1:N_SDR
    t0 = i*dt;
    index_left = floor(t0 / t_max * N)+1;
    %tt(i,index_left:N) = t(index_left:N);
    SDR(i,index_left:N) = - (exp(-tt(i,index_left:N)).*(sin(tt(i,index_left:N))-cos(tt(i,index_left:N))) - ...
      exp(-(tt(i,index_left:N)-t0)).*(sin(tt(i,index_left:N)-t0)-cos(tt(i,index_left:N)-t0)));
    %plot(tt(i,index_left:N), SDR(i,index_left:N));
    plot(tt(i,:), SDR(i,:),'.');
    hold on;
end
%figure

top_surface = zeros(1,N);
%finding the top surface of SDRs
for j = 1:N 
    top_surface(j) = max(SDR(:,j));
end

plot(t,top_surface,'bo')

% ploting envelop of SDRs and fitting
figure
plot(t,top_surface,'bo'); hold on;
index_fit_left = (floor(1.65/dx));
%index_fit_right = (floor(7/dx));
p = polyfit(t(index_fit_left:N),top_surface(index_fit_left:N),16)
%p = polyfit(t,top_surface,16)
y_fit = polyval(p,t(index_fit_left:N));
%y_fit = polyval(p,t);
%plot(t(index_fit_left:index_fit_right),y_fit,'r.')
plot(t(index_fit_left:N),y_fit,'r.')

% Now estimation the error with SDR thickness corrected based on SDRs
% envelope
% along dike-lava interface
figure
yyaxis left
N = 1000;
x0 = linspace(1,10,N);
for i = 1:N  % x0 is the distance off axis for the tip of SDR
x = x0(i);
depth = - (exp(-x).*(sin(x)-cos(x))-...
    exp(-(x-x0(i))).*(sin(x-x0(i))-cos(x-x0(i))));
tan_dip = - (exp(-x).*cos(x)-exp(-(x-x0(i))).*cos(x-x0(i)));

if (x < t(index_fit_left))
    depth_correct = 0;
else
    depth_correct = polyval(p,x);
end

SDR_thick = -(depth_correct - depth);

y = SDR_thick/ tan_dip;
gamma = (-(1+exp(-pi/2)));
y_diff = y - gamma;
err = y_diff/gamma * 100;
plot(x,err,'.'); hold on
grid on
ylim([-10,20])
end
xlabel("X_0/\alpha")
ylabel("Error of \gamma estimation [%]")
%title("Error of \gamma estimation along dike-lava interface")
set(gca,'FontSize',12,'FontWeight','bold') % set the font size of the ticks of the axes


%figure
yyaxis right
% Error in terms of Te
for i = 1:N  % x0 is the distance off axis for the tip of SDR
x = x0(i);
depth = - (exp(-x).*(sin(x)-cos(x))-...
    exp(-(x-x0(i))).*(sin(x-x0(i))-cos(x-x0(i))));
tan_dip = - (exp(-x).*cos(x)-exp(-(x-x0(i))).*cos(x-x0(i)));

if (x < t(index_fit_left))
    depth_correct = 0;
else
    depth_correct = polyval(p,x);
end

SDR_thick = -(depth_correct - depth);

y = SDR_thick/ tan_dip;
gamma = (-(1+exp(-pi/2)));
y_diff = y - gamma;
err = y_diff/gamma * 100;
Te_err = ((err / 100 + 1)^(4/3)-1) * 100;
plot(x,Te_err,'.'); hold on
grid on
ylim([-10,20])
xlim([1,6])
end
xlabel("X_0/\alpha")
ylabel("Error of Te estimation [%]")
title("Error of \gamma and Te estimation along dike-lava interface")
set(gca,'FontSize',12,'FontWeight','bold') % set the font size of the ticks of the axes


% now find out error along a single SDR from Xf
%--------------------------------------------------------------
figure
yyaxis left
N = 1000;
x0 = pi/2;
xx = linspace(x0,x0+0.2,N);
for i = 1:N  % x0 is the distance off axis for the tip of SDR
    x=xx(i);
    depth = - (exp(-x).*(sin(x)-cos(x))-...
        exp(-(x-x0)).*(sin(x-x0)-cos(x-x0)));
    tan_dip = - (exp(-x).*cos(x)-exp(-(x-x0)).*cos(x-x0));

    if (x < t(index_fit_left) )
        depth_correct = 0;
    else
        depth_correct = polyval(p,x);
    end

    SDR_thick = -(depth_correct - depth);

    y = SDR_thick/ tan_dip;
    gamma = (-(1+exp(-pi/2)));
    y_diff = y - gamma;
    err = y_diff/gamma * 100;
    plot(x,err,'.'); hold on
    grid on
    ylim([-24,0])
end
xlabel("X_0/\alpha")
ylabel("Error of \gamma estimation [%]")
%title("Error of \gamma estimation along dike-lava interface")
set(gca,'FontSize',12,'FontWeight','bold') % set the font size of the ticks of the axes


%figure
yyaxis right
% Error in terms of Te
for i = 1:N  % x0 is the distance off axis for the tip of SDR
    x=xx(i);
    depth = - (exp(-x).*(sin(x)-cos(x))-...
        exp(-(x-x0)).*(sin(x-x0)-cos(x-x0)));
    tan_dip = - (exp(-x).*cos(x)-exp(-(x-x0)).*cos(x-x0));

    if (x < t(index_fit_left) )
        depth_correct = 0;
    else
        depth_correct = polyval(p,x);
    end

    SDR_thick = -(depth_correct - depth);

    y = SDR_thick/ tan_dip;
    gamma = (-(1+exp(-pi/2)));
    y_diff = y - gamma;
    err = y_diff/gamma * 100;
    Te_err = ((err / 100 + 1)^(4/3)-1) * 100;
    plot(x,Te_err,'.'); hold on
    grid on
    ylim([-24,0])
end
xlabel("X_0/\alpha")
ylabel("Error of Te estimation [%]")
title("Error of \gamma and Te estimation along a single SDR from Xf")
set(gca,'FontSize',12,'FontWeight','bold') % set the font size of the ticks of the axes

% now still along a single SDR from Xf, but X-axis is depth 
%--------------------------------------------------------------
figure
yyaxis left
N = 1000;
x0 = pi/2;
xx = linspace(x0,x0+0.2,N);
for i = 1:N  % x0 is the distance off axis for the tip of SDR
    x=xx(i);
    depth = - (exp(-x).*(sin(x)-cos(x))-...
        exp(-(x-x0)).*(sin(x-x0)-cos(x-x0)));
    tan_dip = - (exp(-x).*cos(x)-exp(-(x-x0)).*cos(x-x0));

    if (x < t(index_fit_left) )
        depth_correct = 0;
    else
        depth_correct = polyval(p,x);
    end

    SDR_thick = -(depth_correct - depth);

    y = SDR_thick/ tan_dip;
    gamma = (-(1+exp(-pi/2)));
    y_diff = y - gamma;
    err = y_diff/gamma * 100;
    
    w0 = -(1+exp(-pi/2)); % depth at x=x0 = Xf = pi/2
    SDR_real_thick = -(depth_correct - w0);
    thick_err = (abs(SDR_thick)-abs(SDR_real_thick))...
        /abs(SDR_real_thick) * 100;
    plot(thick_err,err,'.'); hold on
    %plot(x,err,'.'); hold on
    grid on
    ylim([-24,0])
    %xlim([30,0])
end
xlabel("thickness underestimation [%]")
ylabel("Error of \gamma estimation [%]")
%title("Error of \gamma estimation along dike-lava interface")
set(gca,'FontSize',12,'FontWeight','bold') % set the font size of the ticks of the axes


%figure
yyaxis right
% Error in terms of Te
for i = 1:N  % x0 is the distance off axis for the tip of SDR
    x=xx(i);
    depth = - (exp(-x).*(sin(x)-cos(x))-...
        exp(-(x-x0)).*(sin(x-x0)-cos(x-x0)));
    tan_dip = - (exp(-x).*cos(x)-exp(-(x-x0)).*cos(x-x0));

    if (x < t(index_fit_left))
        depth_correct = 0;
    else
        depth_correct = polyval(p,x);
    end

    SDR_thick = -(depth_correct - depth);

    y = SDR_thick/ tan_dip;
    gamma = (-(1+exp(-pi/2)));
    y_diff = y - gamma;
    err = y_diff/gamma * 100;
    Te_err = ((err / 100 + 1)^(4/3)-1) * 100;
    
    w0 = -(1+exp(-pi/2)); % depth at x=x0 = Xf = pi/2
    SDR_real_thick = -(depth_correct - w0);
    thick_err = (abs(SDR_thick)-abs(SDR_real_thick))...
        /abs(SDR_real_thick) * 100;
    plot(thick_err,Te_err,'.'); hold on
    %plot(x,err,'.'); hold on
    grid on
    ylim([-24,0])
    %xlim([30,0])
end
xlabel("SDRs thickness underestimation [%]")
ylabel("Error of Te estimation [%]")
title("Error of \gamma and Te estimation along a single SDR from Xf")
set(gca,'FontSize',12,'FontWeight','bold') % set the font size of the ticks of the axes
set(gca, 'XDir','reverse')
xlim([-30,0])

