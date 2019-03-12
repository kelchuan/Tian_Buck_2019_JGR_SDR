%---------------------------------------------------
%ratio_plot_dike_lava_interface_XfSDR_err.m
%---------------------------------------------------
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

yyaxis left
for x0 = 1:.01:5  % x0 is the distance off axis for the tip of SDR
%x = linspace(x0,x0+1,1000);
x = x0;
y = (exp(-x).*(sin(x)-cos(x))-exp(-(x-x0)).*(sin(x-x0)-cos(x-x0)))./(exp(-x).*cos(x)-exp(-(x-x0)).*cos(x-x0));
gamma = (-(1+exp(-pi/2)));
y_diff = y - gamma;
err = y_diff/gamma * 100;
%y = (exp(-x).*(sin(x)-cos(x))+1)./(exp(-x).*cos(x)-1);
plot(x,err,'.'); hold on
%plot(x,y_diff); hold on
%grid on
ylim([-25,25])
end
xlabel("X_0/\alpha")
ylabel("Error of \gamma estimation [%]")
%title("Error of \gamma estimation along dike-lava interface")
set(gca,'FontSize',12,'FontWeight','bold') % set the font size of the ticks of the axes

%figure
yyaxis right
% Error in terms of Te
for x0 = 1:.01:5  % x0 is the distance off axis for the tip of SDR
%x = linspace(x0,x0+1,1000);
x = x0;
y = (exp(-x).*(sin(x)-cos(x))-exp(-(x-x0)).*(sin(x-x0)-cos(x-x0)))./(exp(-x).*cos(x)-exp(-(x-x0)).*cos(x-x0));
gamma = (-(1+exp(-pi/2)));
y_diff = y - gamma;
err = y_diff/gamma * 100;
Te_err = ((err / 100 + 1)^(4/3)-1) * 100;
%y = (exp(-x).*(sin(x)-cos(x))+1)./(exp(-x).*cos(x)-1);
plot(x,Te_err,'.'); hold on
%plot(x,y_diff); hold on
grid on
end
xlabel("X_0/\alpha")
ylabel("Error of Te estimation [%]")
title("Error of \gamma and Te estimation along dike-lava interface")
set(gca,'FontSize',12,'FontWeight','bold') % set the font size of the ticks of the axes

figure
%-----------------------------------------
% 2. Along SDR starting from Xf
%-----------------------------------------
yyaxis left
alpha_dist_along_single_sdr = 0.2;
for x0 = 1:.1:5  % x0 is the distance off axis for the tip of SDR
%x0 = pi/2;
x = linspace(x0,x0+alpha_dist_along_single_sdr,1000);
y = (exp(-x).*(sin(x)-cos(x))-exp(-(x-x0)).*(sin(x-x0)-cos(x-x0)))./(exp(-x).*cos(x)-exp(-(x-x0)).*cos(x-x0));
gamma = (-(1+exp(-pi/2)));
y_diff = y - gamma;
err = y_diff/gamma * 100;
%y = (exp(-x).*(sin(x)-cos(x))+1)./(exp(-x).*cos(x)-1);
plot(x,err,'.'); hold on
%plot(x,y_diff); hold on
grid on
end
ylim([-50,30])
xlabel("X_0/\alpha")
ylabel("Error of \gamma estimation [%]")
%title("Error of \gamma estimation along dike-lava interface")
set(gca,'FontSize',12,'FontWeight','bold') % set the font size of the ticks of the axes

yyaxis right
for x0 = 1:.1:5  % x0 is the distance off axis for the tip of SDR
%x0 = pi/2;
x = linspace(x0,x0+alpha_dist_along_single_sdr,1000);
y = (exp(-x).*(sin(x)-cos(x))-exp(-(x-x0)).*(sin(x-x0)-cos(x-x0)))./(exp(-x).*cos(x)-exp(-(x-x0)).*cos(x-x0));
gamma = (-(1+exp(-pi/2)));
y_diff = y - gamma;
err = y_diff/gamma * 100;
Te_err = ((err / 100 + 1).^(4/3)-1) * 100;
%y = (exp(-x).*(sin(x)-cos(x))+1)./(exp(-x).*cos(x)-1);
plot(x,Te_err,'.'); hold on
%plot(x,y_diff); hold on
grid on
end
%ylim([-25,25])
ylim([-50,30])
xlabel("X_0/\alpha")
ylabel("Error of Te estimation [%]")
title("Error of \gamma and Te estimation along each SDR from dike-lava interface to 0.2 \alpha further away from axis")
set(gca,'FontSize',12,'FontWeight','bold') % set the font size of the ticks of the axes

figure
%-----------------------------------------
% 3. Along each SDR, translating alpha to depth
%-----------------------------------------
%yyaxis left
alpha_dist_along_single_sdr = 0.2;
for x0 = 1:.1:5  % x0 is the distance off axis for the tip of SDR
%x0 = pi/2;
x = linspace(x0,x0+alpha_dist_along_single_sdr,1000);
y = (exp(-x).*(sin(x)-cos(x))-exp(-(x-x0)).*(sin(x-x0)-cos(x-x0)));
w0 = (1+exp(-pi/2)); % depth at x=x0 = Xf = pi/2
y_diff = y - w0;
err = y_diff/w0 * 100;
%y = (exp(-x).*(sin(x)-cos(x))+1)./(exp(-x).*cos(x)-1);
plot(x,err,'.'); hold on
%plot(x,y_diff); hold on
grid on
end
%ylim([-50,30])
xlabel("X_0/\alpha")
ylabel("Error of w0 (depth) estimation [%]")
title("Underestimation of depth along each SDR for 0.2 \alpha")
set(gca,'FontSize',12,'FontWeight','bold') % set the font size of the ticks of the axes
