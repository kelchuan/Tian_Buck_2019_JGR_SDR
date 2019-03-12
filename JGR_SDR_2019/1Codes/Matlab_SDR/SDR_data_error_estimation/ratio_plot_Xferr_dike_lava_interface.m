clear 
clc
close all

%This script estimate based on analytic equation the error of Xf
% when measure along the dike-lava interface
%-----------------------------------------
% 1. Along dike lava interface:  (x0 and x are non-dimensionalized with
% alpha)
%-----------------------------------------

yyaxis left
for x0 = (pi/2-0.2):.0001: (pi/2+0.2) % x0 is the distance off axis for the tip of SDR
Xf = pi/2;
%x = linspace(x0,x0+1,1000);
x_diff = x0 - Xf;
err = x_diff/Xf * 100;
%y = (exp(-x).*(sin(x)-cos(x))+1)./(exp(-x).*cos(x)-1);
plot(x0,err,'.'); hold on
%plot(x,y_diff); hold on
%grid on
ylim([-20,20])
end
xlabel("X_0/\alpha")
ylabel("Error of Xf estimation [%]")
%title("Error of \gamma estimation along dike-lava interface")
set(gca,'FontSize',12,'FontWeight','bold') % set the font size of the ticks of the axes

%figure
yyaxis right
% Error in terms of Te
for x0 = (pi/2-0.2):.0001: (pi/2+0.2)  % x0 is the distance off axis for the tip of SDR
Xf = pi/2;
%x = linspace(x0,x0+1,1000);
x_diff = x0 - Xf;
err = x_diff/Xf * 100;
Te_err = ((err / 100 + 1)^(4/3)-1) * 100;
%y = (exp(-x).*(sin(x)-cos(x))+1)./(exp(-x).*cos(x)-1);
plot(x0,Te_err,'.'); hold on
%plot(x,y_diff); hold on
grid on
end
xlabel("X_0/\alpha")
ylabel("Error of Te estimation [%]")
title("Error of Xf and Te estimation along dike-lava interface with 0.2 \alpha closer or further away from Xf")
set(gca,'FontSize',12,'FontWeight','bold') % set the font size of the ticks of the axes

