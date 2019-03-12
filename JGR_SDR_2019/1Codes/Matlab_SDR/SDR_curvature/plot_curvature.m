clear 
clc
close all

% analytic plot of curvature of SDRs
w0 = 5000;
x0 = 2000;
%x0 = 1000;
alpha = 32500;
x =  linspace(0,200000, 100000);

flag = x > x0;
[c1, ia1, ic1]=unique(flag);

d2wdx2 = -2 * w0 / alpha^2 * (exp(-x / alpha) .* (sin(x/alpha) + cos(x/alpha)) - ...
    exp(-(x-x0) / alpha) .* (sin((x-x0)/alpha) + cos((x-x0)/alpha)));  

d3wdx3 = -4*w0/alpha^3 * (exp(-x / alpha) .* (-sin(x/alpha)) + ...
    exp(-(x-x0) / alpha).*(sin((x-x0)/alpha)));  

plot(x(ia1(2):length(x))/1000, d2wdx2(ia1(2):length(d2wdx2)), '--', x/1000, zeros(size(x)))
%axis tight
xlabel('off-axis distance [km]','Fontsize',26')
ylabel('curvature [1/m]','Fontsize',26')
set(gca,'Fontsize',26','Linewidth',3)
grid on

xm = alpha * atan((1 / tan(x0/alpha) - exp(-x0/alpha)/sin(x0/alpha))^(-1));
% xm is the distance off axis of the location of the maximum curvature 

exp(-xm/alpha)*sin(xm/alpha)

exp(-(xm-x0)/alpha)*sin((xm-x0)/alpha)


figure
for x0 = 100:1000:1.3*alpha
    x0
    xm = alpha * atan((1 / tan(x0/alpha) - exp(-x0/alpha)/sin(x0/alpha))^(-1));
    plot(x0,xm,'r+', x0, (xm - x0), 'b.')
    xlabel('x0')
    ylabel('xm and xm - x0')
    hold on
    daspect([1 1 1])
end
%{
figure
for x0 = 100:500:1.3*alpha
    x0
    xm = alpha * atan((1 / tan(x0/alpha) - exp(-x0/alpha)/sin(x0/alpha))^(-1));
    plot(x0,xm,'r+', x0, (xm - x0), 'b.')
    hold on
    daspect([1 1 1])
end
%}