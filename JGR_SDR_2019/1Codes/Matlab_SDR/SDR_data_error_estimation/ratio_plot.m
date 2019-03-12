clear 
clc
close all
%for x0 = 0.01:.05:6
for x0 = 0.01:.05:6
%x0 = 0.3
x = linspace(x0,x0+1,1000);
%x = linspace(x0,x0+1.4,1000);
y = (exp(-x).*(sin(x)-cos(x))-exp(-(x-x0)).*(sin(x-x0)-cos(x-x0)))./(exp(-x).*cos(x)-exp(-(x-x0)).*cos(x-x0));
%y = (exp(-x).*(sin(x)-cos(x))+1)./(exp(-x).*cos(x)-1);
plot(x,y); hold on
grid on
end
xlabel("x_0/alpha")
ylabel("gamma(x_0)/(alpha/2)")
title("gamma at x_0 with distance off axis normalized by alpha")

