%slope
clear 
clc
close all

%initialization
delta_rho_d = 200;
delta_rho_c = 200;
g = 10;
E = 7.5 * 10^10;          %[Pa] Young's modulus
mu = 0.25;              % Poisson's ratio
x = linspace(0,3,1000);

n = 1001;
Hd_max = 25000;
Te_max = 16000;
Hd_array = linspace(0, Hd_max, n);
Te_array = linspace(0, Te_max, n);
phi = zeros(n,n);
coefficient = 2 * delta_rho_d/delta_rho_c^0.75*(E/(3*g*(1-mu^2)))^-0.25;

for i = 1:n;
    Hd = Hd_array(i);
    for j = 1:n;
        Te = Te_array(j);
        grad_dif = coefficient * Hd/(Te^0.75);
        phi(i,j) = atand(grad_dif);
    end
end

h = pcolor(phi);
set(h, 'EdgeColor', 'none');
colormap( flipud(hsv(1000)) );
colorbar
hold on;

grid on
set(gca,'layer','top')
cp1 = contour(phi,[5 5], 'k-.','linewidth', 3); hold on;
cp2 = contour(phi,[10 10], 'k-.','linewidth', 3);
cp3 = contour(phi,[15 15], 'k-.','linewidth', 3);
cp4 = contour(phi,[20 20], 'k-.','linewidth', 3);
cp5 = contour(phi,[30 30], 'k-.','linewidth', 3);
cp6 = contour(phi,[40 40], 'k-.','linewidth', 3);
cp7 = contour(phi,[50 50], 'k-.','linewidth', 3);
cp8 = contour(phi,[60 60], 'k-.','linewidth', 3);

%legend([cp1,cp2,cp3,cp4,cp5,cp6],'5','10','15','20','30','40')

%plotting model data
% EP cases
EP_Hd = [6,8,10,12,15,18,21];
EP_Te_xf = [2.19,3.11,4.2,5.24,7.12,9.19,10.41];
EP_Te_gamma = [2.2,3.25,4.37,5.58,7.51,9.58,11.76];
EP_Te_avg = 0.5 * (EP_Te_xf + EP_Te_gamma);
EP_Hd = EP_Hd * 1000 / Hd_max * n;
EP_Te_avg = EP_Te_avg * 1000 / Te_max * n;
plot(EP_Te_avg, EP_Hd, 'ro','Markersize', 18,'markerfacecolor','r'); hold on;
%EVP const therm case
EVP_Hd = [5,10,15,20];
EVP_Te_xf = [2.3,4.0,5.8,8.4];
EVP_Te_gamma = [1.6,2.6,3.9,4.4];
EVP_Te_avg = 0.5 * (EVP_Te_xf + EVP_Te_gamma);
EVP_Hd = EVP_Hd * 1000 / Hd_max * n;
EVP_Te_avg = EVP_Te_avg * 1000 / Te_max * n;
plot(EVP_Te_avg, EVP_Hd, 'yv','Markersize', 22,'markerfacecolor','y'); hold on;
%EVP mantle
EVPm_Hd = [15,15,15,15];
EVPm_Te_avg = [4.9,5.6,8.9,11.4];
EVPm_Hd = EVPm_Hd * 1000 / Hd_max * n;
EVPm_Te_avg = EVPm_Te_avg * 1000 / Te_max * n;
plot(EVPm_Te_avg([1 3]), EVPm_Hd([1 3]), 'k^','Markersize', 18,'MarkerEdgeColor','k','markerfacecolor','k');
plot(EVPm_Te_avg([2 4]), EVPm_Hd([2 4]), 'k^','Markersize', 18,'MarkerEdgeColor','k','markerfacecolor','r');
plot(EVPm_Te_avg([3 4]), EVPm_Hd([3 4]), 'go','Markersize', 5,'markerfacecolor','g');



title('Angle \phi with respect to Te and Hd','Fontsize',26')
ylabel('Hd [km]','Fontsize',26')
xlabel('Te [km]','Fontsize',26')

N_interval_Hd = 80;
N_interval_Te = 125;
set(gca,'Fontsize',26','Linewidth',3)
set(gca,'XTick', [1:N_interval_Te:n] ); %This are going to be the only values affected.
set(gca,'XTickLabel',Te_array(1:N_interval_Te:n)/1000 );
set(gca,'YTick', [1:N_interval_Hd:n] ); %This are going to be the only values affected.
set(gca,'YTickLabel',Hd_array(1:N_interval_Hd:n)/1000 );
daspect([2 1.8 1]);
