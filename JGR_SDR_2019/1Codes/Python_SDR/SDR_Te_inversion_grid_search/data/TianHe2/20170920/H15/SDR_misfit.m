
%----------------------------------------------------------------
% Fitting the SDRs with semi-analytic model and reveal Te as a function of
% time 
% Tian 20170417
%----------------------------------------------------------------
function [misfit,particle_x,particle_y,particle_y_semi_analytic,norm_polyfit] = SDR_misfit(Te_guess,n_time,particles_data)
%----------------------------------------------------------------
% 0. Initialization
%----------------------------------------------------------------
%Te_guess = 7000;  % Initial guess of Te [m]
%dTe = 200;           % step for changing Te [m]
dt = 20;           % dt of FLAC models [kyrs]
%n_time = 1;        % number of model time steps, integer; should start
% from 50 for particle02
time = n_time * dt; % [kyrs], multiple of 30 kyrs
time
%time_analytic = time - 1000; % [kyrs] because we are plotting particle 2 which starts to evolve from 1Myr
%time_analytic = time - 2010; % [kyrs] because we are plotting particle 3 which starts to evolve from 2Myr
                             %Should be 2010 because that is when SDR is spreaded
%time_analytic = time - 3000;
%time_analytic = time+2000; %n_time was set to begins with 100, 100*20 = 2000
time_analytic = time; %n_time was set to begins with 100, 100*20 = 2000
%----------------------------------------------------------------
% 1. Get SDR results from semi-analytic at time T = 1Myr
%----------------------------------------------------------------
[X,Ws]=semi_analytic_results(time_analytic, Te_guess, dt);
N_points_sampling = 1000;
dN_sampling = floor(length(X)/N_points_sampling);

%plot(X(1:dN_sampling:length(X))/1000,Ws(1:dN_sampling:length(X)),'ko');
%hold on;

%----------------------------------------------------------------
% 2. Use polyfit to approximate the semi-analytic result
%----------------------------------------------------------------
npoly = 12;   % degree of the fitting polynomial
%[p,S,mu] = polyfit(X(1:dN_sampling:length(X)),Ws(1:dN_sampling:length(X)),npoly);
[p,S] = polyfit(X(1:dN_sampling:length(X)),Ws(1:dN_sampling:length(X)),npoly);


Y = polyval(p,X);  % get the fitting curve Y value
%plot(X(1:dN_sampling:length(X))/1000,Y(1:dN_sampling:length(X)),'m+');
norm_poly=S.normr;
%legend('boundary','SDR0','SDR1','SDR2','SDR3','SDR4','SDR5','SDR6',...
%    'semi-analytic','semi-analytic(sampled)','polyfit')

%ploting the error between polyfit and semi-analytic data
%figure
%plot(X(1:dN_sampling:length(X))/1000,Ws(1:dN_sampling:length(X))-Y(1:dN_sampling:length(X)),'k+');
%ylabel('absolute error')
%xlabel('x [km]')

residual=sum((Ws(1:dN_sampling:length(X))-Y(1:dN_sampling:length(X))).^2);
norm_polyfit = sqrt(residual)/length(Ws(1:dN_sampling:length(X)));
%}

% 3.0 Get data from the FLAC model
% function [particle_x,particle_y] = SDR_get_data(time)
[particle_x,particle_y]=SDR_get_data(time,dt,particles_data);
particle_x = particle_x*1000;   %trans to [m]
particle_y = particle_y*1000;

particle_y_semi_analytic = polyval(p,particle_x);  % get the fitting curve Y value

%{
%plotting  FLAC result versus semi-analytic result
figure
plot(particle_x,particle_y,'ko')   %plot result from FLAC
hold on;
plot(particle_x,particle_y_semi_analytic,'m+') % plot semi-analytic result
axis([0 120000 -30000 1000]); 
title(['' num2str(time-1000) 'kyrs'],'Fontsize',26');
xlabel('distance from the axis [m]','Fontsize',26');
ylabel('depth [m]','Fontsize',26');
set(gca,'Fontsize',22','Linewidth',2)
legend('FLAC','semi-analytic');
grid on;
%}


% calculate the misfit between semi-analytic result and FLAC
dif_model_data = particle_y_semi_analytic - particle_y;
%misfit = sqrt(sum(dif_model_data.^2))/length(particle_x);   %Euclidean norm (wiki)
misfit = sum(dif_model_data.^2)/length(particle_x)/2;   % normal way for gradient descent
%----------------------------------------------------------------
% 4. Use gradient descent to invert the best fit Te for the data
%----------------------------------------------------------------











