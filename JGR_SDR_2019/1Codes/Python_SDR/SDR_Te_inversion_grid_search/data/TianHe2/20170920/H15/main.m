clc
clear
close all
% Main program for inversion to get the best fitting Te for the FLAC model

%load data
%function [particles_2]=SDR_load_data()
[particles_data]=SDR_load_data();

%n_time_start = 68;        % # of time step in FLAC model, should start from 35 for particle02 
n_time_start = 50;  
%n_time_start = 109;        % # of time step in FLAC model, should start from 35 for particle02 
%n_time_end   = 230;       % final # of steps  FLAC ran to  !!*(have to check if the top points reach stagnant stage!)
n_time_end   = 200;
dn_time = 10;              % every 5*dt kyrs   
TIME = (n_time_start:dn_time:(floor((n_time_end-n_time_start)/dn_time)*dn_time+n_time_start));


Te_guess = 20000;    %[m]
dTe = 5;          % step for changing Te [m]
max_num_trial = 100;% max number of times to try Te
tolerance = 20000;      % tolerance for misfit
alpha = 10;         % rate of updating Te_guess in gradient descent method
%alpha = 2;         % rate of updating Te_guess in gradient descent method
%misfit = 1000000*ones(length(TIME),max_num_trial);   %initial misfit
misfit = 10^12*ones(length(TIME),max_num_trial);   %initial misfit
Te_matrix = zeros(length(TIME),max_num_trial);   %initialize Te, goal is to get Te at different TIME
dJdTe = zeros(length(TIME),max_num_trial);
norm_polyfit_matrix = zeros(length(TIME),max_num_trial);
i_time = 1;    %index for time in misfit matrix
%i_time = 37;    %index for time in misfit matrix
%[misfit(i_time,num_trial)] = SDR_misfit(Te_guess,n_time_start);
%fig=figure;
%set(, 'Position', [100, 100, 1049, 895]);
for n_time = n_time_start:dn_time:n_time_end
%parfor i_time = 1:length(TIME)
%    n_time=TIME(i_time)    %print out where it goes to      
    n_time
    i_time
    num_trial = 1; %initial number of trial within each n_time
    temp_misfit = tolerance + 10; %initial temp_misfit to get into while loop
    step_Te = 11;    % Te update step size, equals to alpha * dJdTe...
    % this is to allow getting into the while loop since dJdTe is
    % initialized as zeroes
    while ((temp_misfit > tolerance) && (num_trial <= max_num_trial) ...
            && (step_Te > 2.0))%...
  %          &&(misfit(i_time,(num_trial-1))<misfit(i_time,(num_trial-2))))
        [misfit(i_time,num_trial),particle_x,particle_y,...
            particle_y_semi_analytic,norm_polyfit] = SDR_misfit(Te_guess,n_time,particles_data);
        %[misfit(i_time,num_trial+1),~,~,~] = SDR_misfit(Te_guess+dTe,n_time,particles_data);  %!change here misfit to temp
        %temp = SDR_misfit(Te_guess+dTe,n_time,particles_data);
        temp1 = SDR_misfit(Te_guess-dTe,n_time,particles_data);
        temp2 = SDR_misfit(Te_guess+dTe,n_time,particles_data);
        
        %dJdTe = (misfit(i_time,num_trial+1)-misfit(i_time,num_trial))/dTe;
        %dJdTe(i_time,num_trial) = (misfit(i_time,num_trial+1)-misfit(i_time,num_trial))/dTe;
        
        %dJdTe(i_time,num_trial) = (temp-misfit(i_time,num_trial))/dTe;
        dJdTe(i_time,num_trial) = (temp2-temp1)/dTe/2;
        
        Te_matrix(i_time,num_trial) = Te_guess;  %record the Te for the corresponding misfit
        SDR_misfit_plot(particle_x,particle_y,particle_y_semi_analytic,i_time,num_trial,n_time*20,Te_guess,misfit(i_time,num_trial));
        Te_guess = Te_guess - alpha * dJdTe(i_time,num_trial);
        %step_Te = alpha * dJdTe(i_time,num_trial); %update step_Te
        step_Te = abs(alpha * dJdTe(i_time,num_trial)); %update step_Te ***abs()
        if Te_guess < 0
            fprintf(string('error: negative Te_guess!!!!!!!!!\n'))
            Te_guess = 10000; 
        end
        temp_misfit = misfit(i_time,num_trial)
        norm_polyfit_matrix(i_time,num_trial) = norm_polyfit;
        norm_polyfit
        num_trial
        num_trial = num_trial + 1;
    end
    i_time = i_time + 1;
    save('SDR_20170914_dTe5_alpha2_tolerance20000_particles_1.mat')
end


[M,I]= min(misfit');
TE = zeros(length(TIME),1);
for i=1:length(TIME)
    TE(i) = Te_matrix(i,I(i));
end
figure
plot(TIME,TE/1000,'r+','markersize',18)
grid on;
axis([(min(TIME)-10) (max(TIME)+10) (min(TE/1000)-0.2) (max(TE/1000)+0.2)]);
xlabel('time *20 [kyrs]','Fontsize',22')
ylabel('Te [km]','Fontsize',22')
title('Te evolution with time');
set(gca,'Fontsize',22','Linewidth',2)