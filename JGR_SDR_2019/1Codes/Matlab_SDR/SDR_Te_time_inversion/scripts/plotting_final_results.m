%[M,I]= min(misfit');
%n_time_start = 35;        % # of time step in FLAC model, should start from 35 for particle02 
%n_time_end   = 190;       % final # of steps  FLAC ran to  !!*(have to check if the top points reach stagnant stage!)
%dn_time = 2;              % every 5*dt kyrs   
%TIME = (n_time_start:dn_time:(floor((n_time_end-n_time_start)/dn_time)*dn_time+n_time_start));

misfit_final = misfit;
Te_final = Te_matrix;


%TE = zeros((length(TIME)-1),1);
TE = zeros((length(TIME)),1);
%for i=1:(length(TIME)-1)
for i=1:(length(TIME))
    [M,I,kk]=unique(misfit_final(i,:));
    if (i>=80)
        TE(i) = Te_final(i,I(2));
      MISFIT(i)= misfit_final(i,I(2));
    else
        TE(i) = Te_final(i,I(1));
        MISFIT(i)= misfit_final(i,I(1));
   end
end
%{
for i=1:(length(TIME)-1)
    [M,I,kk]=unique(misfit_final(i,:));
    if i > 41
        TE(i) = Te_final(i,I(2));
        MISFIT(i)= misfit_final(i,I(2));
    else
        TE(i) = Te_final(i,I(1));
        MISFIT(i)= misfit_final(i,I(1));
    end
end
%}
figure
%plot(TIME(1:77),TE/1000,'r+','markersize',18,'linewidth',2.5)
plot(TIME*30,TE/1000,'ko','markersize',12,'linewidth',1.6)
hold on
%plot(TIME(1:77),TE/1000,'b-.','linewidth',2)
plot(TIME*30,TE/1000,'b-.','linewidth',2)
grid on;
%http://www.mathworks.com/help/matlab/ref/plotyy.html#bt405xa-1
%[hAx,hLine1,hLine2] = plotyy(TIME(1:77),TE/1000,TIME(1:77),sqrt(MISFIT))
[hAx,hLine1,hLine2] = plotyy(TIME*30,TE/1000,TIME*30,sqrt(MISFIT))
hLine2.LineStyle = '-.';
ylabel(hAx(2),'Misfit [m]','Fontsize',22') % right y-axis
axis([(min(TIME*30)-120) (max(TIME*30)+120) (min(TE/1000)-1.0) (max(TE/1000)+0.6)]);
%xlabel('time *30 [kyrs]','Fontsize',22')
xlabel('[kyrs]','Fontsize',22')
ylabel('Te [km]','Fontsize',22')
title('Te evolution with time');
set(gca,'Fontsize',22','Linewidth',2)
set(hAx(2),'Fontsize',22','Linewidth',4)
set(hLine2,'Linewidth',2)
legend('Te','','','Misfit')