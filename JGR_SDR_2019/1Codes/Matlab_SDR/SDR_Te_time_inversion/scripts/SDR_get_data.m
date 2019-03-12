function [particle_x,particle_y] = SDR_get_data(time,dt,particles_data)

%plot particles, how surface move with time
%size_particle = size(particles_2);
%n = size_particle(3);
%nn = ceil(33.3333333*time);
%n=nn+1;
n=time/dt;
particle_x = particles_data(:,1,n);
particle_y = particles_data(:,2,n);

%{
figure
for i = n:n
    SDR_2 = plot(particles_2(:,1,i),particles_2(:,2,i),'b.');
    %plot(particles_2(:,1,i),particles_2(:,2,i),'b.');
    hold on
    axis([0 140 -30 1]); 
    title(['' num2str(i*30-1000) 'kyrs'],'Fontsize',26');
    pause(0.05);
    xlabel('distance from the axis [km]','Fontsize',26');
    ylabel('depth [km]','Fontsize',26');
    set(gca,'Fontsize',26','Linewidth',3)
    hold on;
end
%}

