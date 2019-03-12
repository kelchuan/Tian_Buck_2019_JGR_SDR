close all
clc
clear
%plotting the magnetic anomaly
load('data.txt');
x = data(:,2);
y = data(:,1);
z = data(:,3);

num_profiles = 20;
X = reshape(x, [length(x)/num_profiles num_profiles]);
Y = reshape(y, [length(x)/num_profiles num_profiles]);
Z = reshape(z, [length(x)/num_profiles num_profiles]);
Z = -Z;
colormap('jet')
h=pcolor(X,Y,Z)
set(h, 'EdgeColor', 'none');
h.FaceColor='interp';
h.FaceLighting='gouraud';

colorbar;
figure
colormap('jet')
kk = surf(X,Y,Z);
set(kk, 'EdgeColor', 'none');
kk.FaceColor='interp';
kk.FaceLighting='gouraud';
kkh=light('Position',[-600 220 1200])
xlabel('x')
ylabel('y')
colorbar;