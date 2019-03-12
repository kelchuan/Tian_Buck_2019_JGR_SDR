%Taylor expansion of Wt
%https://www.mathworks.com/help/symbolic/taylor.html
clear all
clc 
close all
figure
%x = linspace(0,200,200);
syms x
alpha = 50;
w0 = 5;
f = w0.*(exp(-x/alpha).*(sin(x/alpha)-cos(x/alpha))+1);

%t2 = taylor(f, x, 'Order', 2)
t4 = taylor(f, x, 'Order', 4)
t6 = taylor(f, x, 'Order', 6)
t8 = taylor(f, x, 'Order', 8)
t16 = taylor(f, x, 'Order', 16)
fplot([t4 t6 t8 t16 f]); hold on;
plot(ones(100)* pi/2 * alpha, linspace(-2+w0,2+w0,100)); hold on;
xlim([-5 200])
ylim([-10 20])
grid on

legend('O2','O5','O6','O16','f')


x = linspace(0,200,200)
%plot(x,eval(t16),'r-.')
[pks loc]=findpeaks(eval(t16))
Xf = pi/2 * alpha
findpeaks(eval(t16))
%fplot([t10 f]); hold on
%[xmin,fval] = fminbnd(@t10,0,200);
%plot(xmin,fval,'r*')

%{
plot(x,f,'k-.'); hold on;

yt2 = x/5-x.^2/500;
yt6 = x/5-x.^2/500 + x.^4/7500000 - x.^5/937500000 + x.^6/281250000000;
plot(x,yt2,'r'); hold on;
plot(x,yt6,'b'); hold on;
kk = size(x);
plot(ones(kk(2))* pi/2 * alpha, linspace(-1,10,kk(2)))
%histogram(pi/2 * alpha)
%}