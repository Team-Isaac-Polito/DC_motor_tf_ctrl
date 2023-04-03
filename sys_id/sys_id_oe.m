clc,
clear all,
close all,

csv = importdata('step1.csv');
r1 = csv(1,:)';
y1 = csv(2,:)';
f = 100;
Ts = 1/f;

u = r1(50:200);
y = y1(50:200);
N = length(u);
t = linspace(0,Ts*(N-1),N);

model = oe([y,u],[1,2,1]);
[Num,Den]=tfdata(model);
h = tf(Num,Den,Ts),

% Alternative way
% y = y1;
% u = r1;
% x0 = [0,0,0];
% xopt = fminunc(@(theta) pred_error(theta,y,u),x0);

yp = lsim(h,u);

plot(t,u,'r'), hold on,
plot(t,y,'b'),
plot(t,yp,'k')