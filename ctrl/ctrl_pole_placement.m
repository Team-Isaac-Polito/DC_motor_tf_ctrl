%% controller
f = 100;
Ts = 1/f;
z = tf('z',Ts);
h = 0.4362*z/(z^2 -0.5073*z -0.03655); % from sys id
[zeros,poles,k]=zpkdata(h,'v');
B=k;
A=conv([1 -poles(1)],[1 -poles(2)]);

s_hat=0.01; % [-]
t_s1=0.1; % [s]
zeta=abs(log(s_hat)/sqrt(pi^2+(log(s_hat))^2));
wn=4.6/(zeta*t_s1);

% zgrid(zeta,[]), hold on
% plot(z,0,'rx'),
% plot(p,[0,0],'bo'),

Aplus=[1 -poles(1)];
Aminus=[1 -poles(2)];
Bplus=1;
Bminus=B;

p1c=-zeta*wn+wn*sqrt(1-zeta^2)*j;
p2c=-zeta*wn-wn*sqrt(1-zeta^2)*j;
p3c=-10*zeta*wn;

p1=exp(p1c*Ts);
p2=exp(p2c*Ts);
p3=exp(p3c*Ts);

Am=poly([p1 p2 p3]);
Adioph=conv([1 -1],Aminus);
Bdioph=Bminus;

[R1,S1,Am_check]=dioph_mtx(Adioph,Bdioph,Am);
R=conv([1 -1],R1);
S=conv(Aplus,S1);
C=zpk(tf(S,R,Ts)),

L=minreal(C*h,1e-3);
W=zpk(minreal(L/(1+L),1e-4));

%% 2dof controller
Ttilde=[1 -p3];
kW=1;
kT=kW*polyval(S1,1)/(polyval(Ttilde,1)*dcgain(W));
T1=kT*Ttilde;
F=zpk(tf(T1,S1,Ts)),

%% plot

subplot(2,1,1),
step = 50;
r_sim = [0,step*ones(1,30)];
r_sim = lsim(F,r_sim); % for 2dof controller
y_step = lsim(W,r_sim);
plot(0:Ts:Ts*(length(r_sim)-1),y_step,'x-'),
yline(step,'r'),
xlabel('t [s]'),
ylabel('w_output [rpm]')
title('output')

subplot(2,1,2),
Wu=zpk(minreal(C/(1+L),1e-4));
u_step = lsim(Wu,r_sim);
plot(0:Ts:Ts*(length(r_sim)-1),u_step,'x-'),
yline(step,'r'),
xlabel('t [s]'),
ylabel('w_input [rpm]')
title('input')

sgtitle('STEP RESPONSE')