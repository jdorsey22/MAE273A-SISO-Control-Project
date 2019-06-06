%% Hinf example problem
clearvars
close all
clc

s=tf('s');

%% example 1 
% Gp=(10+s)/s^3;
% Gp=1/s;
Gp = (-64800*s^2 - 2093.4*s - 0.65)/(18000*s*(36*s+1));
sysg=ss(Gp);
[Ag,Bg,Cg,Dg]=ssdata(sysg);

Ag=Ag-0.08*eye(1); % slightly shift A to avoid poles on jw axis

[num,den]=ss2tf(Ag,Bg,Cg,Dg);
Gpn=tf(num,den); % perturbed plant

% %Hinf shaping filter
W1=(s+2)/(2*s+0.001); %have a bandwidth of 1Hz (made it 10Hz so step input response shows tracking)
W2=0.01;%1;
W3=0.5;  
% W3 = makeweight(1/2,1,10);

%Hinf shaping filter
% W1=(s+100)/(2*s+0.001);
% W2=0.0001;%1;
% W3=s/100;
% %Hinf shaping filter
% W1=1/(8*s+1);
% W2=0.5;%1;
% W3=[];
  
%Hinf Controller Computation
ssga_=augtf(Gpn,W1,W2,W3);
[sys3,sscl,GAM]=hinfsyn(ssga_); % 
  
% Hinf Controller
Gc=minreal(tf(sys3));

% results
L=minreal(Gc*Gp);
Y=minreal(Gc/(1+Gc*Gp));
T=minreal(Y*Gp);
S=minreal(1-T);

% plot results
figure(1)
step(T)

figure(2)
bodemag(T)
hold on
bodemag(S)
bodemag(1/W1)
legend('complementary sensitivity','sensitivity','invW1','location','southeast')
hold off

figure(3)
bode(W1)
hold on
bode(S)
title('frequency response')
legend('W1','sensitivity','location','southeast')

figure(4)
w=logspace(-6,3);
% sigma(GAM*sscl,w,'y'), grid
sigma(sscl,ss(GAM)), grid, legend('sscl','GAM')

figure(5)
NP = W1*S;
[numNP, denNP]=tfdata(NP,'v');
[magNP,phaseNP]=bode(numNP,denNP,w);
lmNP=20*log10(magNP);

RS = W3*T;
[numRS, denRS]=tfdata(RS,'v');
[magRS,phaseRS]=bode(numRS,denRS,w);
semilogx(w,magNP,'b',w,magRS,'r'), grid, legend('W1*S','W3*T')

figure(6)
% bode(S,T,W1,'--'), legend('S','T','W1')
bode(S,T,W1,'--',W3,'-.'), legend('S','T','W1','W3')

%% example 2 
Gp=(1-s)/s^3;

sysg=ss(Gp);
[Ag,Bg,Cg,Dg]=ssdata(sysg);

Ag=Ag-0.08*eye(3); % slightly shift A to avoid poles on jw axis

[num,den]=ss2tf(Ag,Bg,Cg,Dg);
Gpn=tf(num,den); % perturbed plant

%Hinf shaping filter
W1=1/(8*s+1);
W2=0.5;%1;
W3=[];
  
%Hinf Controller Computation
ssga_=augtf(Gpn,W1,W2,W3);
[sys3,sscl]=hinfsyn(ssga_); % 
  
% Hinf Controller
Gc=minreal(tf(sys3));

% results
L=minreal(Gc*Gp);
Y=minreal(Gc/(1+Gc*Gp));
T=minreal(Y*Gp);
S=minreal(1-T);

% plot results
figure(4)
step(T)

figure(5)
bodemag(T)
hold on
bodemag(S)
bodemag(1/W1)
legend('complementary sensitivity','sensitivity','invW1')
hold off

figure(6)
bodemag(W1)
hold on
bodemag(S)
title('frequency response')
legend('W1','sensitivity')

%% example 3 
Gp=400/(s^2+2*s+400);

%Hinf shaping filter
W1=(s+350.8)/(s+0.001);
W2=[];
W3=[];
  
%Hinf Controller Computation
ssga_=augtf(Gp,W1,W2,W3);
[sys3,sscl]=hinfsyn(ssga_); % 
  
% Hinf Controller
Gc=minreal(tf(sys3));

% results
L=minreal(Gc*Gp);
Y=minreal(Gc/(1+Gc*Gp));
T=minreal(Y*Gp);
S=minreal(1-T);

% plot results
figure(7)
step(T)

figure(8)
bodemag(T)
hold on
bodemag(S)
bodemag(1/W1)
legend('complementary sensitivity','sensitivity','invW1')
hold off

figure(9)
bodemag(W1)
hold on
bodemag(S)
title('frequency response')
legend('W1','sensitivity')