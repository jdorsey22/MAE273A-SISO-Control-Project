%% Hinf Controller
clearvars
close all
clc

s=tf('s');

Gp = (-64800*s^2 - 2093.4*s - 0.65)/(18000*s*(36*s+1));  %plant
sysg=ss(Gp);
[Ag,Bg,Cg,Dg]=ssdata(sysg);

Ag=Ag-0.08*eye(1); % slightly shift A to avoid poles on jw axis

[num,den]=ss2tf(Ag,Bg,Cg,Dg);
Gpn=tf(num,den); % perturbed plant

% %Hinf shaping filter
W1=(s+43)/(2*s+0.01); %(BW >= 65rad/s(CL at -6dB) so step input resp shows tracking at high freq) *input data at 10Hz
W2=0.05;
% W3=0.5;  
W3 = makeweight(1/2,43,50);  %(low freq gain, cross over freq, high freq gain) 
                             %^^to make it look like a differentiator!
                           
%Hinf Controller Computation
ssga_=augtf(Gpn,W1,W2,W3);
[sys3,sscl,GAM]=hinfsyn(ssga_);  %sys3=controller, sscl=CL tf (w/z), 
  
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

% figure(4)
w=logspace(-6,3);
% % sigma(GAM*sscl,w,'y'), grid
% sigma(sscl,ss(GAM)), grid, legend('sscl','GAM')

figure(4)
NP = W1*S;
[numNP, denNP]=tfdata(NP,'v');
[magNP,phaseNP]=bode(numNP,denNP,w);
lmNP=20*log10(magNP);

RS = W3*T;
[numRS, denRS]=tfdata(RS,'v');
[magRS,phaseRS]=bode(numRS,denRS,w);
semilogx(w,magNP,'b',w,magRS,'r'), grid, legend('W1*S','W3*T')

figure(5)
bode(S,T,W1,'--',W3,'-.'), legend('S','T','W1','W3')

figure(6)
w=logspace(-4,4);
sigma(GAM*sscl,w,'y'),grid   %calculate gamma*sscl
check=GAM/sqrt(2);
lmcheck=20*log10(check);
lmcheck=lmcheck*ones(size(w));
hold on
semilogx(w,lmcheck,'r')
TRP1=-norm(lmcheck,inf)-20*log10(norm(sigma(GAM*sscl),inf));   %Gamma*sscl - Gamma/sqrt(2)

figure(7)
RP=magRS+magNP;
magRP=norm(RP,inf);
lmRP=20*log10(magRP);
magRP=magRP*ones(size(w));
semilogx(w,magRP,'r'),grid
TRP2=0-norm(lmRP,inf);  %difference between 1(0 in dB) and ||WpS| + |WdT|| (we will fail meet this)

figure(8)
bode(T,S)