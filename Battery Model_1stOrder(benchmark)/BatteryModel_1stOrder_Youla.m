%Load files for Simulink Model
load('IV_data_linear.mat')
s = tf('s');

%battery model parameters
Rc =0.015; 
Cc = 2400;
Cbat = 5*3600;
alpha =0.65;
R0=0.01;
Vocv0 = 3.435; %V

%tunning parameters
K = 1;   %gain
zeta = 0.707;     %damping ratio
wn = 0.01;      %natural frequency

%continuous time ss model
A = [-1/(Rc*Cc) 0; 0 0];
B = [1/Cc; -1/Cbat];
C = [-1 alpha];
D = -R0;

SI = [s 0;0 s];
Gp = C*(SI-A)^-1*B+D; %plant
T = minreal(K*wn^2/(s^2+2*zeta*wn*s+wn^2)); %complimentary
Y = minreal(T/Gp); %youla
S = minreal(1-T);  %sensitivity
Gc = Y/S;   %controller
L = Gc*Gp;         %open loop TF
sysTF = Gc*Gp/(1+Gc*Gp); %actual sys TF
% [num,den] = tfdata(Gp) 
[num den] = tfdata(Gc, 'v');
figure(1)
bode(Y,T,S), legend('Y','T','S'), grid on
figure(2)
step(sysTF),legend('System Step Response'), grid on, stepinfo(sysTF)
