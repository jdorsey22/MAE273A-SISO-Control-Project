%Load files for Simulink Model
load('IV_data_linear.mat')
s = tf('s');

%battery model parameters
Rc = 0.015;    %Ohms
Cc = 2400;     %F
Cbat = 5*3600;
alpha =0.65;   
% R0 = 0.01;     %Ohms
R0 = .1;     %Ohms

Vocv0 = 3.435; %V

%tunning parameters
K = 1;         %gain
% zeta = 0.707;  %damping ratio
zeta = 0.25;  %damping ratio

wn = 5;      %natural frequency

%continuous time ss model
A = [-1/(Rc*Cc) 0; 0 0];
B = [1/Cc; -1/Cbat];
C = [-1 alpha];
D = -R0;

A1 = A(1,1);
B1 = B(1,1);
A2 = A(2,2); 
B2 = B(2,1); 
C2 = C(1,2); 

SI = [s 0;0 s];
Gp = C*(SI-A)^-1*B+D;         %plant
T = minreal(K*wn^2/(s^2+2*zeta*wn*s+wn^2)); %complimentary
Y = minreal(T/Gp);            %youla
S = minreal(1-T);             %sensitivity
Gc = Y/S;                     %controller
L = Gc*Gp;                    %open loop TF
sysTF = Gc*Gp/(1+Gc*Gp);      %actual sys TF

[num, den] = tfdata(Gc, 'v'); %get numerator and denominator of Gc tf
%%
figure(1)
bode(Y,T,S), legend('Y','T','S'), grid on
figure(2)
step(sysTF),legend('System Step Response'), grid on, stepinfo(sysTF)
%%
sim('Estimator_Simulink')    %run simulink model