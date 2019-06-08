
clear all
%Load files for Simulink Model
% load('IV_data_linear.mat')
load('FirstOrderTruth_BASELINE_linear.mat')

s = tf('s');

%battery model parameters
Rc = 0.015;    %Ohms
Cc = 2400;     %F
Cbat = 5*3600;
alpha =0.65;   
R0 = 0.01;     %Ohms
Vocv0 = 3.435; %V

%tunning parameters
K = 1;         %gain
% zeta = 0.707;  %damping ratio
zeta = 0.5;  %damping ratio

wn = 75;      %natural frequency

%continuous time ss model
A = [-1/(Rc*Cc) 0; 0 0];
B = [1/Cc; -1/Cbat];
C = [-1 alpha];
D = -R0;

A1 = A(1,1);
B1 = B(1,1);
A2 = A(2,2); 
B2 = B(2,1); 
C2 = alpha;

SI = [s 0;0 s];
Gp = C*(inv(SI-A))*B+D;       %plant
T = minreal(K*wn^2/(s^2+2*zeta*wn*s+wn^2)); %complimentary
Y = minreal(T/Gp);            %youla
S = minreal(1-T);             %sensitivity
Gc = minreal(Y/S);            %controller
L = minreal(Gc*Gp);                    %open loop TF
sysTF = minreal(Gc*Gp/(1+Gc*Gp));      %actual sys TF
[num, den] = tfdata(Gc, 'v'); %get numerator and denominator of Gc tf
% figure(1)
% bode(Y,T,S), legend('Y','T','S'), grid on
% figure(2)
% step(sysTF),legend('System Step Response'), grid on, stepinfo(sysTF)

% sys = ss(A,B,C,D);  %state space model 
% [num,den] = ss2tf(A,B,C,D); %get open loop tf=num/den from state space model
% G = tf(num,den)
% H = G*Gc
% Hcl = feedback(H,1);


%%
out = sim('Estimator_Simulink_trial1') ;   %run simulink model

%%
est_soc = out.SOC_est;
tout = out.tout; 
%%

plot(tout, est_soc)
hold on 
plot(t,SOC_act)


%%


est_soc = est_soc(1:end-1); 

tout = out.tout; 

tout_d = tout(1:end-1);

%%

hold on
plot(t,SOC_act)
plot(tout_d, est_soc)

