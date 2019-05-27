%% KF Battery Model - Felipe Valdez
% clear; clc; close all;

%%%%%%%%%%%%%%% Current Input Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load IV_data_linear

%Battery parameters
C_bat =5*3600; 
R0 = 0.01; %Ohms
Rc = 0.015; %Ohms
Cc = 2400; %F
alpha = 0.65; %V
Vocv0 = 3.435; %V

%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%
Vc = zeros(1,length(t));       %prelocate array
SOC = zeros(1,length(t));     
P = zeros(1,length(t));
SOC_ol = zeros(1,length(t));   
fx = zeros(1,length(t));       

Vc(1) = 0;
SOC(1) = 1;
OCV(1) = Vocv0 + alpha*SOC(1); 
V(1) = OCV(1); %at t=0, V(t)=OVC(t) since there is no V drop
mu = 0;                        %zero mean
Q = 2.5e-7;                    %system noise covariance
R = 10^-4;                     %measurement(voltage) noise covariance
SOC_ol(1) = SOC(1);
dt = 0.1;                      %sampling period

%%%%%%%%%%%%%%%%%%%%%%%%% Kalman Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 2:length(t)
% w = normrnd(mu,Q);
% v = normrnd(mu,R);

%State space 
A = 1;  B = -dt/C_bat;  C = alpha;  D = -R0;

%Open Loop
SOC_ol(k) = A*SOC_ol(k-1)+B*I(k);

Vc(k) = (1-(dt/(Rc*Cc)))*Vc(k-1) + (dt/Cc)*I(k-1);

%%%%%%%%%%%%%%%%%%%%% MODEL PREDICTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
SOC_prev = A*SOC(k-1) + B*I(k-1); %Step 1a: State estimate time update
P_prev = A*P(k-1)*A'+ Q;  %Step 1b: Error covariance time update 
V_est = alpha*SOC_prev - Vc(k) - R0*I(k) + Vocv0; %Step 1c: Estimate system output 

%%%%%%%%%%%%%%%%%%%%% MEASUREMENT UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = P_prev*C'*inv(C*P_prev*C'+R); %Step 2a: Compute Kalman gain 
SOC(k) = SOC_prev + L*(V(k) - V_est); %Step 2b: State estimate measurement update 
P(k) = P_prev - L*C*P_prev; %Step 2c: Error covariance measurement update 

end

figure; plot(t,SOC,t,SOC_act,t,SOC_ol), legend('estimated', 'actual', 'open loop')
xlabel('time'), ylabel('SOC'), title('Felipe Valdez')

%---------------------------------Part c-----------------------------------
syms P_inf

%Algebraic Riccatti equation
R_eqn = A*P_inf*A.' + Q - A*P_inf*C.'*(C*P_inf*C.'+R)^-1*C*P_inf*A.'-P_inf;
sol = solve(R_eqn==0,P_inf); %solve Riccatti eqn
P_inf = eval(sol(2)); %evaluate Riccatti eqn and get solution 2

P_inf_prev = A*P_inf*A.' + Q - A*P_inf*C.'*(C*P_inf*C.'+R)^-1*C*P_inf*A.';

%compare P_inf Riccatti with KF 
Pk_inf = P_inf_prev - P_inf_prev*C.'*(C*P_inf_prev*C.'+R)^-1*C*P_inf_prev
Pk_inf_KF = P(end)

%--------------------------------Part d------------------------------------
e = SOC_act - SOC';
sigma = sqrt(P(end));

interval = linspace(-0.02,0.02,3500);
[f,x] = hist(e,interval); %use hist function and get unnormalized values
fx_pdf = pdf('norm',x,0,sigma); 

figure; plot(x,f,x,fx_pdf,'r');
xlabel('error'), ylabel('probability density')
legend('SOC estimation error','N(0,P)'), title('Felipe Valdez')
figure, plot(t,e), xlabel('time(s)'), ylabel('error'), title('Felipe Valdez')