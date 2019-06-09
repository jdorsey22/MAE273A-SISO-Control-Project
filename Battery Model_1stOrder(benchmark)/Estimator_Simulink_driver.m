%% MAE 273A Project Script 


clear all

List = ["FirstOrderTruth_BASELINE_linear.mat","FirstOrderTruth_R0_25_linear.mat","FirstOrderTruth_RC_25_linear.mat","FirstOrderTruth_Cc_25_linear.mat","FirstOrderTruth_R0_50_linear.mat","FirstOrderTruth_RC_50_linear.mat","FirstOrderTruth_CC_50_linear.mat","FirstOrderTruth_R0_RC_20_linear.mat","FirstOrderTruth_BASELINE_nonlinear.mat","FirstOrderTruth_R0_25_nonlinear.mat","FirstOrderTruth_RC_25_nonlinear.mat","FirstOrderTruth_CC_25_nonlinear.mat","FirstOrderTruth_R0_50_nonlinear.mat","FirstOrderTruth_RC_50_nonlinear.mat","FirstOrderTruth_CC_50_nonlinear.mat","FirstOrderTruth_R0_RC_20_nonlinear.mat"];
NameList = ["Linear Baseline", "R0:25% Linear", "RC:25% Linear","CC:25% Linear","R0:50% Linear","RC:50% Linear","CC:50% Linear","R0:20% RC: 20%Linear","Nonlinear Baseline ","R0:25% Nonlinear","RC:25% Nonlinear","CC:25% Nonlinear","R0:50% Nonlinear","RC:50% Nonlinear","CC:50% Nonlinear", "R0:20% RC:20% Nonlinear" ];


for k = 1:length(List)
    load(List(k))
    
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

out = sim('Estimator_Simulink_trial1') ;   %run simulink model

est_soc = out.SOC_est;
tout = out.tout; 

YOULA_ACTUAL (:,k) = SOC_act; 
YOULA_ESTIMATED(:,k) = est_soc; 

figure() 
hold on 
plot(t,SOC_act); 
plot(tout,est_soc);

title(NameList(k))
xlabel('Time (seconds) ')
ylabel('State of Charge')

end 

%% Nonlinear Estimator Plant with Nonlinear BASELINE data 


List = ["FirstOrderTruth_BASELINE_linear.mat","FirstOrderTruth_R0_25_linear.mat","FirstOrderTruth_RC_25_linear.mat","FirstOrderTruth_Cc_25_linear.mat","FirstOrderTruth_R0_50_linear.mat","FirstOrderTruth_RC_50_linear.mat","FirstOrderTruth_CC_50_linear.mat","FirstOrderTruth_R0_RC_20_linear.mat","FirstOrderTruth_BASELINE_nonlinear.mat","FirstOrderTruth_R0_25_nonlinear.mat","FirstOrderTruth_RC_25_nonlinear.mat","FirstOrderTruth_CC_25_nonlinear.mat","FirstOrderTruth_R0_50_nonlinear.mat","FirstOrderTruth_RC_50_nonlinear.mat","FirstOrderTruth_CC_50_nonlinear.mat","FirstOrderTruth_R0_RC_20_nonlinear.mat"];


load('OCV_table.mat')
load('OCV_slope_table.mat')
load(List(9));
  
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

out = sim('Estimator_Simulink_trial2') ;   %run simulink model

est_soc = out.SOC_est;
tout = out.tout; 

figure() 
hold on 
plot(t,SOC_act); 
plot(tout,est_soc);


%% Kalman Filter 

List = ["FirstOrderTruth_BASELINE_linear.mat","FirstOrderTruth_R0_25_linear.mat","FirstOrderTruth_RC_25_linear.mat","FirstOrderTruth_Cc_25_linear.mat","FirstOrderTruth_R0_50_linear.mat","FirstOrderTruth_RC_50_linear.mat","FirstOrderTruth_CC_50_linear.mat","FirstOrderTruth_R0_RC_20_linear.mat","FirstOrderTruth_BASELINE_nonlinear.mat","FirstOrderTruth_R0_25_nonlinear.mat","FirstOrderTruth_RC_25_nonlinear.mat","FirstOrderTruth_CC_25_nonlinear.mat","FirstOrderTruth_R0_50_nonlinear.mat","FirstOrderTruth_RC_50_nonlinear.mat","FirstOrderTruth_CC_50_nonlinear.mat","FirstOrderTruth_R0_RC_20_nonlinear.mat"];
NameList = ["Linear Baseline", "R0:25% Linear", "RC:25% Linear","CC:25% Linear","R0:50% Linear","RC:50% Linear","CC:50% Linear","R0:20% RC: 20%Linear","Nonlinear Baseline ","R0:25% Nonlinear","RC:25% Nonlinear","CC:25% Nonlinear","R0:50% Nonlinear","RC:50% Nonlinear","CC:50% Nonlinear", "R0:20% RC:20% Nonlinear" ];


for i = 1:length(List)
    
    load(List(i));

R0 = 0.01; 
Rc = 0.015;
Cc = 2400; 
Cbat = (5*3600); 
alpha = .65;
Voc_0 = 3.435; 
SOC0 = 1; 

dt = .1;

%System Dynamics 
%Linear Model
A = [1 0 ; 0 (1-(dt/(Rc*Cc)))]; 
B = [-(dt/Cbat) ; (dt/Cc)]; 
C = [alpha -1];
D = [-R0]; 

%Kalman Model 
Ak = [1]; 
Bk = [-(dt/Cbat)]; 
Ck=[alpha]; 
Dk=[R0]; 

% Kalman Filter 

wk_mean = 0; 
Q = 2.5*10^-7;
vk_mean = 0; 
R = 1*10^-4;
P0 = 0; 

Q = 2.5*10^-9;


% Set Initial Conditions 
P(1) = P0; 
x1(1) = .98; % SOC
x2(1) = 0;  % Vc
x1_hat(1) = .98; 

for k = 2:1:length(t)
    %State Equations: 
    x1(k) = x1(k-1)-(dt/Cbat)*I(k-1); %+normrnd(0,Q);  % SOC 
    x2(k) = (1-(dt/(Rc*Cc)))*x2(k-1)+(dt/Cc)*I(k-1);  % Vc
    
    % Open Loop: 
    x1_ol(k) = Ak*x1(k-1)+Bk*I(k-1);
    
    % Model Prediction: 
    x1_hat_prev(k) = Ak*x1_hat(k-1)+ Bk*I(k-1);
    P_prev = Ak*P(k-1)*Ak' + Q; 

    %Estimated Output: 
    V_hat = alpha*x1_hat_prev(k)-x2(k)- R0*I(k)+ Voc_0;

    % Measurement Update: 
    x1_hat(k) = x1_hat_prev(k) + P_prev*Ck'*inv(Ck*P_prev*Ck'+R)*(V(k)-V_hat);
    P(k) = P_prev -P_prev*Ck'*inv(Ck*P_prev*Ck'+R)*Ck*P_prev;
    
    
end 

KALMAN_ACTUAL(:,i)= SOC_act; 
KALMAN_ESTIMATED(:,i)=x1_hat; 

figure();
hold on 
plot(t,SOC_act)
plot(t,x1_hat)
title('Kalman Filter '+ NameList(i)); 
xlabel('Time (seconds)'); 
ylabel('State of Charge (SOC)'); 
legend('SOC Act','SOC Est');

end 

%%  Extended Kalman Filter 

List = ["FirstOrderTruth_BASELINE_linear.mat","FirstOrderTruth_R0_25_linear.mat","FirstOrderTruth_RC_25_linear.mat","FirstOrderTruth_Cc_25_linear.mat","FirstOrderTruth_R0_50_linear.mat","FirstOrderTruth_RC_50_linear.mat","FirstOrderTruth_CC_50_linear.mat","FirstOrderTruth_R0_RC_20_linear.mat","FirstOrderTruth_BASELINE_nonlinear.mat","FirstOrderTruth_R0_25_nonlinear.mat","FirstOrderTruth_RC_25_nonlinear.mat","FirstOrderTruth_CC_25_nonlinear.mat","FirstOrderTruth_R0_50_nonlinear.mat","FirstOrderTruth_RC_50_nonlinear.mat","FirstOrderTruth_CC_50_nonlinear.mat","FirstOrderTruth_R0_RC_20_nonlinear.mat"];
NameList = ["Linear Baseline", "R0:25% Linear", "RC:25% Linear","CC:25% Linear","R0:50% Linear","RC:50% Linear","CC:50% Linear","R0:20% RC: 20%Linear","Nonlinear Baseline ","R0:25% Nonlinear","RC:25% Nonlinear","CC:25% Nonlinear","R0:50% Nonlinear","RC:50% Nonlinear","CC:50% Nonlinear", "R0:20% RC:20% Nonlinear" ];


for i = 1:length(List)
    
    load(List(i));

R0 = 0.01; 
Rc = 0.015;
Cc = 2400; 
Cbat = (5*3600); 
alpha = .65;
Voc_0 = 3.435; 
SOC0 = 1; 

dt = .1;

wk_mean = 0; 
Q = 2.5*10^-7;
vk_mean = 0; 
R = 1*10^-4;
P0 = 0; 

A_ek = 1 ;

% function val = C_ek(SOC_hat)
%  val = intep1(soc_intpts_OCV_slope,OCV_slope_intpts,SOC_hat);
% end 
E_ek = 1; 
F_ek = 1; 

Ak = 1;
Bk = -dt/Cbat;

load('OCV_table.mat')
load('OCV_slope_table.mat')
% load('IV_data_nonlinear')

R = 1*10^-4;
Q = 2.5*10^-9;

% Set Initial Conditions 
P(1) = P0; 
x1(1) = .98; % SOC
x2(1) = 0;  % Vc
x1_hat(1) = .98; 

for k = 2:length(t)
    x1(k) = x1(k-1)-(dt/Cbat)*I(k-1); %+normrnd(0,Q);  % SOC - Coulomb Counting 
    x2(k) = (1-(dt/(Rc*Cc)))*x2(k-1)+(dt/Cc)*I(k-1);  % Vc

    
    C_ek =  interp1(soc_intpts_OCV_slope' ,OCV_slope_intpts,x1_hat(k-1));
    
    % Model Prediction: 
    x1_hat_prev = Ak*x1_hat(k-1)+Bk*I(k-1);
    P_prev = A_ek*P(k-1)*A_ek'+E_ek*Q*E_ek'; 
    
    % Measurement Update: 
    V_hat = interp1(soc_intpts_OCV' ,OCV_intpts,x1_hat(k-1)) -I(k)*R0-x2(k);
    L = P_prev*C_ek'*inv(C_ek*P_prev*C_ek'+F_ek*R*F_ek');
    
    x1_hat(k) = x1_hat_prev + L*(V(k)-V_hat); 
    P(k) = P_prev -L*C_ek*P_prev;   
    
end 

EXTENDED_KALMAN_ACTUAL(:,i)= SOC_act; 
EXTENDED_KALMAN_ESTIMATED(:,i)=x1_hat; 

figure();
hold on 
plot(t,SOC_act)
plot(t,x1_hat)
title('Extended Kalman Filter' + NameList(i)); 
xlabel('Time (seconds)'); 
ylabel('State of Charge (SOC)'); 

legend('SOC Act','SOC Est');


end 

%% RMS Tabulation 


for i = 1:length(List)


Youla_act_rms(1,i) = rms(YOULA_ACTUAL(:,i));
Youla_est_rms(1,i) = rms(YOULA_ESTIMATED(:,i));

Kalman_act_rms(1,i) = rms(KALMAN_ACTUAL(:,i));  
Kalman_est_rms(1,i) = rms(KALMAN_ESTIMATED(:,i)); 

EKF_act_rms(1,i) = rms(EXTENDED_KALMAN_ACTUAL(:,1)); 
EKF_est_rms(1,i) = rms(EXTENDED_KALMAN_ESTIMATED(:,1)); 


table = [ "",NameList;...
         "Youla Actual RMS", Youla_act_rms;...
         "Youla Estimated RMS", Youla_est_rms;...
         "Youla RMS Error", (Youla_act_rms-Youla_est_rms);... 
         
         "Kalman Filter Actual RMS", Kalman_act_rms;...
         "Kalman Filter Estimated RMS",Kalman_est_rms;...
         "Kalman RMS Error", (Kalman_act_rms-Kalman_est_rms);... 

         "EKF Actual RMS", EKF_act_rms; ...
         "EKF Estimated RMS", EKF_est_rms;...
         "EKF RMS Error", (EKF_act_rms-EKF_est_rms);... 
]; 

end 

