%% Fifth Order Simulation Truth Generation Script: 


%%

%% Set Parameters: 

clear all 

%battery model parameters
Rc = 0.015;    %Ohms
Cc = 2400;     %F
Cbat = 5*3600;
alpha =0.65;   
R0 = 0.01;     %Ohms
Vocv0 = 3.435; %V

% R0 = R0*1.25;
% Rc = Rc*1.25;
% Cc = Cc*1.25;

% R0 = R0*1.5;
% Rc = Rc*1.5;
% Cc = Cc*1.5;


% R0 = R0*1.20;
% Rc = Rc*1.20;


%tunning parameters
K = 1;         %gain
zeta = 0.87;  %damping ratio
wn = 60;      %natural frequency

%continuous time ss model
A = [-1/(Rc*Cc) 0; 0 0];
B = [1/Cc; -1/Cbat];
C = [-1 alpha];
D = -R0;


Tau1 = Cc*Rc; 

dt = .1; 

% System Dynamics

% Linear State Dynamics: Dual Polarity Model 

% Continuous Time Model: 
A_c = [0       0         ; ...
     0  (-1/(Rc*Cc))     ]; 
B_c = [(-1/Cbat); (1/Cc)]; 
C_c = [alpha -1 ];
D_c = [-R0]; 


Ad = [1      0        ; ...
     0 exp(-dt/Tau1) ]; 
Bd = [(-dt/Cbat); (Rc)*(1-exp(-dt/Tau1))]; 
Cd = C_c; 
Dd = D_c; 


% Load Battery Measurements 
load('OCV_table.mat')
load('OCV_slope_table.mat')
load('IV_data_nonlinear.mat')
% load('IV_data_linear.mat')


%% State/Output Simulation with Process/Measurement Noise (Truth) 

P(1) = 0;           % Covariance 
x1(1) = .98;          % SOC - Battery Fully Charged 
x2(1) = 0;          % Vc1
% x3(1) = 0;          % Vc2


for k = 2:1:length(t)
    
    x1(k) = Ad(1,1)*x1(k-1) + Bd(1,1)*I(k-1)+ normrnd(0,.00045); % soc
    x2(k) = Ad(2,2)*x2(k-1) + Bd(2,1)*I(k-1) +normrnd(0,.0002); % Vc1
%     x3(k) = Ad(3,3)*x3(k-1) + Bd(3,1)*I(k-1)+ normrnd(0,.0002); % Vc2
    
    V_truth(k) = interp1(soc_intpts_OCV',OCV_intpts,x1(k-1)) - I(k-1)*R0 - x2(k-1)+ normrnd(0,.05);
%     V_truth(k) = Vocv0 + alpha*x1(k) - I(k-1)*R0 - x2(k-1)+ normrnd(0,.05);

end 
%%
figure()
hold on 
plot(t,x1)
plot(t,SOC_act)
legend('Simulated Truth','Lin_SOC_act');

figure(); 
plot(t,V_truth)

%% NaN Determination for Preventing Interpolation Fuckups

% thing = isnan(V_truth);
thing = isnan(x1);
% thing = isnan(V);


counter =0; 
for k=1:length(t)
    if thing(k)==1
       display('yes') 
       
       counter = counter +1;
       
    end     
end 

counter



%% NAN Correction 

for u = 1:length(t)
    if thing(u) == 1
        V_truth(u)  = 0; 
        
    end
    
    
    
end 

%% FINAL SOC Correction: 

for O = 1:length(t)
    if SOC_act(O) >1
        SOC_act(O) = 1; 
        
    end 
    
    
end 




%% Relabel New truth Data

clear V SOC_act 

V = V_truth; 
SOC_act = x1;

V = V'; 
SOC_act = SOC_act'; 

%%

save('FirstOrderTruth_BASELINE_nonlinear.mat','V','SOC_act','t','I'); 


%% 
hold on 
plot(t,SOC_act)




%%
    List = ["FirstOrderTruth_BASELINE_linear.mat","FirstOrderTruth_R0_25_linear.mat","FirstOrderTruth_RC_25_linear.mat","FirstOrderTruth_Cc_25_linear.mat","FirstOrderTruth_R0_50_linear.mat","FirstOrderTruth_RC_50_linear.mat","FirstOrderTruth_CC_50_linear.mat","FirstOrderTruth_R0_RC_20_linear.mat","FirstOrderTruth_BASELINE_nonlinear.mat","FirstOrderTruth_R0_25_nonlinear.mat","FirstOrderTruth_RC_25_nonlinear.mat","FirstOrderTruth_CC_25_nonlinear.mat","FirstOrderTruth_R0_50_nonlinear.mat","FirstOrderTruth_RC_50_nonlinear.mat","FirstOrderTruth_CC_50_nonlinear.mat","FirstOrderTruth_R0_RC_20_nonlinear.mat"];


for j = 1:length(List)
    
    clearvars -except List j



%battery model parameters
Rc = 0.015;    %Ohms
Cc = 2400;     %F
Cbat = 5*3600;
alpha =0.65;   
R0 = 0.01;     %Ohms
Vocv0 = 3.435; %V

switch j
    
    case 1 
        % Do Nothing 
        
    case 2 
        R0 = R0*1.25;
    case 3 
        Rc = Rc*1.25;

    case 4
        Cc = Cc*1.25;

    case 5 
        R0 = R0*1.5;

    case 6 
        Rc = Rc*1.5;

    case 7 
        Cc = Cc*1.5;

    case 8 
        R0 = R0*1.20;
        Rc = Rc*1.20;

    case 9 
        % Do nothing 

    case 10 
        R0 = R0*1.25;

    case 11
        Rc = Rc*1.25;
        
    case 12
        Cc = Cc*1.25;
        
    case 13
        R0 = R0*1.5;
        
    case 14
        Rc = Rc*1.5;
    
    case 15 
        Cc = Cc*1.5;
        
    case 16 
        R0 = R0*1.20;
        Rc = Rc*1.20;
    otherwise 
        % Do nothing 
    
% R0 = R0*1.25;
% Rc = Rc*1.25;
% Cc = Cc*1.25;
% 
% R0 = R0*1.5;
% Rc = Rc*1.5;
% Cc = Cc*1.5;
% 
% 
% R0 = R0*1.20;
% Rc = Rc*1.20;

end 


%tunning parameters
K = 1;         %gain
zeta = 0.87;  %damping ratio
wn = 60;      %natural frequency

%continuous time ss model
A = [-1/(Rc*Cc) 0; 0 0];
B = [1/Cc; -1/Cbat];
C = [-1 alpha];
D = -R0;


Tau1 = Cc*Rc; 

dt = .1; 

% System Dynamics

% Linear State Dynamics: Dual Polarity Model 

% Continuous Time Model: 
A_c = [0       0         ; ...
     0  (-1/(Rc*Cc))     ]; 
B_c = [(-1/Cbat); (1/Cc)]; 
C_c = [alpha -1 ];
D_c = [-R0]; 


Ad = [1      0        ; ...
     0 exp(-dt/Tau1) ]; 
Bd = [(-dt/Cbat); (Rc)*(1-exp(-dt/Tau1))]; 
Cd = C_c; 
Dd = D_c; 


% Load Battery Measurements 
load('OCV_table.mat')
load('OCV_slope_table.mat')

if j <=8 
    load('IV_data_linear.mat')
elseif j>=9
    load('IV_data_nonlinear.mat')   
end 


P(1) = 0;           % Covariance 
x1(1) = .98;          % SOC - Battery Fully Charged 
x2(1) = 0;          % Vc1


for k = 2:1:length(t)
    
    x1(k) = Ad(1,1)*x1(k-1) + Bd(1,1)*I(k-1)+ normrnd(0,.00045); % soc
    x2(k) = Ad(2,2)*x2(k-1) + Bd(2,1)*I(k-1) +normrnd(0,.0002); % Vc1
    if j <=8 
    V_truth(k) = Vocv0 + alpha*x1(k) - I(k-1)*R0 - x2(k-1)+ normrnd(0,.05);

    elseif j>=9
    V_truth(k) = interp1(soc_intpts_OCV',OCV_intpts,x1(k-1)) - I(k-1)*R0 - x2(k-1)+ normrnd(0,.05);

    end 

end 

    
clear V SOC_act 

V = V_truth; 
SOC_act = x1;

V = V'; 
SOC_act = SOC_act'; 



save(List(j),'V','SOC_act','t','I'); 
    
    
    
    
end 