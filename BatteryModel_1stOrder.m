% syms s Rc Cc Cbat alpha R0
% syms s
s = tf('s');
Rc =0.015; 
Cc = 2400;
Cbat = 5*3600;
alpha =0.65;
R0=0.01;

K = -18000/.65;   %gain
zeta = 0.707;     %damping ratio
wn = 0.0005;         %natural frequency

%continuous time ss model
A = [-1/(Rc*Cc) 0; 0 0] ;
B = [1/Cc; -1/Cbat];
C = [-1 alpha];
D = -R0;
SI = [s 0;0 s];
% SI_inv = (SI - A)^-1
% C_SIinv = C*SI_inv
% C_SIinv_B =C_SIinv*B
% final= C_SIinv_B+D
% Gp = -R0 - Rc/(Rc*Cc*s+1) - alpha/(Cbat*s) %hand computed

Gp = C*(SI-A)^-1*B+D; %plant
Y = K*wn^2*s/(s^2+2*zeta*wn*s+wn^2); %youla
T = minreal(Gp*Y); %complimentary
S = minreal(1-T);  %sensitivity
Gc = Y/(1-Y*Gp);   %controller
L = Gc*Gp;         %open loop TF
sysTF = Gc*Gp/(1+Gc*Gp); %actual sys TF

figure(1)
bode(Y,T,S), legend('Y','T','S')
figure(2)
step(sysTF), stepinfo(sysTF)

