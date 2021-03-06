%% Set parameters using Numberical data in the paper

m1=0.300;%kg
m2=0.075; %kg
L1=0.278;%m
L2=0.3%m
l1=0.150;%m
l2=0.148;%m
g=9.81; %m/sec^2
J1=2.48*10^-2; %Kg.m2
J0=J1 + (m1*((l1)^2)) + m2*((L1)^2);
J2_1=3.89*10^-3;  %Kg.m2
J2=J2_1+(m2*((l2)^2));%Kg.m2
b1=1*10^-4; % Nms
b2=2.8*10^-4;%Nms


A1=-b1*J2;
A2=m2*L1*l2*b2;
A3=-((J2)^2);
A4=-0.5*((J2)^2)*m2*L1*l2;
A5=J2*m2*L1*l2;
A6=J2;
A7=-m2*L1*l2;
A8=0.5*((m2)^2)*((L2)^2)*(l1);

B1=m2*L1*l2*b1;
B2=-J0*b2;
B3=-J2*b2;
B4=m2*L1*l2*J2;
B5=0.5*J0*J2;
B6=0.5*((J2)^2);
B7=-0.5*((m2)^2)*((L1)^2)*((l2)^2);
B8=-m2*L1*l2;
B9=J0;
B10=J2;
B11=-m2*l2*J0;
B12=-m2*l2*J2;

C1=-km/Lm;
C2=-Rm/Lm;
C3=1/Lm;
C4=Km;


D1=J0*J2;
D2=(J2)^2;
D3=((m2)^2)*((L1)^2)*((l2)^2);

H1=g*m2^2*l2^2*L1/(D1+D3);
H2=A1/(D1+D3);
H3=-A1/(D1+D3);
H4=g*m2*l2*J0/(D1+D3);
H5=-B1/(D1+D3);
H6=B2/(D1+D3);

H7=A6/(D1+D3);
H8=-B8/(D1+D3);

theta1_allowed_A=(180/pi)*1/180;
theta2_allowed_A=(180/pi)*1/1;


A_MAT=[0 0 1 0;0 0 0 1;0 H1 H2 H3;0 H4 H5 H6];
B_MAT=[0;0;H7;H8];  
R=1;
C_MAT=[1 1 0 0];
Q_Senstivity=[theta1_allowed_A 0;0 theta2_allowed_A];

Q=C'*Q_Senstivity'*C*Q_Senstivity;
F=lqr(A,B,Q,R);


