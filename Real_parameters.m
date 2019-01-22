
J_M = 6.75e-6;          % Inertia Motor 6.75e-6
J_Enc1 = 6e-14;         % Inertia Encoder Theta1
J_Enc2 = 0.1e-6;       % Inertia Encoder Theta2
J_m_b = 3.8281e-6;      % Inertia m_b about theta2 (= 1/2 * m_b *r_b^2)

r_b=0.035/2;

m_sens = 0.084;
m_horzArm=0.263;
m_r = m_horzArm - m_sens;
m_b=0.025;

l_a=.177;
l_b=.007;
l_c = 0.01;
l_r=0.29;

% Use cylinder for m_a:
m_a=.260;
r_a=.015;

m_pm=0.0383;
m_p=0.0107;
l_pm=.198;
l_p=.233;

l_s=l_r+l_b-l_a+(l_c/2);

m1=m_sens+m_r+m_b;
m2=m_pm+m_p;
Jsh=m_sens*(l_s)^2;
Jmah=0.5*m_a*(r_a)^2
Jmrh=((1/12)*(m_r)*(l_r)^2)+(((m_r)*((-l_r/2)+l_a-l_b)^2));
Jmbh=(m_b*(l_a)^2);
J1h=Jsh+Jmrh+Jmbh+J_Enc1;
Jmph=(1/3)*m_p*(l_pm+r_b)^2    %% Ignored
Jmpmh=m_pm*(l_p)^2;
J2=J_m_b+Jmpmh+J_Enc2  %% J2h
J0=J_M+J1h+(m2*(l_a)^2)+Jmah
L1=sqrt(J1h/m1);
l2=sqrt(J2/m2);
b1=1*10^-4;
b2=2.8*10^-4;
g=9.81;
xi=[0 0.1 0 0];

km=0.09;
Rm=7.8;
Lm=0.005;

A1=-b1*J2;
A2=m2*L1*l2*b2;
A3=-((J2)^2);
A4=-0.5*((J2)^2)*m2*L1*l2;
A5=J2*m2*L1*l2;
A6=J2;
A7=-m2*L1*l2;
A8=0.5*((m2)^2)*((l2)^2)*(L1);

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
C3=1/Lm; % %%%
C4=km;


D1=J0*J2;
D2=(J2)^2;
D3=((m2)^2)*((L1)^2)*((l2)^2);
%% A, B and C Matrix for the Model

A31=0;
A32=g*m2^2*l2^2*L1/(D1+D3);
A33=A1/(D1+D3);
A34=-b2*m2*l2*L1/(D1+D3);
A41=0;
A42=g*m2*l2*J0/(D1+D3);
A43=-B1/(D1+D3);
A44=B2/(D1+D3);

B31=J2/(D1+D3);
B32=m2*L1*l2/(D1+D3);
B41=m2*L1*l2/(D1+D3);
B42=J0/(D1+D3);



H7=A6/(D1+D3);
H8=-B8/(D1+D3);

theta1_allowed_A=(180/pi)*1/180;
theta2_allowed_A=(180/pi)*1/.1;
E0=0;
K_swing=-100;

%%A_MAT=[0 0 1 0;0 0 0 1;A31 A32 A33 A34;A41 A42 A43 A44]
A_MAT=[0 0 1 0;0 0 0 1;0 24.4662 -0.0062 -0.0664;0 65.7367 -0.0051 -0.17]
B_MAT=[0;0;3.2;3.1];
R=0.01;
C_MAT=[1 0 0 0;0 1 0 0];
Q_Senstivity=[theta1_allowed_A 0;0 theta2_allowed_A];


Q=C_MAT'*(Q_Senstivity)'*Q_Senstivity*C_MAT;
Q=[.1 0 0 0;0 100000 0 0;0 0 .1 0;0 0 0 0.1];
%%Q=10*eye(4);
%%F_MAT=-lqr(A_MAT,B_MAT,Q,R);
F_MAT=-[0.16 -320 1.3 -16];
Q1=[1000 0 0 0;0 10 0 0;0 0 1 0;0 0 0 1];
F_MAT2=-lqr(A_MAT,B_MAT,Q1,R);
Open_Loop_Poles=eig(A_MAT);
Closed_Loop_Poles=eig(A_MAT+(B_MAT*F_MAT));
L_MAT=place(A_MAT',C_MAT',[12*Closed_Loop_Poles(1) 12*Closed_Loop_Poles(2) 10*Closed_Loop_Poles(3) 10*Closed_Loop_Poles(4)])';

