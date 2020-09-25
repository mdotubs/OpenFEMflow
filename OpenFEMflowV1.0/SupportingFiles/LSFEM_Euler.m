clear all
close all
clc

%% Defining symbols

syms L1 L2
syms x y x1 x2 x3 y1 y2 y3 
syms r u v p
syms r1 u1 v1 p1 r2 u2 v2 p2 r3 u3 v3 p3
syms r0 u0 v0 p0
syms drdx dudx dvdx dpdx drdy dudy dvdy dpdy
syms dt

%% coordinate definition

L3 = 1 - L1 - L2;


%% Defining shape functions - three nodes element for pressure

S1 = L1;
S2 = L2;
S3 = L3;


%% derivatives

dx_dL1 = x1 - x3;
dx_dL2 = x2 - x3;
dy_dL1 = y1 - y3;
dy_dL2 = y2 - y3;

J = [dx_dL1  dy_dL1;
     dx_dL2  dy_dL2];

% Ji = inv(J);
detJ = (x1 - x3)*(y2 - y3) - (x2 - x3)*(y1 - y3);

dS1 = J\[diff(S1,L1); diff(S1,L2)];  dS1dx = dS1(1);  dS1dy = dS1(2);
dS2 = J\[diff(S2,L1); diff(S2,L2)];  dS2dx = dS2(1);  dS2dy = dS2(2);
dS3 = J\[diff(S3,L1); diff(S3,L2)];  dS3dx = dS3(1);  dS3dy = dS3(2);

%% INS equation

g = 1.4;


A1 = [u  r  0  0
      0  u  0  1/r  
      0  0  u  0
      0  g*p 0  u];
  
A2 = [v  0  r  0  
      0  v  0  0
      0  0  v  1/r
      0  0  g*p  v];
  
E = eye(4);


Ln{1} = E*S1 + dt*A1*dS1dx + dt*A2*dS1dy;
Ln{2} = E*S2 + dt*A1*dS2dx + dt*A2*dS2dy;
Ln{3} = E*S3 + dt*A1*dS3dx + dt*A2*dS3dy;

L = [Ln{1} Ln{2} Ln{3}];

% f = -dt*(A1*[drdx;dudx;dvdx;dpdx] +A2*[drdy;dudy;dvdy;dpdy]);
f = [r; u; v; p];

%% Matrices

SK = transpose(L)*L;
K = int(int(SK*detJ,L2,0,1-L1),L1,0,1);    

SF = transpose(L)*f;
F = int(int(SF*detJ,L2,0,1-L1),L1,0,1);   


%%

U = [r1; u1; v1; p1; r2; u2; v2; p2; r3; u3; v3; p3];
R = K*U-F;

dRdr = diff(R,r);
dRdu = diff(R,u);
dRdv = diff(R,v);
dRdp = diff(R,p);

dRdr1 = diff(R,r1) + dRdr*S1;
dRdr2 = diff(R,r2) + dRdr*S2;
dRdr3 = diff(R,r3) + dRdr*S3;

dRdu1 = diff(R,u1) + dRdu*S1;
dRdu2 = diff(R,u2) + dRdu*S2;
dRdu3 = diff(R,u3) + dRdu*S3;

dRdv1 = diff(R,v1) + dRdv*S1;
dRdv2 = diff(R,v2) + dRdv*S2;
dRdv3 = diff(R,v3) + dRdv*S3;

dRdp1 = diff(R,p1) + dRdp*S1;
dRdp2 = diff(R,p2) + dRdp*S2;
dRdp3 = diff(R,p3) + dRdp*S3;

dRdU1 = [dRdr1 dRdu1 dRdv1 dRdp1];
dRdU2 = [dRdr2 dRdu2 dRdv2 dRdp2];
dRdU3 = [dRdr3 dRdu3 dRdv3 dRdp3];


%% Error for mesh refinement

Ln{1} = A1*dS1dx + A2*dS1dy;
Ln{2} = A1*dS2dx + A2*dS2dy;
Ln{3} = A1*dS3dx + A2*dS3dy;

L = [Ln{1} Ln{2} Ln{3}];

SE = transpose(L)*L;
E = int(int(SE*detJ,L2,0,1-L1),L1,0,1);    
E = transpose(U)/detJ*E*U;

