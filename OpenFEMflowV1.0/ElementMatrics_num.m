function [R, K] = ElementMatrics_num(X,Y,U,dt,CFL)

[SK1, SF1] =  IntegralArg(X,Y,U,dt,CFL,1/6,1/6);
[SK2, SF2] =  IntegralArg(X,Y,U,dt,CFL,2/3,1/6);
[SK3, SF3] =  IntegralArg(X,Y,U,dt,CFL,1/6,2/3);


K = 1/6 * (SK1 + SK2 + SK3);
F = 1/6 * (SF1 + SF2 + SF3);

R = K*U - F;

end

function [SK, SF] = IntegralArg(X,Y,U,dt,CFL,L1,L2)
x1 = X(1);
x2 = X(2);
x3 = X(3);
y1 = Y(1);
y2 = Y(2);
y3 = Y(3);

r1 = U(1);  u1 = U(2);  v1 = U(3);  p1 = U(4); 
r2 = U(5);  u2 = U(6);  v2 = U(7);  p2 = U(8);  
r3 = U(9);  u3 = U(10); v3 = U(11); p3 = U(12);


L3 = 1 - L1 - L2;

S1 = L1;
S2 = L2;
S3 = L3;

u = u1*S1 + u2*S2 + u3*S3;
v = v1*S1 + v2*S2 + v3*S3;
p = p1*S1 + p2*S2 + p3*S3;
r = r1*S1 + r2*S2 + r3*S3;

if CFL~=0
    l1 = sqrt((x1-x2)^2+(y1-y2)^2);
    l2 = sqrt((x3-x2)^2+(y3-y2)^2);
    l3 = sqrt((x1-x3)^2+(y1-y3)^2);

    L = min([l1 l2 l3]);

    V = sqrt(u^2+v^2);
    a = 1;
    dt =max(CFL*L/sqrt(V^2+a^2),dt);
end

dx_dL1 = x1 - x3;
dx_dL2 = x2 - x3;
dy_dL1 = y1 - y3;
dy_dL2 = y2 - y3;

J = [dx_dL1  dy_dL1;
     dx_dL2  dy_dL2];

detJ = (x1 - x3)*(y2 - y3) - (x2 - x3)*(y1 - y3);

dS1 = J\[1; 0];  dS1dx = dS1(1);  dS1dy = dS1(2);
dS2 = J\[0; 1];  dS2dx = dS2(1);  dS2dy = dS2(2);
dS3 = J\[-1; -1];  dS3dx = dS3(1);  dS3dy = dS3(2);

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
f = [r; u; v; p];

SK = transpose(L)*L;
SK = SK*detJ;

SF = transpose(L)*f;
SF = SF*detJ;


end

 