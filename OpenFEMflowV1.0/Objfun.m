function [F, G] = Objfun(Xi)

global Mach Re Options
global X0 Sol0 Sol t0

Cl_t = Sol0.Cl;
Cd0 = Sol0.Cd;
Cm0 = Sol0.Cm;

Xi = Xi.*X0;

Alpha = Xi(end);
X = Xi(1:end-1);

[Cl, Cd, Cm, ~, Sol, dCl_dAlpha, dCd_dAlpha, dCm_dAlpha, ~, ~, ~, dCl_dX, dCd_dX, dCm_dX] = FEMflow(Mach, Alpha, 1, Options, X',Sol);
[t, dt_dX] = AirfoilThickness(X);

c1 = Cl/Cl_t - 1;
dc1_dX = dCl_dX/Cl_t;
dc1_dAlpha = dCl_dAlpha/Cl_t;

c2 = 1 - Cm/Cm0;
dc2_dX = -dCm_dX/Cm0;
dc2_dAlpha = -dCm_dAlpha/Cm0;

c3 = t - t0;
dc3_dX = dt_dX';
dc3_dAlpha = zeros(size(c3,2),1);

F = [Cd/Cd0; c1; c2; c3'];
G = [dCd_dX/Cd0 dCd_dAlpha/Cd0;
     dc1_dX     dc1_dAlpha; 
     dc2_dX     dc2_dAlpha;
     dc3_dX     dc3_dAlpha].*[X0'; X0'; X0'; X0'; X0'];


