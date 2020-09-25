function [F, dFdX, C, Ceq, dCdX, dCeqdX] = CallFEMflow(Xi,Adjoint)

global Mach Options
global X0 Sol0 Sol t0

Cl_t = Sol0.Cl;
Cd0 = Sol0.Cd;
Cm0 = Sol0.Cm;

Xi = Xi.*X0;

Alpha = Xi(end);
X = Xi(1:end-1);

[Cl, Cd, Cm, ~, Sol, dCl_dAlpha, dCd_dAlpha, dCm_dAlpha, ~, ~, ~, dCl_dX, dCd_dX, dCm_dX] = FEMflow(Mach, Alpha, Adjoint, Options, X,Sol);
[t, dt_dX] = AirfoilThickness(X);

F = Cd/Cd0;
c1 = Cl/Cl_t - 1;
c2 = Cm/Cm0 - 1;
c3 = t0 - t;
C = [c2 c3];
Ceq = c1;

if Adjoint ==1
        
    dFdX = [dCd_dX/Cd0' dCd_dAlpha/Cd0].*X0;

    dc1_dX = dCl_dX/Cl_t;
    dc1_dAlpha = dCl_dAlpha/Cl_t;

    dc2_dX = dCm_dX/Cm0;
    dc2_dAlpha = dCm_dAlpha/Cm0;

    dc3_dX = -dt_dX';
    dc3_dAlpha = zeros(size(c3,2),1);

    dCdX =[dc2_dX     dc2_dAlpha;
           dc3_dX     dc3_dAlpha].*[X0; X0; X0];
             
    dCeqdX = [dc1_dX dc1_dAlpha].*X0; 
           
else
    dFdX = [];
    dCdX = [];
    dCeqdX = [];
end
