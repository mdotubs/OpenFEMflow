function [t, dt] = AirfoilThickness(A)

global xt
Au = A(1:length(A)/2);
Al = A(length(A)/2+1:end);

% Nx = 100;
% x=[1 1-sin(pi/2/Nx:pi/2/Nx:pi/2)];


%x = [0.378 0.38];
%x = [0.30 0.301];

[yu, ~, dydAu] =CSTairfoil(Au,xt);
[yl, ~, dydAl] =CSTairfoil(Al,xt);

t = yu-yl;

dt = [dydAu; zeros(size(dydAu))] - [zeros(size(dydAl)); dydAl];

