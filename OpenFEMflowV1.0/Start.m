%% VERSION CONTROL --------------------------------------------------------
% Author    : Ali Elham
% Date      : 2017


close all
clc

global Mach Options
global X0
global Sol0 Sol t0 xt

%% Setting inputs

Mach = 0.75;

Alpha = 1;

Adjoint = 0;

Options.MESH = 'NACA0012';  % name of the mesh fine without .msh! the mesh file should be placed in the MEsh folder
 
Options.Parallel = 0;         % 0 for no parallel computing; 1 for parallel computing

Options.rho = 1;              % scaled density
Options.P = 0.7143;           % scaled pressure pressure  

Options.dt  = 0.1;
Options.CFL = 1.0;

Options.iter_max = 400;                 % maximum number of iterations
Options.Newton_Switch_time = 2.0;       % Swith from Picard to Newton iteration after this (computational) time
Options.Newton_Switch_ratio = 0.5;      % Switch from Picard to Newton iteration when the residulas are reduced less than this ratio
Options. tol = 1e-9;

Options.XM = 0.25;  % X position for Cm calculation
Options.YM = 0;     % Y position for Cm calculation

Options.nMAD = 0;   % number of Mesh adaptation
Options.ep = 40;    % percentage of elements to be adapted
Options.tMAD = 2;   % typr of mesh adaptation 2 for dividing each element into two elements 3 for dividing into three elements


Options.MeshDef = 1;   % 0 for no mesh deformation, 1 for mesh deformationg
Options.Display = 1;   % 0 for not displaying anything; 1 for displaying the iteration and results 
Options.Plot = 1;      % 0 for no plot; 1 for plotting the mesh, convergence and results

%% Airfoil geometry using CST

 X0 = [0.1295    0.1370    0.1938    0.1903    0.2053   -0.1306   -0.1321   -0.2297   -0.0722    0.0388]; % RAE2822
% X0 = [0.1309    0.1238    0.1663    0.1387    0.1702    0.1989    0.1864  0.1943    0.1906    0.2155   -0.1310   -0.1300   -0.1600   -0.1185 -0.2359   -0.1279   -0.0730   -0.1165 -0.0078    0.0565]; % RAE2822
% X0 = [0.1718    0.1495    0.1637    0.1191    0.1690   -0.1718   -0.1495   -0.1637   -0.1191   -0.1690]; % NACA0012
% X0 =  [0.1743    0.1488    0.1914    0.1061    0.1962    0.1393    0.0969    0.2261    0.0642    0.2153  -0.1743    -0.1488    -0.1914    -0.1061    -0.1962    -0.1393    -0.0969    -0.2261    -0.0642    -0.2153]; % NACA0012
% X0 = [0.1628    0.2025    0.2612    0.2412    0.2411   -0.1023   -0.0539   -0.1752    0.0184  -0.0383]; % NACA 64A412

%% Running FEMflow

%[Cl, Cd, Cm, CP, Sol0, dCl_dAlpha, dCd_dAlpha, dCm_dAlpha, dCl_dMach, dCd_dMach, dCm_dMach, dCl_dX, dCd_dX, dCm_dX] = FEMflow(Mach, Alpha, Adjoint, Options, X0);
[Cl, Cd, Cm, CP, Sol0] = FEMflow(Mach, Alpha, Adjoint, Options, X0);  % add Sol as the last input for FEMflow, in case you would like to start from an existing solution

Sol = Sol0;

xt = [0.25 0.65];
t0 = AirfoilThickness(X0);
 
%% Sensitivity analysis check

% dx = 1e-6;
% for n=1:20
%     Xf = X0;
%     Xf(n) = Xf(n) + dx;
%     [Clf, Cdf, Cmf] = FEMflow(Mach, Re, Alpha, 0, Options, Xf);
%     
%     dCl_dX_f(n) = (Clf-Cl)/dx;
%     dCd_dX_f(n) = (Cdf-Cd)/dx;
%     dCm_dX_f(n) = (Cmf-Cm)/dx;
% end

%% Optimization using build in SQP algorithm

% Options.Plot = 0;
% 
% Xu = [3*X0(1:5) 1/3*X0(6:10)];
% Xl = [1/3*X0(1:5) 3*X0(6:10)];
% 
% iter_max_sqp = 50;
% 
% tic
% [Xopt, Fopt, X_iter] = SQP(X0,Alpha,iter_max_sqp);
% t_sqp = toc;
% 
% Options.Plot = 1;
% [Clopt, Cdopt, Cmopt] = FEMflow(Mach, Xopt(11), 0, Options, Xopt(1:10),Sol);

%% Optimization using MATLAB FMINCON

Sol = Sol0;

Options.Plot = 0;
X0 = [X0 1];  % adding Alpha to DV
 
x0 = [ones(1,length(X0)-1) Alpha];

ub = [3*ones(1,length(X0)-1) 5];
lb = [1/3*ones(1,length(X0)-1) -2];

plots = {@optimplotx @optimplotfval @optimplotfirstorderopt @optimplotconstrviolation}; 
options = optimset('GradObj','on','GradConstr','on','algorithm','sqp','display','iter','TolFun',1e-4,'TolCon',1e-4,'MaxIter',50,'MaxFunEvals',100,'PlotFcns',plots);

tic
[Xopt, Fopt] = RunOptimization(x0, ub, lb,options);
Xopt = Xopt.*X0;
t_fmincon = toc;
save('Sol_fmincon.mat','Sol');

Options.Plot = 1;
[Clopt, Cdopt, Cmopt] = FEMflow(Mach, Xopt(end), 0, Options, Xopt(1:end-1),Sol);

%% Optimization using SNOPT

% Options.Plot = 0;
% 
% Sol = Sol0;
% X0 = [X0'; 1];  % adding Alpha to DV
% 
% snscreen on
% 
% snprint('MDO.out'); 
% MDO.spc = which('MDO.spc');
% snspec(MDO.spc);
% 
% flow = [0; 0; 0;  0;   0];
% fupp = [1; 0; inf; inf; inf];
% 
% x0 = [ones(length(X0)-1,1); Alpha];
% 
% ub = [1.5*ones(length(X0)-1,1); 5];
% lb = [2/3*ones(length(X0)-1,1); -2];
% 
% [iGfun,jGvar] = find(ones(length(flow),length(X0)));
% iAfun = [];  jAvar = [];  A     = [];
% 
% tic
% [x,f,inform,xmul,Fmul] = snopt(x0,lb, ub,flow,fupp,'Objfun',A, iAfun, jAvar, iGfun, jGvar);
% t_snopt = toc;
% save('Sol_snopt.mat','Sol');
