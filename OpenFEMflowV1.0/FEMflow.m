function [Cl, Cd, Cm, CP, Sol, dCl_dAlpha, dCd_dAlpha, dCm_dAlpha, dCl_dMach, dCd_dMach, dCm_dMach, dCl_dX, dCd_dX, dCm_dX] = FEMflow(Mach, Alpha, Adjoint_Switch, Options, X, Sol)

global dt CFL 
global nodes element nDOF BC theta Surf
global P Parallel
global XM YM
global Adjoint

%% settings

Parallel = Options.Parallel;
rho = Options.rho;
P = Options.P;
dt = Options.dt;
CFL = Options.CFL;


XM = Options.XM;
YM = Options.YM;

nMAD = Options.nMAD+1;
ep = Options.ep;


MESH = Options.MESH;

%% flow

Uinf = Mach*cosd(Alpha);  % scaled velocity
Vinf = Mach*sind(Alpha);  

    
%% Settings

Adjoint = 0;
    
iter_max = Options.iter_max;
Newton_Switch_time = Options.Newton_Switch_time;
Newton_Switch_ratio = Options.Newton_Switch_ratio;
tol = Options.tol;
    
%% Mesh    
addpath('./Mesh')

if Options.Display ==1
display('Reading mesh ...');
end

if exist('Sol','var')
    nodes = Sol.nodes;
    element = Sol.element;
    BC = Sol.BC;
    theta = Sol.theta;
    Surf = Sol.Surf;
    U = Sol.U;
else
    [nodes, element, BC, theta, Surf] = ReadMesh([MESH '.msh']);
end

nDOF = size(nodes,1)*4;

if Options.Display ==1
display(['Number of nodes: ' num2str(size(nodes,1)) ', Number of elements: ' num2str(size(element,1))]);
end

%% Mesh Deformation

if Options.MeshDef ==1
    if Options.Display ==1
        display('Deforming mesh ...');
    end
    
    [nodes, theta, ~, dYnodes_dX, dtheta_dX] = MeshDef(nodes,element,X);

end

if Options.Plot ==1
   MeshPlot(nodes,element);
end

%% Initializing and applying BC

if Options.Display ==1
display('Setting boundary conditions ...');
end

if ~exist('U','var')
    U = zeros(nDOF,1);
    U(1:4:end) = rho;
    U(2:4:end) = Uinf;
    U(3:4:end) = Vinf;
    U(4:4:end) = P;
end

for mad= 1: nMAD
    if mad~=1 
        Uxy = U;
        Uxy(Walldof) = T*U(Walldof);
        if Options.Display ==1
        display('Automatic Mesh Refinement ...');
        end
        
        if Options.tMAD ==2
            [nodes, element, U] = MeshAdaptation2(nodes,element,Uxy,ep);
        elseif Options.tMAD ==3
            [nodes, element, U] = MeshAdaptation3(nodes,element,Uxy,ep);
        end
        
        if Options.MeshDef ==1
             [nodes, theta, ~, dYnodes_dX, dtheta_dX] = MeshDef(nodes,element,X);
        end
         if Options.Display ==1
             display(['Number of nodes: ' num2str(size(nodes,1)) ', Number of elements: ' num2str(size(element,1))]);
         end
        
%         if mad==2
%            Newton_Switch_time = Newton_Switch_time/2;
%         end
    end

   
    nDOF = size(nodes,1)*4;
    FixedDOF = zeros(1,length(BC(BC(:,2)==1))+length(BC(BC(:,2)==2))*4+length(BC(BC(:,2)==3)));
    iDOF = 1;
    Iu = zeros(1,length(BC(BC(:,2)==2)));
    Iv = zeros(1,length(BC(BC(:,2)==2)));
    iI = 1;

    Walldof = [];
    T = sparse(nDOF,nDOF); 

    for i=1:size(BC,1)
        if BC(i,2) ==1 % WALL: Un = 0

%           changing [u;v] to [ut; un]
%           U(4*(BC(i,1)-1)+2) = U(4*(BC(i,1)-1)+2)*cos(theta(theta(:,1)==BC(i,1),2))+U(4*(BC(i,1)-1)+3)*sin(theta(theta(:,1)==BC(i,1),2));    % U --> Ut

            U(4*(BC(i,1)-1)+3) = 0;   % V --> Un = 0   

            T(4*(BC(i,1)-1)+2:4*(BC(i,1)-1)+3,4*(BC(i,1)-1)+2:4*(BC(i,1)-1)+3) = [cos(theta(theta(:,1)==BC(i,1),2)) -sin(theta(theta(:,1)==BC(i,1),2)); ...
                                                                                  sin(theta(theta(:,1)==BC(i,1),2))  cos(theta(theta(:,1)==BC(i,1),2))];

            FixedDOF(iDOF) = 4*(BC(i,1)-1)+3;  % Un is fixed
            iDOF = iDOF + 1;

            Walldof(end+1:end+2) = [4*(BC(i,1)-1)+2; 4*(BC(i,1)-1)+3];


        elseif BC(i,2) ==2 % Inflow

            U(4*(BC(i,1)-1)+1) = rho;
            U(4*(BC(i,1)-1)+2) = Uinf;
            U(4*(BC(i,1)-1)+3) = Vinf;
            U(4*(BC(i,1)-1)+4) = P;
            FixedDOF(iDOF) = 4*(BC(i,1)-1)+1;
            FixedDOF(iDOF+1) = 4*(BC(i,1)-1)+2;
            FixedDOF(iDOF+2) = 4*(BC(i,1)-1)+3;
            FixedDOF(iDOF+3) = 4*(BC(i,1)-1)+4;
            Iu(iI) = 4*(BC(i,1)-1)+2;  % inlet U velocity DOF
            Iv(iI) = 4*(BC(i,1)-1)+3;  % inlet V celocity DOF              

            iDOF = iDOF + 4;
            iI  = iI + 1;

        elseif BC(i,2) ==3 % Outflow

            U(4*(BC(i,1)-1)+4) = P;
            FixedDOF(iDOF) = 4*(BC(i,1)-1)+4;   
            iDOF = iDOF + 1;

        end
    end

    Walldof = Walldof';
    T = T(Walldof,Walldof);

    ActiveDOF = setdiff(1:nDOF,FixedDOF);
    
           
    %% Solving using Picard method
    if Options.Display ==1
    display('Solving the prolem ...');
    end
    
    dU = zeros(size(U));
    Method = 'Picard';
    Switch_time = Newton_Switch_time;
    time = 0;
    Res = [];
    Error = [];

       
    for iter=1: iter_max 
              
       %tic    
       [R, dRdU] = Matrices_Parallel(U,Method);
       
       dU(ActiveDOF) = -dRdU(ActiveDOF,ActiveDOF)\R(ActiveDOF);   
       U = U + dU; 

       Error(iter)  = norm(dU);
       Res(iter) = norm(R(ActiveDOF));
       
       if iter==1
           Res_Picard = Res(1);
       end
             
       [Cl, Cd, Cm] = Forces(U, Alpha, Mach);
            
                
        time = time + dt;
               
        if Options.Display ==1
            if iter == 1 || round(iter/20) == iter/20
                fprintf('\n');
                fprintf('%s     %s         %s           %s           %s        %s\n', 'Iter', 'Residual', 'Cl','Cd','Cm ','Method');
            end      
            fprintf('%d .      %f .   %f .   %f .   %f .    %s\n',iter, Res(end), Cl, Cd, Cm, Method);
        end
        
               
        if Res(end) <=tol
            break;
        elseif Res(end) > 10*Res(1)
           if strcmpi(Method,'Newton')
               if Options.Display ==1
                display('Newton iteration is diverging! Switching back to Picard Iteration.');   
               end
                U = Upicard;
                time = time_Picard;
                
                Method = 'Picard';
                Switch_time =  time+Newton_Switch_time;
           else
               if Options.Display ==1
                display('Solution is diverging! Stopping the iteration.');
               end
                break;
            end
        elseif strcmpi(Method,'Picard') && ( time >=Switch_time ||  Res(end)/Res_Picard <=Newton_Switch_ratio)
            Method = 'Newton';
            Upicard = U;
            time_Picard = time;
            Res_Picard = Res(end);
        end
    end
end


if Res(end) <=tol
    if Options.Display ==1
        display(['Convergence is achieved at ' num2str(iter) ' iterations with the norm of residual equal to ' num2str(Res(end))]);
        display(['Cl = ' num2str(Cl) '  Cd = ' num2str(Cd) '  Cm = ' num2str(Cm)]);
    end
        ExitFlag =1;
else
    if Options.Display ==1
        display('Convergence is NOT achieved. Increase the number of iterations.');
    end
        ExitFlag =0;
end

if Options.Plot ==1
    if nMAD>1
        MeshPlot(nodes,element);
    end

    figure
    hold on
    set(gca,'fontsize',14);
    plot(log10(Error),'linewidth',1.5);%-log10(Error(1)));
    plot(log10(Res),'linewidth',1.5); %-log10(Res(1)));
    grid
    xlabel('Iteration');
    ylabel('Log(R), Log(\Delta U)');
    legend('Norm of \Delta U','Norm of Residual');
    hold off
end
% save('Res2.mat','Res');


Uxy = U;
Uxy(Walldof) = T*U(Walldof); 

CP = PostProc(Uxy,Mach,Options);

Sol.nodes = nodes;
Sol.element = element;
Sol.U = U;
Sol.BC = BC;
Sol.theta = theta;
Sol.Surf = Surf;
Sol.CP = CP;
Sol.Alpha = Alpha;
Sol.Cl = Cl;
Sol.Cd = Cd;
Sol.Cm = Cm;
Sol.ExitFlag = ExitFlag;

%% Sensitivity analysis

if Adjoint_Switch ==1

    display('Adjoint sensitivity analysis ...')
    Adjoint = 1;
    [~, dRdU, ~, dRdYnodes, dRdtheta] = Matrices_Parallel(U,'Newton');
    [Cl, Cd, Cm, pCl_pAlpha, pCd_pAlpha, pCl_pMach, pCd_pMach, pCm_pMach, dCl_dU, dCd_dU, dCm_dU, ~, ~, ~, pCl_pYnodes, pCd_pYnodes, pCm_pYnodes] = Forces(U, Alpha, Mach);


    % **************** Adjoint vectors *************************
    Lambda_cl = dRdU(ActiveDOF,ActiveDOF)'\dCl_dU(ActiveDOF);
    Lambda_cd = dRdU(ActiveDOF,ActiveDOF)'\dCd_dU(ActiveDOF);
    Lambda_cm = dRdU(ActiveDOF,ActiveDOF)'\dCm_dU(ActiveDOF);


    % ****************** d_dAlpha **************************
    dU_dAlpha = sparse([Iu;Iv],ones(length(Iu)+length(Iv),1),[-Mach*sind(Alpha)*ones(size(Iu)); Mach*cosd(Alpha)*ones(size(Iv))],nDOF,1);
    dR_dAlpha = dRdU*dU_dAlpha;

    dCl_dAlpha = pCl_pAlpha - Lambda_cl'*dR_dAlpha(ActiveDOF);
    dCd_dAlpha = pCd_pAlpha - Lambda_cd'*dR_dAlpha(ActiveDOF);
    dCm_dAlpha = - Lambda_cm'*dR_dAlpha(ActiveDOF);
    dCl_dAlpha = dCl_dAlpha * pi/180;
    dCd_dAlpha = dCd_dAlpha * pi/180;
    dCm_dAlpha = dCm_dAlpha * pi/180;

    % ***************** d_dMach **************************
%     dU_dMach = sparse([Iu;Iv],ones(length(Iu)+length(Iv),1),[cosd(Alpha)*ones(size(Iu)); sind(Alpha)*ones(size(Iv))],nDOF,1);
%     dR_dMach = dRdU*dU_dMach;
% 
%     dCl_dMach = pCl_pMach - Lambda_cl'*dR_dMach(ActiveDOF);
%     dCd_dMach = pCd_pMach - Lambda_cd'*dR_dMach(ActiveDOF);
%     dCm_dMach = pCm_pMach - Lambda_cm'*dR_dMach(ActiveDOF);

%   d_dMach is not needed for optimization / if you need it uncomment the lines above
    dCl_dMach = [];
    dCd_dMach = [];
    dCm_dMach = [];

    % **************** d_dX ******************************
    dCl_dYnodes = pCl_pYnodes - Lambda_cl'*dRdYnodes(ActiveDOF,:);
    dCd_dYnodes = pCd_pYnodes - Lambda_cd'*dRdYnodes(ActiveDOF,:);
    dCm_dYnodes = pCm_pYnodes - Lambda_cm'*dRdYnodes(ActiveDOF,:);

    dCl_dtheta = - Lambda_cl'*dRdtheta(ActiveDOF,:);
    dCd_dtheta = - Lambda_cd'*dRdtheta(ActiveDOF,:);
    dCm_dtheta = - Lambda_cm'*dRdtheta(ActiveDOF,:);

    dCl_dX = dCl_dYnodes*dYnodes_dX + dCl_dtheta*dtheta_dX;
    dCd_dX = dCd_dYnodes*dYnodes_dX + dCd_dtheta*dtheta_dX;
    dCm_dX = dCm_dYnodes*dYnodes_dX + dCm_dtheta*dtheta_dX;

else
   dCl_dAlpha = [];
   dCd_dAlpha = [];
   dCm_dAlpha = [];
   dCl_dMach = [];
   dCd_dMach = [];
   dCm_dMach = [];
   dCl_dX = [];
   dCd_dX = [];
   dCm_dX = [];
end
